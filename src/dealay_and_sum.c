#include "../include/setting.h"
#include "../include/fft_path.h"
#include "../include/gcc_phat_delay.h"
#include "../include/tf_gsc.h"
#include "../include/mmse_lsa.h"
#include "../include/signal_restore.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// 计算互相关，估计第二个信号相对于第一个信号的延迟（样本数）
// 返回值：延迟 d，满足 y[n] ≈ x[n-d] (d>0表示y滞后于x)
int estimate_delay(const int16_t* x, uint32_t len_x,
                   const int16_t* y, uint32_t len_y,
                   int max_delay) {
    // 限制最大搜索范围不超过信号长度
    int max_lag = max_delay;
    if (max_lag > (int)len_x - 1) max_lag = len_x - 1;
    if (max_lag > (int)len_y - 1) max_lag = len_y - 1;
    if (max_lag < 0) max_lag = 0;

    int best_delay = 0;
    double max_corr = -1e9;

    // 时域互相关：对于每个可能的延迟 k，计算 x[0..N-1] 与 y[k..k+N-1] 的点积
    // 为了简单，只使用两个信号重叠部分的最大长度
    // 实际使用重叠长度取决于延迟，但为了效率，我们取固定长度，忽略边界效应
    // 更精确的做法是对于每个 k，取有效重叠部分
    // 这里简化：使用所有重叠样本（动态计算）
    for (int k = -max_lag; k <= max_lag; k++) {
        double corr = 0.0;
        int start_x = (k < 0) ? -k : 0;
        int start_y = (k > 0) ? k : 0;
        int overlap_len = len_x - start_x;
        if (len_y - start_y < overlap_len) overlap_len = len_y - start_y;
        if (overlap_len <= 0) continue;

        for (int i = 0; i < overlap_len; i++) {
            corr += (double)x[start_x + i] * y[start_y + i];
        }
        if (corr > max_corr) {
            max_corr = corr;
            best_delay = k;
        }
    }
    return best_delay;
}




// 延迟求和，两个信号各自延迟 delay1 和 delay2 样本（支持子样点精度），求和并平均
int16_t* delay_sum(const int16_t* data1, uint32_t len1, float delay1,
                   const int16_t* data2, uint32_t len2, float delay2,
                   uint32_t* outLen)
{
    // 找到最小延迟，将所有延迟调整为非负
    float minDelay = (delay1 < delay2) ? delay1 : delay2;
    float adjDelay1 = delay1 - minDelay;
    float adjDelay2 = delay2 - minDelay;

    // 输出长度 = max(ceil(adjDelay) + len, ...)
    uint32_t outLen1 = len1 + (uint32_t)ceilf(adjDelay1);
    uint32_t outLen2 = len2 + (uint32_t)ceilf(adjDelay2);
    *outLen = (outLen1 > outLen2) ? outLen1 : outLen2;

    int16_t* outData = (int16_t*)calloc(*outLen, sizeof(int16_t));
    if (!outData) {
        fprintf(stderr, "Memory allocation failed\n");
        return NULL;
    }

    // 求和（线性插值处理子样点延迟）
    for (uint32_t i = 0; i < *outLen; i++)
    {
        float sum = 0.0f;

        // 通道1
        float pos1 = (float)i - adjDelay1;
        if (pos1 >= 0.0f && pos1 < (float)len1 - 1.0f) {
            int idx = (int)pos1;
            float frac = pos1 - (float)idx;
            sum += (1.0f - frac) * (float)data1[idx] + frac * (float)data1[idx + 1];
        } else if (pos1 >= 0.0f && pos1 < (float)len1) {
            sum += (float)data1[(int)pos1];
        }

        // 通道2
        float pos2 = (float)i - adjDelay2;
        if (pos2 >= 0.0f && pos2 < (float)len2 - 1.0f) {
            int idx = (int)pos2;
            float frac = pos2 - (float)idx;
            sum += (1.0f - frac) * (float)data2[idx] + frac * (float)data2[idx + 1];
        } else if (pos2 >= 0.0f && pos2 < (float)len2) {
            sum += (float)data2[(int)pos2];
        }

        int16_t val = (int16_t)(sum / 2.0f);
        outData[i] = val;
    }

    return outData;
}

/**
 * @brief 频域延迟补偿+波束成形（帧循环+重叠相加）
 *
 * 对输入全长音频做分帧处理：
 *   每帧 Hamming 加窗 → FFT → 相位旋转(补偿浮点延迟)
 *   → DSB 频域相加 → 300~3400Hz 带通掩码 → IFFT
 *   → 重叠相加合成完整时域输出
 */
int16_t* freq_domain_beamform(GccPhatContext* ctx,
                              const int16_t* data1, uint32_t len1,
                              const int16_t* data2, uint32_t len2,
                              float delay, int fs, uint32_t* outLen)
{
    int fft_size  = ctx->fft_size;      // 1024
    int input_len = ctx->input_len;     // 2048 = 2*fft_size
    int nbins     = ctx->complex_len;   // fft_size + 1
    int hop       = fft_size / 2;       // 50% 重叠

    // 人声带通边界 + 过渡段
    float bin_hz = (float)fs / (float)input_len;
    int bin_low  = (int)ceilf(300.0f / bin_hz);
    int bin_high = (int)floorf(3400.0f / bin_hz);

    uint32_t min_len = (len1 < len2) ? len1 : len2;
    if (min_len == 0) { *outLen = 0; return NULL; }

    // 帧数计算
    int nframes = (min_len + hop - 1) / hop;
    if (nframes < 1) nframes = 1;
    uint32_t total_out = (nframes - 1) * hop + fft_size;

    // 浮点累加缓冲 + 归一化权重累加 + VAD 权重累加
    float* out_float = (float*)calloc(total_out, sizeof(float));
    float* out_norm  = (float*)calloc(total_out, sizeof(float));
    float* out_vad   = (float*)calloc(total_out, sizeof(float));
    if (!out_float || !out_norm || !out_vad) {
        free(out_float); free(out_norm); free(out_vad);
        *outLen = 0; return NULL;
    }

    // 自适应帧级 VAD 阈值跟踪
    double max_frame_energy_db = -100.0;  // 对数域平滑最大值

    // 初始化 TF-GSC（频域自适应噪声对消）
    // mu=0.025 单语音 NLMS 步长, alpha=0.92 功率平滑, beta=0.4 频点VAD门限
    // leak=0.9995 语音帧缓慢衰减 W，减少人声泄漏抵消
    // smooth_factor=0.7 BM泄漏维纳增益平滑（快速跟踪）
    TfGscContext* tf_gsc = tf_gsc_init(nbins, 0.025f, 0.92f, 0.4f, 0.9995f, 0.70f);
    if (!tf_gsc) {
        free(out_float); free(out_norm);
        *outLen = 0;
        return NULL;
    }

    // 配置 W 硬幅值钳位 ±0.15，防止权重过大吞噬人声
    tf_gsc->W_max = 0.15f;

    // 配置 BM 泄漏抑制：语音基频 80-800Hz G_U ≥ 0.9
    tf_gsc->gU_low_bin  = (int)ceilf(80.0f / bin_hz);
    tf_gsc->gU_high_bin = (int)floorf(800.0f / bin_hz);
    tf_gsc->gU_min      = 0.9f;

    // 模块 A: 时频扩散度掩膜 — 默认启用（统计量始终平滑，掩膜始终生效）
    // 模块 B: FBF 前置直达语音增强 — 开启，从源头抬高直达人声基底
    tf_gsc->gfb_enabled = 1;

    // 模块 C: ANC 输出人声下限钳位 — 开启，保留最低 35% FBF 人声
    tf_gsc->clamp_enabled   = 1;
    tf_gsc->clamp_min_ratio = 0.35f;

    // 初始化后置 MMSE-LSA（降 TF-GSC 残留的不相关底噪）
    MmseLsaCtx lsa;
    lsa_init(&lsa, nbins);

    // 初始化信号恢复模块（补偿波束成形对期望语音的损伤）
    // lambda=0.88 余弦相似度平滑, eta=0.10/0.03 三段式补偿,
    // beta_low=10(低SNR冻结), beta_high=5(高SNR适度增强), SNR门限=15dB
    SignalRestoreContext *restore = restore_init(nbins, 0.88f, 0.10f, 0.03f, 10.0f, 5.0f, 15.0f);

    // ── MSC 时频掩膜状态（双麦幅度平方相干） ──
    const float msc_alpha = 0.85f;                    // 功率谱时间平滑因子
    float* msc_S11 = (float*)calloc(nbins, sizeof(float));  // 平滑 |X1|²
    float* msc_S22 = (float*)calloc(nbins, sizeof(float));  // 平滑 |X2|²
    float* msc_S12_re = (float*)calloc(nbins, sizeof(float)); // 平滑 Re(X1·X2*)
    float* msc_S12_im = (float*)calloc(nbins, sizeof(float)); // 平滑 Im(X1·X2*)
    float* msc_gain  = (float*)malloc(nbins * sizeof(float));
    if (!msc_S11 || !msc_S22 || !msc_S12_re || !msc_S12_im || !msc_gain) {
        free(msc_S11); free(msc_S22); free(msc_S12_re); free(msc_S12_im); free(msc_gain);
        free(out_float); free(out_norm); free(out_vad);
        *outLen = 0; return NULL;
    }
    for (int k = 0; k < nbins; k++) msc_gain[k] = 1.0f;  // 初始=1

    for (int f = 0; f < nframes; f++)
    {
        int start = f * hop;
        int remain = (int)min_len - start;
        int frame_len = (remain < fft_size) ? remain : fft_size;
        if (frame_len <= 0) break;

        // 清空输入缓冲
        memset(ctx->in1, 0, input_len * sizeof(float));
        memset(ctx->in2, 0, input_len * sizeof(float));

        // Hamming 加窗
        for (int i = 0; i < frame_len; i++) {
            float w = ctx->window_coeffs[i];
            ctx->in1[i] = (float)data1[start + i] * w;
            ctx->in2[i] = (float)data2[start + i] * w;
        }

        // FFT
        fftwf_execute(ctx->plan_fwd1);
        fftwf_execute(ctx->plan_fwd2);

        // 相位旋转补偿亚采样延迟（复制 mic2 频谱用于旋转）
        fftwf_complex mic2_align[nbins];
        memcpy(mic2_align, ctx->out2, sizeof(fftwf_complex) * nbins);
        for (int k = 0; k < nbins; k++) {
            float omega_k = 2.0f * (float)3.14159265358979323846f * k / (float)input_len;
            float phase = omega_k * delay;
            float c = cosf(phase), s = sinf(phase);
            float yr = mic2_align[k][0], yi = mic2_align[k][1];
            mic2_align[k][0] = yr * c + yi * s;
            mic2_align[k][1] = -yr * s + yi * c;
        }

        // TF-GSC 频域自适应噪声对消
        //   同时输出: Y_FBF = S (FBF 主波束), Y_U = Ñ (净化噪声参考)
        fftwf_complex y_fbf[nbins], y_u[nbins];
        tf_gsc_process_frame(tf_gsc, ctx->out1, mic2_align, y_fbf, y_u, ctx->cs, nbins);

        // 信号恢复模块：根据余弦相似度选择麦克信号，补偿 FBF 对期望语音的损伤
        //   Y_bar_FBF = Y_FBF + G_X · X_selected
        //   最终输出叠加: Y_GSC + (Y_bar_FBF - Y_FBF) = Y_GSC + G_X · X_selected
        if (restore) {
            fftwf_complex y_restored[nbins];
            restore_process_frame(restore, y_fbf, y_u, ctx->out1, mic2_align, y_restored, nbins);
            for (int k = 0; k < nbins; k++) {
                ctx->cs[k][0] += y_restored[k][0] - y_fbf[k][0];
                ctx->cs[k][1] += y_restored[k][1] - y_fbf[k][1];
            }
        }

        // ── MSC 双麦相干后滤波（抑制全向扩散噪声） ──
        //  MSC = |E{X1·X2*}|² / (E{|X1|²}·E{|X2|²})
        //  MSC>0.8→G=1, MSC<0.5→G=0.25(-12dB), 线性渐变
        for (int k = 0; k < nbins; k++) {
            float re_x1 = ctx->out1[k][0], im_x1 = ctx->out1[k][1];
            float re_x2 = mic2_align[k][0], im_x2 = mic2_align[k][1];
            float p1 = re_x1*re_x1 + im_x1*im_x1;
            float p2 = re_x2*re_x2 + im_x2*im_x2;
            float cr = re_x1*re_x2 + im_x1*im_x2;   // Re(X1·X2*)
            float ci = im_x1*re_x2 - re_x1*im_x2;   // Im(X1·X2*)

            msc_S11[k] = msc_alpha * msc_S11[k] + (1.0f - msc_alpha) * p1;
            msc_S22[k] = msc_alpha * msc_S22[k] + (1.0f - msc_alpha) * p2;
            msc_S12_re[k] = msc_alpha * msc_S12_re[k] + (1.0f - msc_alpha) * cr;
            msc_S12_im[k] = msc_alpha * msc_S12_im[k] + (1.0f - msc_alpha) * ci;

            float S11 = msc_S11[k], S22 = msc_S22[k];
            float S12_mag2 = msc_S12_re[k]*msc_S12_re[k] + msc_S12_im[k]*msc_S12_im[k];
            float denom = S11 * S22 + 1e-10f;
            float msc = S12_mag2 / denom;
            if (msc > 1.0f) msc = 1.0f;

            // MSC→增益映射：0.5以下衰减12dB，0.8以上保留
            float g_msc;
            if (msc < 0.5f)
                g_msc = 0.25f;
            else if (msc < 0.8f)
                g_msc = 0.25f + (msc - 0.5f) / 0.3f * 0.75f;
            else
                g_msc = 1.0f;

            // 二阶平滑，防止音乐噪声（0.5递归平滑）
            msc_gain[k] = 0.5f * msc_gain[k] + 0.5f * g_msc;
            ctx->cs[k][0] *= msc_gain[k];
            ctx->cs[k][1] *= msc_gain[k];
        }

        // 后置 MMSE-LSA 降噪：消除 TF-GSC 无法处理的不相关底噪
        lsa_process_frame(&lsa, ctx->cs);

        // ── 帧级 VAD 检测（基于 GSC 输出能量） ──
        double frame_energy_db = -100.0;
        {
            double e = 0.0;
            for (int k = 0; k < nbins; k++)
                e += (double)(ctx->cs[k][0]*ctx->cs[k][0] + ctx->cs[k][1]*ctx->cs[k][1]);
            e /= (double)nbins;
            frame_energy_db = 10.0 * log10(e + 1e-30);
        }

        // 跟踪全局能量峰值，带 6dB/s 慢衰减（防止静音后误判）
        if (frame_energy_db > max_frame_energy_db)
            max_frame_energy_db = frame_energy_db;
        else
            max_frame_energy_db -= 0.15;  // 约 6dB/s @ 50fps

        int vad_active = (frame_energy_db > max_frame_energy_db - 25.0) ? 1 : 0;

        // VAD 门控带通掩膜：VAD=1 时压制带外噪声，VAD=0 时全通
        if (vad_active) {
            // 300~3400Hz 带通掩码（仅语音帧启用）
            for (int k = 0; k < nbins; k++) {
                if (k < bin_low || k > bin_high) {
                    ctx->cs[k][0] = 0.0f;
                    ctx->cs[k][1] = 0.0f;
                }
            }
        }

        // IFFT
        fftwf_execute(ctx->plan_inv);

        // 重叠相加 + VAD 权重累积
        for (int i = 0; i < fft_size; i++) {
            int pos = start + i;
            if (pos < (int)total_out) {
                out_float[pos] += ctx->xcorr[i];
                out_norm[pos]  += ctx->window_coeffs[i % fft_size];
                if (vad_active)
                    out_vad[pos] += ctx->window_coeffs[i % fft_size];
            }
        }
    }

    // 销毁 TF-GSC 上下文
    tf_gsc_destroy(tf_gsc);

    // 归一化：xcorr = input_len * 2ch * signal * window_overlap_sum
    // 所以 signal = xcorr / (input_len * 2 * norm_sum)
    int16_t* out = (int16_t*)malloc(total_out * sizeof(int16_t));
    if (!out) { free(out_float); free(out_norm); *outLen = 0; return NULL; }

    /* ========== 输出时域滤波器链 ========== */

    // 1. 六阶 Butterworth 高通 fc=80Hz（3 级级联双二阶），提升低频衰减斜率
    //    原二阶 @80Hz 仅 12dB/oct，六阶达 36dB/oct，有效压制 50~80Hz 共振
    //    系数由模拟极点经双线性变换预计算，fs=16kHz
    //    第1级 (σ=0.2588)
    float hp1_b0 =  0.9917f, hp1_b1 = -1.9834f, hp1_b2 =  0.9917f;
    float hp1_a1 = -1.9829f, hp1_a2 =  0.9839f;
    float hp1_x1=0, hp1_x2=0, hp1_y1=0, hp1_y2=0;
    //    第2级 (σ=0.7071)
    float hp2_b0 =  0.9780f, hp2_b1 = -1.9561f, hp2_b2 =  0.9780f;
    float hp2_a1 = -1.9556f, hp2_a2 =  0.9566f;
    float hp2_x1=0, hp2_x2=0, hp2_y1=0, hp2_y2=0;
    //    第3级 (σ=0.9659)
    float hp3_b0 =  0.9703f, hp3_b1 = -1.9406f, hp3_b2 =  0.9703f;
    float hp3_a1 = -1.9401f, hp3_a2 =  0.9411f;
    float hp3_x1=0, hp3_x2=0, hp3_y1=0, hp3_y2=0;

    // 2. 自适应 50Hz 窄带陷波（Q=15），消除工频谐波
    float notch_Q = 15.0f;
    float w0_50 = 2.0f * (float)M_PI * 50.0f / (float)fs;
    float alpha_50 = sinf(w0_50) / (2.0f * notch_Q);
    float a0_50 = 1.0f + alpha_50;
    float n50_b0 = 1.0f / a0_50;
    float n50_b1 = -2.0f * cosf(w0_50) / a0_50;
    float n50_b2 = 1.0f / a0_50;
    float n50_a1 = n50_b1;    /* b1/a0 == -2cos(w0)/a0 */
    float n50_a2 = (1.0f - alpha_50) / a0_50;
    float n50_x1 = 0.0f, n50_x2 = 0.0f, n50_y1 = 0.0f, n50_y2 = 0.0f;

    // 3. 自适应 100Hz 窄带陷波（Q=15），消除工频二次谐波
    float w0_100 = 2.0f * (float)M_PI * 100.0f / (float)fs;
    float alpha_100 = sinf(w0_100) / (2.0f * notch_Q);
    float a0_100 = 1.0f + alpha_100;
    float n100_b0 = 1.0f / a0_100;
    float n100_b1 = -2.0f * cosf(w0_100) / a0_100;
    float n100_b2 = 1.0f / a0_100;
    float n100_a1 = n100_b1;
    float n100_a2 = (1.0f - alpha_100) / a0_100;
    float n100_x1 = 0.0f, n100_x2 = 0.0f, n100_y1 = 0.0f, n100_y2 = 0.0f;

    // 4. 高频搁架（High-shelf, fc=1.2kHz, +6dB），还原齿音/辅音清晰度
    float hs_gain_db = 6.0f;
    float A_hs = powf(10.0f, hs_gain_db / 40.0f);  /* A = 10^(dB/40) */
    float w0_hs = 2.0f * (float)M_PI * 1200.0f / (float)fs;
    float alpha_hs = sinf(w0_hs) * 0.70710678f;    /* S=1 → alpha = sin(w0)/√2 */
    float cos_hs = cosf(w0_hs);
    float sqrtA = sqrtf(A_hs);
    float hs_b0 = A_hs * ((A_hs+1.0f) + (A_hs-1.0f)*cos_hs + 2.0f*sqrtA*alpha_hs);
    float hs_b1 = -2.0f*A_hs * ((A_hs-1.0f) + (A_hs+1.0f)*cos_hs);
    float hs_b2 = A_hs * ((A_hs+1.0f) + (A_hs-1.0f)*cos_hs - 2.0f*sqrtA*alpha_hs);
    float hs_a0 = (A_hs+1.0f) - (A_hs-1.0f)*cos_hs + 2.0f*sqrtA*alpha_hs;
    float hs_a1 = 2.0f * ((A_hs-1.0f) - (A_hs+1.0f)*cos_hs);
    float hs_a2 = (A_hs+1.0f) - (A_hs-1.0f)*cos_hs - 2.0f*sqrtA*alpha_hs;
    /* 归一化到 a0 */
    float hs_b0n = hs_b0 / hs_a0, hs_b1n = hs_b1 / hs_a0, hs_b2n = hs_b2 / hs_a0;
    float hs_a1n = hs_a1 / hs_a0, hs_a2n = hs_a2 / hs_a0;
    float hs_x1 = 0.0f, hs_x2 = 0.0f, hs_y1 = 0.0f, hs_y2 = 0.0f;

    for (uint32_t i = 0; i < total_out; i++) {
        float norm = (out_norm[i] > 1e-6f) ? out_norm[i] : 1.0f;
        float vad_ratio = out_vad[i] / (norm + 1e-10f);
        float val = out_float[i] / ((float)input_len * norm);

        // ① 六阶高通 80Hz（3 级级联），压制 50~80Hz 共振
        float out_hp1 = hp1_b0*val + hp1_b1*hp1_x1 + hp1_b2*hp1_x2
                       - hp1_a1*hp1_y1 - hp1_a2*hp1_y2;
        hp1_x2 = hp1_x1; hp1_x1 = val;
        hp1_y2 = hp1_y1; hp1_y1 = out_hp1;

        float out_hp2 = hp2_b0*out_hp1 + hp2_b1*hp2_x1 + hp2_b2*hp2_x2
                       - hp2_a1*hp2_y1 - hp2_a2*hp2_y2;
        hp2_x2 = hp2_x1; hp2_x1 = out_hp1;
        hp2_y2 = hp2_y1; hp2_y1 = out_hp2;

        float out_hp3 = hp3_b0*out_hp2 + hp3_b1*hp3_x1 + hp3_b2*hp3_x2
                       - hp3_a1*hp3_y1 - hp3_a2*hp3_y2;
        hp3_x2 = hp3_x1; hp3_x1 = out_hp2;
        hp3_y2 = hp3_y1; hp3_y1 = out_hp3;
        val = out_hp3;

        // ② 50Hz 窄带陷波 —— 消除工频基波
        float out_n50 = n50_b0*val + n50_b1*n50_x1 + n50_b2*n50_x2
                      - n50_a1*n50_y1 - n50_a2*n50_y2;
        n50_x2 = n50_x1; n50_x1 = val;
        n50_y2 = n50_y1; n50_y1 = out_n50;
        val = out_n50;

        // ③ 100Hz 窄带陷波 —— 消除工频二次谐波
        float out_n100 = n100_b0*val + n100_b1*n100_x1 + n100_b2*n100_x2
                       - n100_a1*n100_y1 - n100_a2*n100_y2;
        n100_x2 = n100_x1; n100_x1 = val;
        n100_y2 = n100_y1; n100_y1 = out_n100;
        val = out_n100;

        // ④ 1.2kHz 高频搁架 +6dB —— VAD 门控：仅语音段提升齿音
        if (vad_ratio > 0.1f) {
            float out_hs = hs_b0n*val + hs_b1n*hs_x1 + hs_b2n*hs_x2
                         - hs_a1n*hs_y1 - hs_a2n*hs_y2;
            hs_x2 = hs_x1; hs_x1 = val;
            hs_y2 = hs_y1; hs_y1 = out_hs;
            val = out_hs;
        }
        // 静音段：旁路高搁架，不更新 IIR 状态，不放大高频底噪

        // 幅度恢复 + 钳位
        val *= 1.5f;
        if (val >  32767.0f) val =  32767.0f;
        if (val < -32768.0f) val = -32768.0f;
        out[i] = (int16_t)val;
    }

    free(msc_S11); free(msc_S22); free(msc_S12_re); free(msc_S12_im); free(msc_gain);

    free(out_vad);

    free(out_float);
    free(out_norm);
    *outLen = total_out;
    return out;
}


/* ========== 三个延迟估计函数（各自独立封装） ========== */

/**
 * @brief 方法1：时域互相关延迟估计
 */
static float estimate_delay_time_domain(const int16_t* data1, uint32_t len1,
                                         const int16_t* data2, uint32_t len2)
{
    int max_delay = 5000;
    if ((int)len1 - 1 < max_delay) max_delay = (int)len1 - 1;
    if ((int)len2 - 1 < max_delay) max_delay = (int)len2 - 1;
    if (max_delay < 0) max_delay = 0;
    float delay = (float)estimate_delay(data1, len1, data2, len2, max_delay);
    printf("Time-domain estimated delay: %.4f samples\n", delay);
    return delay;
}

/**
 * @brief 方法2：FFT-PHAT (legacy) 延迟估计
 */
static float estimate_delay_fft_phat(const int16_t* data1, uint32_t len1,
                                      const int16_t* data2, uint32_t len2)
{
    uint32_t min_len = (len1 < len2) ? len1 : len2;

    float* f1 = (float*)malloc(len1 * sizeof(float));
    float* f2 = (float*)malloc(len2 * sizeof(float));
    if (!f1 || !f2) {
        fprintf(stderr, "Memory allocation failed\n");
        free(f1); free(f2);
        return 0.0f;
    }
    for (uint32_t i = 0; i < len1; i++) f1[i] = (float)data1[i] / 32768.0f;
    for (uint32_t i = 0; i < len2; i++) f2[i] = (float)data2[i] / 32768.0f;

    int fft_size = 512;
    int window = (min_len < (uint32_t)fft_size) ? min_len : (uint32_t)fft_size;
    int margin = 200;
    int peak_num = 1;
    int* delays = (int*)malloc(sizeof(int));
    float* peak_values = (float*)malloc(sizeof(float));
    if (!delays || !peak_values) {
        fprintf(stderr, "Memory allocation failed\n");
        free(f1); free(f2); free(delays); free(peak_values);
        return 0.0f;
    }

    FFT_Real_Gcc_Path(delays, peak_values, &peak_num,
                      f2, f1, margin, window, fft_size);

    float estimated_delay = 0.0f;
    if (peak_num > 0) {
        estimated_delay = (float)delays[0];
        printf("FFT-PHAT estimated delay: %.4f samples (confidence: %.6f)\n",
               estimated_delay, peak_values[0]);
    } else {
        printf("FFT-PHAT failed, falling back to time-domain method\n");
        estimated_delay = estimate_delay(data1, len1, data2, len2, 5000);
    }

    free(f1); free(f2); free(delays); free(peak_values);
    return estimated_delay;
}

/**
 * @brief 方法3：GCC-PHAT (context-based) 延迟估计
 *
 * 将 int16_t 输入转为 float 后委托给 gcc_phat_compute，
 * gcc_phat_compute 内部完成：
 *   DC去除 → Hamming加窗 → FFT → PHAT归一化 → 300~3400Hz带通掩码层1
 *   → IFFT → 峰值搜索 → 亚采样抛物线插值
 * 返回后 ctx->out1/out2 保留FFT结果，供 freq_domain_beamform 做频域波束成形。
 */
float estimate_delay_gcc_phat(const int16_t* data1, uint32_t len1,
                                 const int16_t* data2, uint32_t len2,
                                 int sample_rate, GccPhatContext** out_ctx)
{
    // 选用1024点FFT，可按需调整
    int fft_size = 1024;
    GccPhatContext* ctx = malloc(sizeof(GccPhatContext));
    if(gcc_phat_init(ctx, fft_size) != 0)
    {
        free(ctx);
        *out_ctx = NULL;
        return NAN;
    }

    // 将int16转float临时数组
    uint32_t copy_len = (len1 < len2) ? len1 : len2;
    int window_len = (copy_len < fft_size) ? copy_len : fft_size;

    float* f1 = (float*)malloc(window_len * sizeof(float));
    float* f2 = (float*)malloc(window_len * sizeof(float));
    if(!f1 || !f2)
    {
        free(f1); free(f2);
        gcc_phat_destroy(ctx);
        free(ctx);
        *out_ctx = NULL;
        return NAN;
    }
    for(int i = 0; i < window_len; i++)
    {
        f1[i] = (float)data1[i];
        f2[i] = (float)data2[i];
    }

    // gcc_phat_compute 内部：DC去除/加窗/FFT/PHAT/带通掩码/IFFT/峰值搜索+亚采样插值
    float estimated_delay = gcc_phat_compute(ctx, f1, f2, window_len, fft_size / 2, sample_rate);

    free(f1);
    free(f2);

    // ctx->out1/out2 保留FFT结果，供 freq_domain_beamform 使用
    *out_ctx = ctx;
    return estimated_delay;
}

/**
 * @brief 延迟估计入口
 *
 * 调用具体方法估算两路信号间延迟。
 * 将 int16_t 转换为 float 的逻辑封装在各方法内部，main 无需关心。
 *
 * == 三种方法接口 ==
 *
 *   // 1. 时域互相关
 *   // float delay = estimate_delay_time_domain(data1, len1, data2, len2);
 *
 *   // 2. FFT-PHAT (legacy)
 *   // float delay = estimate_delay_fft_phat(data1, len1, data2, len2);
 *
 *   // 3. GCC-PHAT (当前使用)
 *   // float delay = estimate_delay_gcc_phat(data1, len1, data2, len2, sample_rate);
 *
 * 
 */
float estimate_delay_interactive(const int16_t* data1, uint32_t len1,
                                 const int16_t* data2, uint32_t len2,
                                 int sample_rate, GccPhatContext** out_ctx)
{
    /* 当前使用 ---- GCC-PHAT ----------------------------------- */
    float estimated_delay = estimate_delay_gcc_phat(data1, len1, data2, len2, sample_rate, out_ctx);
    /* ---------------------------------------------------------- */

    /*
    // 备选1：时域互相关
    float estimated_delay = estimate_delay_time_domain(data1, len1, data2, len2);
    */

    /*
    // 备选2：FFT-PHAT (legacy)
    float estimated_delay = estimate_delay_fft_phat(data1, len1, data2, len2);
    */

    printf("Estimated relative delay: %.4f samples (positive = second signal lags)\n", estimated_delay);
    return estimated_delay;
}
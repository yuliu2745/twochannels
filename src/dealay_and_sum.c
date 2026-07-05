#include "../include/setting.h"
#include "../include/fft_path.h"
#include "../include/gcc_phat_delay.h"

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

    // 人声带通边界
    float bin_hz = (float)fs / (float)input_len;
    int bin_low  = (int)ceilf(300.0f / bin_hz);
    int bin_high = (int)floorf(3400.0f / bin_hz);

    uint32_t min_len = (len1 < len2) ? len1 : len2;
    if (min_len == 0) { *outLen = 0; return NULL; }

    // 帧数计算
    int nframes = (min_len + hop - 1) / hop;
    if (nframes < 1) nframes = 1;
    uint32_t total_out = (nframes - 1) * hop + fft_size;

    // 浮点累加缓冲 + 归一化权重累加
    float* out_float = (float*)calloc(total_out, sizeof(float));
    float* out_norm  = (float*)calloc(total_out, sizeof(float));
    if (!out_float || !out_norm) {
        free(out_float); free(out_norm);
        *outLen = 0; return NULL;
    }

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
            float omega = 2.0f * (float)M_PI * k / (float)input_len;
            float phase = omega * delay;
            float c = cosf(phase), s = sinf(phase);
            float yr = mic2_align[k][0], yi = mic2_align[k][1];
            mic2_align[k][0] = yr * c + yi * s;
            mic2_align[k][1] = -yr * s + yi * c;
        }

        // DSB 频域相加 + 300~3400Hz 带通掩码（层2）
        for (int k = 0; k < nbins; k++) {
            float r = ctx->out1[k][0] + mic2_align[k][0];
            float i = ctx->out1[k][1] + mic2_align[k][1];
            if (k < bin_low || k > bin_high) {
                ctx->cs[k][0] = 0.0f;
                ctx->cs[k][1] = 0.0f;
            } else {
                ctx->cs[k][0] = r;
                ctx->cs[k][1] = i;
            }
        }

        // IFFT
        fftwf_execute(ctx->plan_inv);

        // 重叠相加：取 xcorr[0 : fft_size]（后半段是循环混叠）
        for (int i = 0; i < fft_size; i++) {
            int pos = start + i;
            if (pos < (int)total_out) {
                out_float[pos] += ctx->xcorr[i];
                out_norm[pos]  += ctx->window_coeffs[i % fft_size];
            }
        }
    }

    // 归一化：xcorr = input_len * 2ch * signal * window_overlap_sum
    // 所以 signal = xcorr / (input_len * 2 * norm_sum)
    int16_t* out = (int16_t*)malloc(total_out * sizeof(int16_t));
    if (!out) { free(out_float); free(out_norm); *outLen = 0; return NULL; }

    for (uint32_t i = 0; i < total_out; i++) {
        float norm = (out_norm[i] > 1e-6f) ? out_norm[i] : 1.0f;
        float gain = 1.5f;
        float val = out_float[i] * gain / ((float)input_len * 2.0f * norm);
        if (val >  32767.0f) val =  32767.0f;
        if (val < -32768.0f) val = -32768.0f;
        out[i] = (int16_t)val;
    }

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
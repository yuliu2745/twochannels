#include "mmse_lsa.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void lsa_init(MmseLsaCtx* lsa, int nbins)
{
    memset(lsa, 0, sizeof(MmseLsaCtx));
    lsa->bin_count = nbins;
    lsa->alpha_noise = 0.80f;           // 静音帧：快速跟踪纯噪声基底
    lsa->alpha_noise_speech = 0.98f;    // 语音帧：极慢更新，不被融合带入噪声拉高
    lsa->alpha_speech = 0.85f;
    lsa->vad_snr_thresh = 4.0f;         // 高门限，仅强语音划入 speech 帧
    lsa->max_attenuation = 0.03f;       // -30dB 安全下限
    // hs_bin_* 全由 memset 置 0（默认关闭）
    lsa->initialized = 0;
}

void lsa_process_frame(MmseLsaCtx* lsa, fftwf_complex* spec, const fftwf_complex* ref_snr)
{
    int nbins = lsa->bin_count;
    float frame_pow[MAX_FFT_BINS];     // spec 功率（用于噪声跟踪）
    float ref_pow[MAX_FFT_BINS];       // ref_snr 功率（用于 SNR 计算，不被融合注入抬高）
    float total_snr = 0.0f;

    // 两个基准分别取功率
    for (int k = 0; k < nbins; k++)
    {
        float re_s = spec[k][0], im_s = spec[k][1];
        frame_pow[k] = re_s * re_s + im_s * im_s;
        float re_r = ref_snr[k][0], im_r = ref_snr[k][1];
        ref_pow[k] = re_r * re_r + im_r * im_r;
    }

    // 首帧用 ref_snr 初始化噪声谱（首帧无融合注入问题）
    if (!lsa->initialized)
    {
        for (int k = 0; k < nbins; k++)
            lsa->noise_pow[k] = ref_pow[k];
        lsa->initialized = 1;
        return;
    }

    // 全局 VAD：用 ref_pow（不被融合注入影响）
    for (int k = 0; k < nbins; k++)
    {
        float snr = ref_pow[k] / (lsa->noise_pow[k] + 1e-10f);
        total_snr += snr;
    }
    float avg_snr = total_snr / nbins;

    // 语音→静音快速过渡：γ_avg 从 ≥1 跌落 <1 后，连续 6 帧强制 α=0.75
    if (lsa->initialized) {
        if (avg_snr < 1.0f && lsa->prev_avg_snr >= 1.0f) {
            lsa->noise_fast_frames = 6;
            lsa->shelf_hangover = 4;    // 搁架强制关闭 4 帧
            lsa->shelf_fadein   = 0;
        }
    }

    // 噪声跟踪系数 α
    float alpha_n;
    if (lsa->noise_fast_frames > 0) {
        alpha_n = 0.75f;                     // 语音→静音快速追底
    } else if (avg_snr < 0.8f) {
        alpha_n = lsa->alpha_noise;          // 0.80
    } else if (avg_snr > 1.5f) {
        alpha_n = lsa->alpha_noise_speech;   // 0.98
    } else {
        float t = (avg_snr - 0.8f) / 0.7f;
        alpha_n = lsa->alpha_noise + t * (lsa->alpha_noise_speech - lsa->alpha_noise);
    }
    for (int k = 0; k < nbins; k++)
    {
        lsa->noise_pow[k] = alpha_n * lsa->noise_pow[k]
                          + (1.0f - alpha_n) * frame_pow[k];
    }

    // 起音检测（慢启动 TEST B）：avg_snr 从静默跳升到 >2.0 判为起音
    if (lsa->initialized) {
        if (avg_snr > 2.0f && lsa->prev_avg_snr < 1.2f)
            lsa->onset_frames = 8;
    }
    lsa->prev_avg_snr = avg_snr;

    /* ---- 搁架 Hangover 帧计数 ---- */
    if (lsa->shelf_hangover > 0) {
        lsa->shelf_hangover--;
        if (lsa->shelf_hangover == 0)
            lsa->shelf_fadein = 4;      // 开始 4 帧线性淡入
    }
    if (lsa->shelf_fadein > 0)
        lsa->shelf_fadein--;

    /* 搁架帧级使能 + 淡出因子 */
    int shelf_frame_active;
    float shelf_fade_factor = 1.0f;
    if (lsa->shelf_hangover > 0) {
        shelf_frame_active = 0;                    // 强制全关
    } else if (lsa->shelf_fadein > 0) {
        shelf_frame_active = 1;
        shelf_fade_factor = (float)(4 - lsa->shelf_fadein) / 4.0f;  // 0→1 线性淡入
    } else {
        shelf_frame_active = (avg_snr >= 1.0f) ? 1 : 0;
    }

    // ---- Phase 1: 计算原始增益（逐频点） ----
    float raw_gain[MAX_FFT_BINS];
    for (int k = 0; k < nbins; k++)
    {
        float sig_pow = ref_pow[k];
        float noi_pow = lsa->noise_pow[k] + 1e-10f;
        float gamma_k = sig_pow / noi_pow;

        // 线性渐变增益：γ≤1.0→0.06, 1.0<γ<3.0 线性渐变至 1.0, γ≥3.0→1.0
        float gain;
        if (gamma_k > 3.0f) {
            gain = 1.0f;
        } else if (gamma_k > 1.0f) {
            float t = (gamma_k - 1.0f) / 2.0f;
            gain = 0.06f + t * (1.0f - 0.06f);
        } else {
            gain = 0.06f;
        }

        // 高频搁架 — 三段平滑过渡（γ>1.5 子带生效，帧级 hangover 控制）
        float shelf_boost = 1.0f;
        if (gamma_k > 1.5f && lsa->hs_bin_500 > 0 && k >= lsa->hs_bin_500 && k < lsa->hs_bin_3400) {
            if (k < lsa->hs_bin_1000) {
                float t = (float)(k - lsa->hs_bin_500)
                        / (float)(lsa->hs_bin_1000 - lsa->hs_bin_500);
                float sin_t = sinf(t * (float)M_PI * 0.5f);
                shelf_boost = (1.0f + 0.8f * sin_t * sin_t);
            } else if (k < lsa->hs_bin_2500) {
                float t = (float)(k - lsa->hs_bin_1000)
                        / (float)(lsa->hs_bin_2500 - lsa->hs_bin_1000);
                shelf_boost = (1.8f + t * 0.2f);
            } else {
                float t = (float)(k - lsa->hs_bin_2500)
                        / (float)(lsa->hs_bin_3400 - lsa->hs_bin_2500);
                float sin_t = sinf((1.0f - t) * (float)M_PI * 0.5f);
                shelf_boost = (1.0f + 1.0f * sin_t * sin_t);
            }
        }
        /* 帧级静音锁 + 淡出 Hangover */
        if (!shelf_frame_active) {
            shelf_boost = 1.0f;
        } else if (shelf_fade_factor < 1.0f) {
            shelf_boost = 1.0f + (shelf_boost - 1.0f) * shelf_fade_factor;
        }
        gain *= shelf_boost;

        // 安全钳位
        if (gain < lsa->max_attenuation)
            gain = lsa->max_attenuation;

        raw_gain[k] = gain;
    }

    // ---- Phase 2: 频域子带间 3 点邻域滑动平均（合并总增益后平滑） ----
    float freq_sm_gain[MAX_FFT_BINS];
    freq_sm_gain[0] = raw_gain[0] * 0.75f + raw_gain[1] * 0.25f;
    for (int k = 1; k < nbins - 1; k++) {
        freq_sm_gain[k] = 0.2f * raw_gain[k-1] + 0.6f * raw_gain[k] + 0.2f * raw_gain[k+1];
    }
    freq_sm_gain[nbins-1] = raw_gain[nbins-2] * 0.25f + raw_gain[nbins-1] * 0.75f;

    // ---- Phase 3: 帧间增益平滑（非对称 + 高频收紧） ----
    for (int k = 0; k < nbins; k++)
    {
        float gain = freq_sm_gain[k];

        float gain_sm;
        if (gain > lsa->gain_smooth_prev[k]) {
            gain_sm = 0.2f * lsa->gain_smooth_prev[k] + 0.8f * gain;  // 快攻（始终不变）
        } else {
            // 分频段帧间平滑：低频保瞬态，中高频削颗粒感，高频彻底压毛刺
            float freq_k = (float)k / (float)nbins * 8000.0f;  // 近似 bin_hz，16kHz/2
            float release_alpha;
            if (freq_k > 3000.0f) {
                release_alpha = 0.90f;   // >3kHz 最强平滑，彻底压制毛刺
            } else if (freq_k > 1000.0f) {
                release_alpha = 0.85f;   // 1~3kHz 加强平滑，颗粒感集中区
            } else {
                release_alpha = 0.70f;   // <1kHz 不变，保瞬态
            }
            gain_sm = release_alpha * lsa->gain_smooth_prev[k] + (1.0f - release_alpha) * gain;
        }
        lsa->gain_smooth_prev[k] = gain_sm;
        gain = gain_sm;

        spec[k][0] *= gain;
        spec[k][1] *= gain;
    }

    if (lsa->onset_frames > 0)
        lsa->onset_frames--;
    if (lsa->noise_fast_frames > 0)
        lsa->noise_fast_frames--;
}

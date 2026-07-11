#include "mmse_lsa.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void lsa_init(MmseLsaCtx* lsa, int nbins)
{
    memset(lsa, 0, sizeof(MmseLsaCtx));
    lsa->bin_count = nbins;
    lsa->alpha_noise = 0.80f;           // 静音帧：快速跟踪纯噪声基底
    lsa->alpha_noise_speech = 0.92f;    // 语音帧：缓慢更新，不跟踪人声波动
    lsa->alpha_speech = 0.85f;
    lsa->vad_snr_thresh = 4.0f;         // 高门限，仅强语音划入 speech 帧
    lsa->max_attenuation = 0.03f;       // -30dB 安全下限
    lsa->hs_bin_ramp_start = 0;     // 默认关闭（0=不生效）
    lsa->hs_bin_full = 0;
    lsa->initialized = 0;
}

void lsa_process_frame(MmseLsaCtx* lsa, fftwf_complex* spec)
{
    int nbins = lsa->bin_count;
    float frame_pow[MAX_FFT_BINS];
    float total_snr = 0.0f;

    for (int k = 0; k < nbins; k++)
    {
        float re = spec[k][0];
        float im = spec[k][1];
        frame_pow[k] = re * re + im * im;
    }

    // 首帧用自身功率初始化噪声谱
    if (!lsa->initialized)
    {
        for (int k = 0; k < nbins; k++)
            lsa->noise_pow[k] = frame_pow[k];
        lsa->initialized = 1;
        return;
    }

    // 全局 VAD：平均 SNR
    for (int k = 0; k < nbins; k++)
    {
        float snr = frame_pow[k] / (lsa->noise_pow[k] + 1e-10f);
        total_snr += snr;
    }
    float avg_snr = total_snr / nbins;
    int is_speech_frame = (avg_snr > lsa->vad_snr_thresh) ? 1 : 0;

    // 按帧类型选择噪声更新速度（每帧都更新，避免静音过长导致噪声基底过时）
    float alpha_n = is_speech_frame ? lsa->alpha_noise_speech : lsa->alpha_noise;
    for (int k = 0; k < nbins; k++)
    {
        lsa->noise_pow[k] = alpha_n * lsa->noise_pow[k]
                          + (1.0f - alpha_n) * frame_pow[k];
    }

    // 逐子带 SNR 分档增益（核心：三档动态，无全局固定倍率）
    for (int k = 0; k < nbins; k++)
    {
        float sig_pow = frame_pow[k];
        float noi_pow = lsa->noise_pow[k] + 1e-10f;

        float gamma_k = sig_pow / noi_pow;          // 后验 SNR

        float gain;
        if (gamma_k > 3.0f) {
            // 高 SNR：语音主导，完全保留人声
            gain = 1.0f;
        } else if (gamma_k > 1.0f) {
            // 中 SNR：弱辅音 / 尾音，轻微降噪保留共振峰
            gain = 0.85f;
        } else {
            // 低 SNR：纯噪声，深度压制
            gain = 0.12f;
        }

        // 高频搁架渐变提亮 900→1200Hz（仅 γ>1 语音子带生效）
        if (gamma_k > 1.0f && lsa->hs_bin_ramp_start > 0 && k >= lsa->hs_bin_ramp_start) {
            if (k >= lsa->hs_bin_full) {
                gain *= 2.0f;
            } else {
                float t = (float)(k - lsa->hs_bin_ramp_start)
                        / (float)(lsa->hs_bin_full - lsa->hs_bin_ramp_start);
                gain *= (1.0f + t);  // 1.0 → 2.0
            }
        }

        // 安全钳位
        if (gain < lsa->max_attenuation)
            gain = lsa->max_attenuation;

        spec[k][0] *= gain;
        spec[k][1] *= gain;
    }
}

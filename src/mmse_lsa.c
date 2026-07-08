#include "mmse_lsa.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void lsa_init(MmseLsaCtx* lsa, int nbins)
{
    memset(lsa, 0, sizeof(MmseLsaCtx));
    lsa->bin_count = nbins;
    lsa->alpha_noise = 0.95f;      // 平滑更新避免语音泄漏
    lsa->alpha_speech = 0.85f;
    lsa->vad_snr_thresh = 2.0f;    // 中门限（原 3.0），弱辅音帧不被误判为噪声
    lsa->max_attenuation = 0.10f;  // -20dB，温和降噪，减少吞音
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

    // 仅静音帧更新噪声功率谱
    if (!is_speech_frame)
    {
        for (int k = 0; k < nbins; k++)
        {
            lsa->noise_pow[k] = lsa->alpha_noise * lsa->noise_pow[k]
                              + (1.0f - lsa->alpha_noise) * frame_pow[k];
        }
    }

    // MMSE-LSA 增益
    for (int k = 0; k < nbins; k++)
    {
        float sig_pow = frame_pow[k];
        float noi_pow = lsa->noise_pow[k] + 1e-10f;

        float gamma_k = sig_pow / noi_pow;          // 后验 SNR
        float xi_k = gamma_k;                        // 先验 SNR（简化，后置降噪够用）

        float v = gamma_k * xi_k / (1.0f + xi_k);
        float gain;

        if (v > 50.0f) {
            gain = 1.0f;
        } else if (xi_k < 1e-6f) {
            gain = lsa->max_attenuation;
        } else {
            gain = (xi_k / (1.0f + xi_k)) * expf(-0.5f * v);
        }

        // 钳位最小增益，防止过度抑制
        if (gain < lsa->max_attenuation)
            gain = lsa->max_attenuation;

        // 静音帧额外压低
        if (!is_speech_frame)
            gain *= 0.5f;

        spec[k][0] *= gain;
        spec[k][1] *= gain;
    }
}

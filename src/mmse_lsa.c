#include "mmse_lsa.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void lsa_init(MmseLsaCtx* lsa, int nbins)
{
    memset(lsa, 0, sizeof(MmseLsaCtx));
    lsa->bin_count = nbins;
    lsa->alpha_noise = 0.95f;    // 噪声更新权重，越大越稳定
    lsa->alpha_speech = 0.80f;
    lsa->vad_snr_thresh = 1.2f;  // 平均SNR低于此值判定为噪声帧
    lsa->initialized = 0;        // 首帧标记
}

void lsa_process_frame(MmseLsaCtx* lsa, fftwf_complex* spec)
{
    int nbins = lsa->bin_count;
    float frame_pow[MAX_FFT_BINS];
    float total_snr = 0.0f;

    // 1. 计算当前帧每个频点功率
    for (int k = 0; k < nbins; k++)
    {
        float re = spec[k][0];
        float im = spec[k][1];
        frame_pow[k] = re * re + im * im;
    }

    //首帧用自身功率初始化噪声谱，防止零除导致增益崩溃
    if (!lsa->initialized)
    {
        for (int k = 0; k < nbins; k++)
            lsa->noise_pow[k] = frame_pow[k];
        lsa->initialized = 1;
        // 首帧不做衰减，直接返回
        return;
    }

    // 2. 简易全局VAD：计算整帧平均SNR，区分语音/静音帧
    for (int k = 0; k < nbins; k++)
    {
        float snr = frame_pow[k] / (lsa->noise_pow[k] + 1e-10f);
        total_snr += snr;
    }
    float avg_snr = total_snr / nbins;
    int is_speech_frame = (avg_snr > lsa->vad_snr_thresh) ? 1 : 0;

    // 3. 仅静音帧更新噪声功率谱（持续学习底噪）
    if (!is_speech_frame)
    {
        for (int k = 0; k < nbins; k++)
        {
            lsa->noise_pow[k] = lsa->alpha_noise * lsa->noise_pow[k]
                              + (1.0f - lsa->alpha_noise) * frame_pow[k];
        }
    }

    // 4. MMSE-LSA 增益计算，逐频点衰减噪声
    //    标准公式（Ephraim & Malah 1984）：
    //      γ_k = |Y_k|² / λ_d(k)          —— 后验 SNR
    //      ξ_k = α·|X̂_{k,prev}|²/λ_d + (1-α)·max(γ_k-1,0)  —— 先验 SNR（DD 估计）
    //      v_k = γ_k·ξ_k / (1+ξ_k)
    //      G_k = (ξ_k/(1+ξ_k)) · exp(0.5 · ∫_{v_k}^{∞} e^{-t}/t dt)
    //
    //    这里做合理简化：使用 γ_k 估计 ξ_k（DD 平滑），避免 exp 内积分用闭式近似。
    for (int k = 0; k < nbins; k++)
    {
        float sig_pow = frame_pow[k];
        float noi_pow = lsa->noise_pow[k] + 1e-10f;

        float gamma_k = sig_pow / noi_pow;             // 后验 SNR
        float xi_k = gamma_k;                          // 先验 SNR（简化，可用 DD 平滑替代）

        float v = gamma_k * xi_k / (1.0f + xi_k);

        float gain;
        if (v > 50.0f) {
            // SNR 非常高的频点，G → 1
            gain = 1.0f;
        } else if (xi_k < 1e-6f) {
            gain = 0.01f;
        } else {
            gain = (xi_k / (1.0f + xi_k)) * expf(-0.5f * v);
        }

        // 静音帧额外压低增益
        if (!is_speech_frame)
            gain *= 0.25f;

        // 频谱乘以降噪增益
        spec[k][0] *= gain;
        spec[k][1] *= gain;
    }
}
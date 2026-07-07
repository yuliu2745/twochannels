#ifndef MMSE_LSA_H
#define MMSE_LSA_H

#include <fftw3.h>
#include <string.h>
#include <math.h>

#define MAX_FFT_BINS 2049

typedef struct
{
    float noise_pow[MAX_FFT_BINS];  // 噪声功率谱
    float alpha_noise;              // 噪声更新平滑系数（静音帧慢更新）
    float alpha_speech;             // 语音平滑系数
    float vad_snr_thresh;           // VAD全局SNR阈值
    float max_attenuation;          // 最大衰减(dB 线性), 0.1=-20dB
    int bin_count;                  // 当前FFT有效频点数量
    int initialized;                // 首帧噪声功率已初始化标记
} MmseLsaCtx;

void lsa_init(MmseLsaCtx* lsa, int nbins);
void lsa_process_frame(MmseLsaCtx* lsa, fftwf_complex* spec);

#endif

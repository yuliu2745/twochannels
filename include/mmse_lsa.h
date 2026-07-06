#ifndef MMSE_LSA_H
#define MMSE_LSA_H

#include <fftw3.h>
#include <string.h>
#include <math.h>

#define MAX_FFT_BINS 2049 // 适配你input_len=2048，complex_len=1025，预留余量

typedef struct
{
    float noise_pow[MAX_FFT_BINS];  // 噪声功率谱
    float alpha_noise;              // 噪声更新平滑系数（静音帧慢更新）
    float alpha_speech;             // 语音平滑系数
    float vad_snr_thresh;           // VAD全局SNR阈值
    int bin_count;                  // 当前FFT有效频点数量
    int initialized;                // 首帧噪声功率已初始化标记
} MmseLsaCtx;

// 初始化LSA上下文
void lsa_init(MmseLsaCtx* lsa, int nbins);
// 单帧频域MMSE-LSA处理，输入输出均为复数频谱ctx->cs
void lsa_process_frame(MmseLsaCtx* lsa, fftwf_complex* spec);

#endif
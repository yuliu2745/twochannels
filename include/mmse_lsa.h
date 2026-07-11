#ifndef MMSE_LSA_H
#define MMSE_LSA_H

#include <fftw3.h>
#include <string.h>
#include <math.h>

#define MAX_FFT_BINS 2049

typedef struct
{
    float noise_pow[MAX_FFT_BINS];  // 噪声功率谱
    float alpha_noise;              // 噪声更新平滑系数（静音帧 0.80，快速跟踪）
    float alpha_noise_speech;       // 语音帧噪声更新平滑系数（0.92，稳定基底）
    float alpha_speech;             // 语音平滑系数
    float vad_snr_thresh;           // VAD全局SNR阈值（4.0，仅强语音判为speech帧）
    float max_attenuation;          // 最低增益钳位(线性) 0.03=-30dB，安全下限
    int bin_count;                  // 当前FFT有效频点数量
    int hs_bin_ramp_start;          // 高频搁架渐升起始bin（900Hz），γ>1时生效
    int hs_bin_full;                // 高频搁架恒定增益起始bin（1200Hz），gain=2.0
    int initialized;                // 首帧噪声功率已初始化标记
} MmseLsaCtx;

void lsa_init(MmseLsaCtx* lsa, int nbins);
void lsa_process_frame(MmseLsaCtx* lsa, fftwf_complex* spec);

#endif

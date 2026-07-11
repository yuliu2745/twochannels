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
    float gain_smooth_prev[MAX_FFT_BINS]; // 逐频点增益平滑历史（帧间递归 0.7）
    int bin_count;                  // 当前FFT有效频点数量
    int hs_bin_500;                 // 高频搁架缓升起始bin（500Hz），γ>1时生效
    int hs_bin_1000;                // 正弦平方缓升结束bin（1000Hz），gain=1.8
    int hs_bin_2500;                // 线性抬升结束bin（2500Hz），gain=2.2
    int hs_bin_3400;                // 正弦缓降截止bin（3400Hz），gain=1.0
    float prev_avg_snr;             // 前一帧 avg_snr，用于起音检测 & 过渡检测
    int    onset_frames;            // 起音慢启动剩余帧数
    int    noise_fast_frames;       // 语音→静音快速过渡剩余帧数（α=0.75）
    int    shelf_hangover;          // 搁架强制关闭剩余帧数（语音转静音后保持静音）
    int    shelf_fadein;            // 搁架淡入剩余帧数（hangover结束后逐步恢复）
    int initialized;                // 首帧噪声功率已初始化标记
} MmseLsaCtx;

void lsa_init(MmseLsaCtx* lsa, int nbins);
void lsa_process_frame(MmseLsaCtx* lsa, fftwf_complex* spec, const fftwf_complex* ref_snr);

#endif

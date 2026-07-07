#ifndef GSC_H
#define GSC_H

#include <fftw3.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/**
 * @brief GSC (Generalized Sidelobe Canceller) context
 *
 * Frequency-domain GSC that operates on aligned mic spectra:
 *   S[k] = (X1[k] + X2_rotated[k]) / 2  — speech reference (DSB)
 *   N[k] = (X1[k] - X2_rotated[k]) / 2  — noise reference (blocking matrix)
 *   Y[k] = S[k] - W[k] * N[k]           — GSC output
 *
 * Adaptive filter W[k] updated via leaky NLMS per frequency bin.
 */
typedef struct {
    fftwf_complex *W;       /* 复数自适应权重 (complex_len) */
    float         *Pn;      /* 噪声参考平滑功率 (complex_len) */
    float         *G_smooth; /* 逐频点维纳增益(平滑后)，用于BM泄漏抑制 */
    float         *S_floor;  /* 语音参考功率噪声底 (complex_len)，缓慢跟踪最小值 */
    float          mu;      /* NLMS 步长 (典型值 0.01~0.03) */
    float          alpha;   /* 功率平滑因子 (典型值 0.9) */
    float          vad_thresh; /* 频点VAD能量比阈值 (β, 典型值 0.4~3.0) */
    float          leak;    /* 语音帧 W 泄漏因子 (0.999~1.0), 1.0=不泄漏 */
    float          smooth_factor; /* G_U 平滑因子(0.8~0.95)，越大G变化越慢 */
    int            complex_len; /* 频点数 */
    int            initialized; /* Pn 首帧已初始化标记 */
} GscContext;

/**
 * @brief 初始化 GSC 上下文
 *
 * @param complex_len  FFT 复数频点数 (fft_size + 1)
 * @param mu           NLMS 步长
 * @param alpha        噪声功率平滑因子
 * @param vad_thresh   频点 VAD 能量比阈值，|S|² > β·|N|² 判为有语音
 * @param leak         语音帧 W 泄漏因子 (0.999~1.0)，1.0=不泄漏
 * @return 初始化成功的指针，失败返回 NULL
 */
GscContext* gsc_init(int complex_len, float mu, float alpha, float vad_thresh, float leak, float smooth_factor);

/**
 * @brief 对一帧频域数据执行 GSC 处理
 *
 * 在每个频点计算：
 *   Speech ref:  S = (X1 + X2) / 2
 *   Noise ref:   N = (X1 - X2) / 2
 *   Y = S - W * N                          (GSC 输出)
 *   更新 W 用泄漏 NLMS（仅当该频点 VAD 判为无语音时）
 *
 * @param gsc   GSC 上下文
 * @param X1    通道1 频谱 (complex_len)
 * @param X2    通道2 对齐后频谱 (complex_len)
 * @param Y     输出频谱 (complex_len)，可与 X1/X2 重叠
 * @param nbins 实际处理的频点数（建议传 complex_len）
 */
void gsc_process_frame(GscContext *gsc,
                       const fftwf_complex *X1,
                       const fftwf_complex *X2,
                       fftwf_complex *Y,
                       int nbins);

/**
 * @brief 销毁 GSC 上下文，释放所有资源
 */
void gsc_destroy(GscContext *gsc);

#endif /* GSC_H */

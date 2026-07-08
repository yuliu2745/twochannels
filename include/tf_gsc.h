#ifndef TF_GSC_H
#define TF_GSC_H

#include <fftw3.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/**
 * @brief TF_GSC (Transfer Function GSC) context
 *
 * Frequency-domain TF-GSC that operates on aligned mic spectra:
 *   S[k] = (X1[k] + X2_rotated[k]) / 2  — speech reference (DSB)
 *   N[k] = (X1[k] - X2_rotated[k]) / 2  — noise reference (blocking matrix)
 *   Y[k] = S[k] - W[k] * N[k]           — TF-GSC output
 *
 * Adaptive filter W[k] updated via leaky NLMS per frequency bin.
 * BM leakage suppression via per-bin Wiener gain G_U.
 */
typedef struct {
    fftwf_complex *W;       /* 复数自适应权重 (complex_len) */
    float         *Pn;      /* 噪声参考平滑功率 (complex_len) */
    float         *G_smooth; /* 逐频点维纳增益(平滑后)，用于BM泄漏抑制 */
    float         *S_floor;  /* 语音参考功率噪声底 (complex_len)，缓慢跟踪最小值 */
    float          mu;      /* NLMS 基准步长 (典型值 0.01~0.03) */
    float          alpha;   /* 功率平滑因子 (典型值 0.9) */
    float          vad_thresh; /* 频点VAD能量比阈值 (β, 典型值 0.4~3.0) */
    float          leak;    /* 语音帧 W 泄漏因子 (0.999~1.0), 1.0=不泄漏 */
    float          smooth_factor; /* G_U 平滑因子(0.8~0.95)，越大G变化越慢 */
    int            complex_len; /* 频点数 */
    int            initialized; /* Pn 首帧已初始化标记 */
    float          prev_fbf_energy; /* 前一帧 FBF 平滑能量，用于语音起振检测 */
    int            onset_frames;    /* 剩余起振保护帧数（检测到起振后计数） */
    float          onset_thresh;    /* 起振能量比阈值（默认 3.0） */

    /* AMC 自适应模式控制器 */
    float          amc_mu;      /* 当前帧有效 NLMS 步长（由 AMC 根据场景设定） */
    float          W_max;       /* 权重幅值钳位上界 |W[k]| < W_max (推荐 0.12~0.18) */
    float          coh_smooth;  /* 帧级双麦相干系数平滑值 (0~1) */

    /* BM 泄漏抑制：语音基频段增益保护 */
    int            gU_low_bin;   /* 保护起始 bin (80Hz) */
    int            gU_high_bin;  /* 保护结束 bin (800Hz) */
    float          gU_min;       /* 最小增益钳位 (推荐 0.9) */

    /* ============================================================ */
    /*  模块 A: 时频扩散度掩膜（预过滤 BM 输入，削减混响多径泄露）       */
    /* ============================================================ */
    float         *diff_P1;       /* 平滑 |X1|² (complex_len) */
    float         *diff_P2;       /* 平滑 |X2|² (complex_len) */
    float         *diff_C12_re;   /* 平滑 Re(X1·X2*) (complex_len) */
    float         *diff_C12_im;   /* 平滑 Im(X1·X2*) (complex_len) */
    float          diff_alpha;    /* 功率谱平滑因子 (推荐 0.85) */
    float          diff_thresh;   /* 扩散度判高阈值 (0.6) — >此值×0.4 */
    float          diff_suppress; /* 高扩散抑制增益 (0.4) */

    /* ============================================================ */
    /*  模块 B: FBF 前置直达语音增强                                    */
    /* ============================================================ */
    int            gfb_enabled;   /* 启用 G_FB(k) 增强标志 */

    /* ============================================================ */
    /*  模块 C: ANC 输出人声下限钳位（VAD 门控，最后一道防线）           */
    /* ============================================================ */
    int            clamp_enabled;     /* 启用标志 */
    float          clamp_min_ratio;   /* 下限比例 (0.35) */
    float          fbf_peak_db;       /* FBF 能量峰值跟踪 (dB) */
    int            clamp_vad;         /* 当前帧内部 VAD (0/1) */
} TfGscContext;

/**
 * @brief 初始化 TF-GSC 上下文
 *
 * @param complex_len  FFT 复数频点数 (fft_size + 1)
 * @param mu           NLMS 步长
 * @param alpha        噪声功率平滑因子
 * @param vad_thresh   频点 VAD 能量比阈值，|S|² > β·|N|² 判为有语音
 * @param leak         语音帧 W 泄漏因子 (0.999~1.0)，1.0=不泄漏
 * @param smooth_factor BM泄漏抑制增益平滑因子
 * @return 初始化成功的指针，失败返回 NULL
 */
TfGscContext* tf_gsc_init(int complex_len, float mu, float alpha, float vad_thresh, float leak, float smooth_factor);

/**
 * @brief 对一帧频域数据执行 TF-GSC 处理
 *
 * 在每个频点计算：
 *   Speech ref:  S = (X1 + X2) / 2          (FBF 输出)
 *   Noise ref:   N = (X1 - X2) / 2          (BM 输出)
 *   Ñ = G_U · N                             (BM 泄漏抑制后的净化噪声)
 *   Y = S - W · Ñ                           (TF-GSC 输出)
 *   更新 W 用泄漏 NLMS（仅当该频点 VAD 判为无语音时）
 *
 * @param tf_gsc TF-GSC 上下文
 * @param X1    通道1 频谱 (complex_len)
 * @param X2    通道2 对齐后频谱 (complex_len)
 * @param Y_FBF [out] FBF 输出 S（可为 NULL，传 NULL 不输出）
 * @param Y_U   [out] 净化噪声参考 Ñ（可为 NULL，传 NULL 不输出）
 * @param Y     [out] TF-GSC 最终输出 (complex_len)
 * @param nbins 实际处理的频点数（建议传 complex_len）
 */
void tf_gsc_process_frame(TfGscContext *tf_gsc,
                          const fftwf_complex *X1,
                          const fftwf_complex *X2,
                          fftwf_complex *Y_FBF,
                          fftwf_complex *Y_U,
                          fftwf_complex *Y,
                          int nbins);

/**
 * @brief 销毁 TF-GSC 上下文，释放所有资源
 */
void tf_gsc_destroy(TfGscContext *tf_gsc);

#endif /* TF_GSC_H */

#ifndef SIGNAL_RESTORE_H
#define SIGNAL_RESTORE_H

#include <fftw3.h>

/**
 * @brief 信号恢复模块上下文
 *
 * 三段式余弦相似度补偿：
 *   SD=|S_X1_Y_FBF - S_X2_Y_FBF| > η      → 选择匹配度高的单麦
 *   η_low < SD < η                         → (X1+X2)/2 混叠补偿，绝不闭锁
 *   SD < η_low                             → 关闭补偿，仅保留 FBF 输出
 *
 * 分 SNR 自适应 beta：
 *   帧 SNR > snr_thresh_db  → beta=beta_high (5~6)，适度增强
 *   帧 SNR ≤ snr_thresh_db  → beta=beta_low (10)，锁死增益，不回填底噪
 */
typedef struct {
    /* 参数 */
    float lambda;           /* 余弦相似度时间平滑因子 (0~1), 推荐 0.88~0.90 */
    float eta;              /* 上阈值：高于此选择单麦 (推荐 0.08~0.12) */
    float eta_low;          /* 下阈值：低于此关闭补偿 (推荐 0.02~0.03) */
    float beta_low;         /* 低SNR帧 G_X 分母权重 (推荐 10.0) */
    float beta_high;        /* 高SNR帧 G_X 分母权重 (推荐 5.0~6.0) */
    float snr_thresh_db;    /* SNR 门限 (dB)，高于此使用 beta_high (推荐 15.0) */
    int   sample_rate;      /* 音频采样率 (Hz)，用于计算 bin_hz */
    float bin_hz;           /* 每个频点的频率宽度 (Hz/bin) */

    /* 帧级平滑余弦相似度 (状态) */
    float S_X1_Y_FBF;
    float S_X2_Y_FBF;
    float S_X1_Y_U;
    float S_X2_Y_U;

    int complex_len;
    int initialized;
} SignalRestoreContext;

/**
 * @brief 初始化信号恢复模块
 *
 * @param complex_len  频点数
 * @param lambda       平滑因子 (推荐 0.88)
 * @param eta          上决策阈值 (推荐 0.10)
 * @param eta_low      下决策阈值 (推荐 0.03)
 * @param beta_low     低SNR G_X分母权重 (推荐 10.0)
 * @param beta_high    高SNR G_X分母权重 (推荐 5.0)
 * @param snr_thresh_db SNR 门限 dB (推荐 15.0)
 * @return  SignalRestoreContext*  失败返回 NULL
 */
SignalRestoreContext* restore_init(int complex_len, float lambda, float eta, float eta_low,
                                   float beta_low, float beta_high, float snr_thresh_db,
                                   int sample_rate);

/**
 * @brief 对一帧频域数据执行信号恢复
 *
 * 步骤:
 *   1. 帧级余弦相似度 S(Xp, Y_FBF) 和 S(Xp, Y_U)
 *   2. 帧 SNR = 10·log₁₀(|Y_FBF|²/|Y_U|²) → 选择 beta_low 或 beta_high
 *   3. 三段式决策选择 X_selected
 *   4. 逐 bin 增益 G_X, 输出 Ȳ_FBF = Y_FBF + G_X · X_selected
 *
 * @param ctx     信号恢复上下文
 * @param Y_FBF   FBF 输出 S (X1+X2)/2
 * @param Y_U     净化噪声参考 Ñ
 * @param X1      通道 1 频谱
 * @param X2      通道 2 对齐后频谱
 * @param Y_out   补偿后输出 Ȳ_FBF
 * @param nbins   有效频点数
 */
void restore_process_frame(SignalRestoreContext *ctx,
                           const fftwf_complex *Y_FBF,
                           const fftwf_complex *Y_U,
                           const fftwf_complex *X1,
                           const fftwf_complex *X2,
                           fftwf_complex *Y_out,
                           int nbins);

/**
 * @brief 销毁信号恢复模块
 */
void restore_destroy(SignalRestoreContext *ctx);

#endif /* SIGNAL_RESTORE_H */

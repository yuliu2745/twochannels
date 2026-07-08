#ifndef SIGNAL_RESTORE_H
#define SIGNAL_RESTORE_H

#include <fftw3.h>

/**
 * @brief 信号恢复模块上下文
 *
 * 根据两路麦克信号与 FBF/Ȳᵤ 的余弦相似度选择"最安全"的麦克信号，
 * 以增益 G_X 叠加回 FBF 输出，补偿波束成形对期望语音的损伤。
 *
 * 参考文献对应 Fig.2 的决策逻辑。
 */
typedef struct {
    /* 参数 */
    float lambda;           /* 余弦相似度时间平滑因子 (0~1) */
    float eta;              /* 选择决策阈值 (推荐 0.1) */
    float beta;             /* G_X 分母中 |Y_U|² 的权重 (>1 偏保守) */

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
 * @param lambda       平滑因子 (推荐 0.7~0.9)
 * @param eta          决策阈值 (推荐 0.1)
 * @param beta         G_X 保守系数 (推荐 10.0)
 * @return  SignalRestoreContext*  失败返回 NULL
 */
SignalRestoreContext* restore_init(int complex_len, float lambda, float eta, float beta);

/**
 * @brief 对一帧频域数据执行信号恢复
 *
 * 步骤:
 *   1. 帧级余弦相似度 S(Xp, Y_FBF) 和 S(Xp, Y_U)
 *   2. 决策逻辑选择 X_selected ∈ {X1, X2, (X1+X2)/2, 0}
 *   3. 逐 bin 增益 G_X, 输出 Ȳ_FBF = Y_FBF + G_X · X_selected
 *
 * @param ctx     信号恢复上下文
 * @param Y_FBF   FBF 输出 S (X1+X2)/2
 * @param Y_U     净化噪声参考 Ȳᵤ
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
 * @brief 销毁信号恢复上下文
 */
void restore_destroy(SignalRestoreContext *ctx);

#endif /* SIGNAL_RESTORE_H */

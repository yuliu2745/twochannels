#ifndef DELAY_AND_SUM_H
#define DELAY_AND_SUM_H

#include "setting.h"
#include "gcc_phat_delay.h"

//时域接口
int16_t* delay_sum(const int16_t* data1, uint32_t len1, float delay1,
                   const int16_t* data2, uint32_t len2, float delay2,
                   uint32_t* outLen);

//频域相位旋转
int16_t* freq_domain_beamform(GccPhatContext* ctx,
                              const int16_t* data1, uint32_t len1,
                              const int16_t* data2, uint32_t len2,
                              float delay, int fs, uint32_t* outLen);

int estimate_delay(const int16_t* x, uint32_t len_x, const int16_t* y, uint32_t len_y, int max_delay);

/**
 * @brief 交互式延迟估计：菜单选择时域/FFT-PHAT/GCC-PHAT 三种算法
 * @param data1       信号1 (int16_t)
 * @param len1        信号1 长度
 * @param data2       信号2 (int16_t)
 * @param len2        信号2 长度
 * @param sample_rate 采样率（仅用于参数提示）
 * @return 估计的延迟（样点数），正=信号2滞后，负=信号2超前，0=失败
 *         子样点精度（GCC-PHAT 返回浮点结果）
 */
float estimate_delay_gcc_phat(const int16_t* data1, uint32_t len1,
                                 const int16_t* data2, uint32_t len2,
                                 int sample_rate, GccPhatContext** out_ctx);

/**
 * @brief 延迟估计入口（交互式选择算法）
 *
 * 封装具体延迟估计算法调用，当前使用 GCC-PHAT (context-based)。
 *
 * @param data1       信号1 (int16_t)
 * @param len1        信号1 长度
 * @param data2       信号2 (int16_t)
 * @param len2        信号2 长度
 * @param sample_rate 采样率
 * @param out_ctx     输出 GccPhatContext 指针（内含 FFT 结果供后续频域波束成形）
 * @return 估计的延迟（样点数），含亚样点精度
 */
float estimate_delay_interactive(const int16_t* data1, uint32_t len1,
                                 const int16_t* data2, uint32_t len2,
                                 int sample_rate, GccPhatContext** out_ctx);

#endif
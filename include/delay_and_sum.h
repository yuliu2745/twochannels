#ifndef DELAY_AND_SUM_H
#define DELAY_AND_SUM_H

#include "setting.h"

int16_t* delay_sum(const int16_t* data1, uint32_t len1, int delay1, const int16_t* data2, uint32_t len2, int delay2, uint32_t* outLen);
int estimate_delay(const int16_t* x, uint32_t len_x, const int16_t* y, uint32_t len_y, int max_delay);

/**
 * @brief 交互式延迟估计：菜单选择时域/FFT-PHAT/GCC-PHAT 三种算法
 * @param data1       信号1 (int16_t)
 * @param len1        信号1 长度
 * @param data2       信号2 (int16_t)
 * @param len2        信号2 长度
 * @param sample_rate 采样率（仅用于参数提示）
 * @return 估计的延迟（样点数），正=信号2滞后，负=信号2超前，0=失败
 */
int estimate_delay_interactive(const int16_t* data1, uint32_t len1,
                                const int16_t* data2, uint32_t len2,
                                int sample_rate);

#endif
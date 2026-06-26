#ifndef DELAY_AND_SUM_H
#define DELAY_AND_SUM_H

#include "setting.h"

int16_t* delay_sum(const int16_t* data1, uint32_t len1, int delay1, const int16_t* data2, uint32_t len2, int delay2, uint32_t* outLen);
int estimate_delay(const int16_t* x, uint32_t len_x, const int16_t* y, uint32_t len_y, int max_delay);

#endif
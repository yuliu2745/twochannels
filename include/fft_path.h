#ifndef FFT_PATH_H
#define FFT_PATH_H

#include "setting.h"

// FFT-PHAT延迟估计函数
float FFT_Real_Gcc_Path(int *delays, float *peak_values, int *peak_num, 
                        float *channel2, float *channel1, int margin, 
                        int window, int fft_size);

// 查找互相关峰值
float find_nbest_maximums(int *delays, float *values, int amount_max, 
                         float *xcorr_value, int margin, int fftwindow_size);

#endif
#ifndef MERGE_AUDIO_H
#define MERGE_AUDIO_H

#include "setting.h"

// 重采样函数（简单的线性插值）
int16_t* resample(const int16_t* input, uint32_t inputLength, uint32_t inputRate, uint32_t outputRate, uint32_t* outputLength);

// 应用延迟到信号
int16_t* apply_delay(const int16_t* input, uint32_t length, int delay_samples, uint32_t* outputLength);

// 计算TDOA延迟样本数
int calculate_tdoa_delay(float mic_distance, float sound_speed, float angle_degrees, uint32_t sample_rate);

// 合并两个音频信号
int16_t* merge_signals(const int16_t* signal1, const int16_t* signal2, uint32_t length);

// 生成测试用的双通道音频
int generate_stereo_with_tdoa(const char* voiceFile, const char* sineFile,
                              const char* outputFile, float mic_distance,
                              float sound_speed, float angle_voice, float angle_noise);

#endif

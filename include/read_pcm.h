#ifndef READ_PCM_H
#define READ_PCM_H

#include "setting.h"

// 读取裸 PCM 文件（16-bit, mono, little-endian）
int16_t* read_pcm(const char* filename, uint32_t* numSamples);

// 写入裸 PCM 文件（16-bit, mono, little-endian）
int write_pcm(const char* filename, const int16_t* data, uint32_t numSamples);

#endif

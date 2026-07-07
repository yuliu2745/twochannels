#ifndef READ_PCM_H
#define READ_PCM_H

#include "setting.h"

/**
 * @struct audio_t
 * @brief 统一音频数据结构，兼容 WAV 和 PCM
 */
typedef struct {
    int16_t*  data;            // 音频样本数据
    uint32_t  num_samples;     // 样本数
    uint32_t  sample_rate;     // 采样率（Hz）
    uint16_t  bits_per_sample; // 位深度
    uint16_t  num_channels;    // 声道数
} audio_t;

// 读取裸 PCM 文件（16-bit, mono, little-endian）
int16_t* read_pcm(const char* filename, uint32_t* numSamples);

// 写入裸 PCM 文件（16-bit, mono, little-endian）
int write_pcm(const char* filename, const int16_t* data, uint32_t numSamples);

/**
 * @brief 统一加载音频文件（自动检测 WAV / PCM）
 * @param filename          文件路径
 * @param force_sample_rate PCM 文件时使用的采样率（WAV 时忽略）
 * @return audio_t* 成功，NULL 失败（需调用 audio_free 释放）
 *
 * 检测规则：读取文件前 4 字节，若为 "RIFF" 则按 WAV 加载，否则按 PCM 加载。
 * PCM 文件必须指定 force_sample_rate（>0）。
 */
audio_t* audio_load(const char* filename, int force_sample_rate);

/// @brief 释放 audio_load 返回的资源
void audio_free(audio_t* a);

#endif

#ifndef SETTING_H
#define SETTING_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

// WAV文件头结构
typedef struct {
    char     chunkID[4];       // "RIFF"
    uint32_t chunkSize;        // 文件大小-8
    char     format[4];        // "WAVE"
    char     subchunk1ID[4];   // "fmt "
    uint32_t subchunk1Size;    // 16 for PCM
    uint16_t audioFormat;      // 1 for PCM
    uint16_t numChannels;
    uint32_t sampleRate;
    uint32_t byteRate;
    uint16_t blockAlign;
    uint16_t bitsPerSample;
    char     subchunk2ID[4];   // "data"
    uint32_t subchunk2Size;    // 数据大小
} WAVHeader;

// 音频文件输出目录
#define AUDIO_OUTPUT_DIR "audio_files/"

// 默认音频参数
#define DEFAULT_MIC_DISTANCE 0.2f    // 麦克风间距（米）
#define DEFAULT_SOUND_SPEED 343.0f   // 声速（米/秒）
#define DEFAULT_SAMPLE_RATE 16000    // 采样率（Hz）

int do_work(int argc, char* argv[]);

// API函数声明
int split_stereo_api(const char* stereoFile, const char* leftFile, const char* rightFile);
int generate_stereo_with_tdoa_api(const char* voiceFile, const char* sineFile,
                                  const char* outputFile, float mic_distance,
                                  float sound_speed, float angle_voice, float angle_noise);

#endif
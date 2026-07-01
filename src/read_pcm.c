#include "../include/read_pcm.h"
#include "../include/readwav.h"
#include <string.h>

/**
 * @brief 读取裸 PCM 文件（16-bit, mono, little-endian）
 * @param filename  文件路径
 * @param numSamples 输出参数：样本数
 * @return 成功返回 int16_t* 音频数据，失败返回 NULL
 *
 * 裸 PCM 文件不包含 WAV 头，文件数据即为连续的 16-bit 采样值。
 */
int16_t* read_pcm(const char* filename, uint32_t* numSamples) {
    FILE* fp = fopen(filename, "rb");
    if (!fp) {
        fprintf(stderr, "Cannot open PCM file: %s\n", filename);
        return NULL;
    }

    // 获取文件大小
    fseek(fp, 0, SEEK_END);
    long fileSize = ftell(fp);
    fseek(fp, 0, SEEK_SET);

    if (fileSize <= 0) {
        fprintf(stderr, "Empty PCM file: %s\n", filename);
        fclose(fp);
        return NULL;
    }

    if (fileSize % 2 != 0) {
        fprintf(stderr, "Warning: PCM file has odd byte size (%ld), may be incomplete\n", fileSize);
    }

    // 计算样本数（16-bit = 2 bytes per sample）
    *numSamples = fileSize / (uint32_t)sizeof(int16_t);

    int16_t* data = (int16_t*)malloc(fileSize);
    if (!data) {
        fprintf(stderr, "Memory allocation failed for %s\n", filename);
        fclose(fp);
        return NULL;
    }

    size_t bytesRead = fread(data, 1, fileSize, fp);
    if (bytesRead != (size_t)fileSize) {
        fprintf(stderr, "Failed to read PCM data: %s (read %zu/%ld bytes)\n",
                filename, bytesRead, fileSize);
        free(data);
        fclose(fp);
        return NULL;
    }

    fclose(fp);

    printf("  Loaded PCM: %s (%u samples, %ld bytes)\n",
           filename, *numSamples, (unsigned long)fileSize);

    return data;
}

/**
 * @brief 写入裸 PCM 文件（16-bit, mono, little-endian）
 */
int write_pcm(const char* filename, const int16_t* data, uint32_t numSamples) {
    FILE* fp = fopen(filename, "wb");
    if (!fp) {
        fprintf(stderr, "Cannot create PCM file: %s\n", filename);
        return 0;
    }

    size_t bytesToWrite = numSamples * sizeof(int16_t);
    if (fwrite(data, 1, bytesToWrite, fp) != bytesToWrite) {
        fprintf(stderr, "Failed to write PCM data to %s\n", filename);
        fclose(fp);
        return 0;
    }

    fclose(fp);
    printf("  Wrote PCM: %s (%u samples)\n", filename, numSamples);
    return 1;
}

/**
 * @brief 统一加载音频文件，自动检测 WAV / PCM 格式
 *
 * 检测规则：读取文件前 4 字节，若为 "RIFF" 则按 WAV 加载（采样率从 WAV 头读取），
 * 否则按 16-bit 单声道裸 PCM 加载（需通过 force_sample_rate 指定采样率）。
 */
audio_t* audio_load(const char* filename, int force_sample_rate)
{
    if (!filename) return NULL;

    audio_t* a = (audio_t*)calloc(1, sizeof(audio_t));
    if (!a) return NULL;

    /* 打开文件检查魔数 */
    FILE* fp = fopen(filename, "rb");
    if (!fp) {
        fprintf(stderr, "ERROR: cannot open file: %s\n", filename);
        free(a);
        return NULL;
    }

    char magic[4] = {0};
    if (fread(magic, 1, 4, fp) < 4) {
        fclose(fp);
        free(a);
        return NULL;
    }
    fclose(fp);

    if (memcmp(magic, "RIFF", 4) == 0) {
        /* ---------- WAV 格式 ---------- */
        WAVHeader hdr;
        a->data = read_wav(filename, &hdr, &a->num_samples);
        if (!a->data) {
            free(a);
            return NULL;
        }
        a->sample_rate     = hdr.sampleRate;
        a->bits_per_sample = hdr.bitsPerSample;
        a->num_channels    = hdr.numChannels;

        printf("  [WAV] %s: %u Hz, %u-bit, %u ch, %u samples (%.2f s)\n",
               filename, a->sample_rate, a->bits_per_sample, a->num_channels,
               a->num_samples, (double)a->num_samples / a->sample_rate);
    } else {
        /* ---------- 裸 PCM 格式 ---------- */
        if (force_sample_rate <= 0) {
            fprintf(stderr, "ERROR: PCM file '%s' requires specifying sample_rate.\n", filename);
            fprintf(stderr, "       Usage: %s input1.pcm input2.pcm output.wav sample_rate\n", "beamforming.exe");
            free(a);
            return NULL;
        }
        a->data = read_pcm(filename, &a->num_samples);
        if (!a->data) {
            free(a);
            return NULL;
        }
        a->sample_rate     = (uint32_t)force_sample_rate;
        a->bits_per_sample = 16;
        a->num_channels    = 1;

        printf("  [PCM] %s: %u Hz, 16-bit, mono, %u samples (%.2f s)\n",
               filename, a->sample_rate, a->num_samples,
               (double)a->num_samples / a->sample_rate);
    }

    return a;
}

void audio_free(audio_t* a)
{
    if (!a) return;
    free(a->data);
    free(a);
}

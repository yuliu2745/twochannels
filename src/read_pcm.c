#include "../include/read_pcm.h"

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

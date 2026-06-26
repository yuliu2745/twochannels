#include "../include/setting.h"

// 读取WAV文件，返回样本数据（int16_t*），并通过参数返回样本数和头信息
int16_t* read_wav(const char* filename, WAVHeader* header, uint32_t* numSamples) {
    FILE* fp = fopen(filename, "rb");
    if (!fp) {
        fprintf(stderr, "Cannot open file %s\n", filename);
        return NULL;
    }

    // 读取头
    if (fread(header, sizeof(WAVHeader), 1, fp) != 1) {
        fprintf(stderr, "Failed to read WAV header: %s\n", filename);
        fclose(fp);
        return NULL;
    }

    // 验证是否为PCM格式
    if (header->audioFormat != 1) {
        fprintf(stderr, "Unsupported audio format (PCM only): %s\n", filename);
        fclose(fp);
        return NULL;
    }
    if (header->bitsPerSample != 16) {
        fprintf(stderr, "Unsupported bit depth (16-bit only): %s\n", filename);
        fclose(fp);
        return NULL;
    }
    // 支持单声道和立体声
    if (header->numChannels != 1 && header->numChannels != 2) {
        fprintf(stderr, "Unsupported channels: %d (only mono and stereo supported)\n", header->numChannels);
        fclose(fp);
        return NULL;
    }

    // 计算样本数
    *numSamples = header->subchunk2Size / (header->bitsPerSample / 8);
    int16_t* data = (int16_t*)malloc(header->subchunk2Size);
    if (!data) {
        fprintf(stderr, "Memory allocation failed\n");
        fclose(fp);
        return NULL;
    }

    // 读取数据
    if (fread(data, header->subchunk2Size, 1, fp) != 1) {
        fprintf(stderr, "Failed to read audio data: %s\n", filename);
        free(data);
        fclose(fp);
        return NULL;
    }

    fclose(fp);
    return data;
}

// 写入WAV文件
int write_wav(const char* filename, const WAVHeader* templateHeader, const int16_t* data, uint32_t numSamples) {
    FILE* fp = fopen(filename, "wb");
    if (!fp) {
        fprintf(stderr, "Cannot create file %s\n", filename);
        return 0;
    }

    // 创建输出头，基于模板但修改数据大小
    WAVHeader outHeader = *templateHeader;
    outHeader.subchunk2Size = numSamples * (outHeader.bitsPerSample / 8);
    outHeader.byteRate = outHeader.sampleRate * outHeader.numChannels * (outHeader.bitsPerSample / 8);
    outHeader.blockAlign = outHeader.numChannels * (outHeader.bitsPerSample / 8);
    outHeader.chunkSize = 36 + outHeader.subchunk2Size; // 文件大小-8

    // 写入头
    if (fwrite(&outHeader, sizeof(WAVHeader), 1, fp) != 1) {
        fprintf(stderr, "Failed to write WAV header\n");
        fclose(fp);
        return 0;
    }

    // 写入数据
    if (fwrite(data, outHeader.subchunk2Size, 1, fp) != 1) {
        fprintf(stderr, "Failed to write audio data\n");
        fclose(fp);
        return 0;
    }

    fclose(fp);
    return 1;
}
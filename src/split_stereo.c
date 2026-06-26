#include "../include/setting.h"
#include "../include/file_utils.h"

// 分离立体声文件为两个单通道文件
int split_stereo(const char* stereoFile, const char* leftFile, const char* rightFile) {
    FILE* fp = fopen(stereoFile, "rb");
    if (!fp) {
        fprintf(stderr, "Cannot open stereo file: %s\n", stereoFile);
        return 0;
    }

    // 读取WAV头
    WAVHeader header;
    if (fread(&header, sizeof(WAVHeader), 1, fp) != 1) {
        fprintf(stderr, "Failed to read WAV header: %s\n", stereoFile);
        fclose(fp);
        return 0;
    }

    // 验证是否为立体声
    if (header.numChannels != 2) {
        fprintf(stderr, "Not a stereo file (channels: %d)\n", header.numChannels);
        fclose(fp);
        return 0;
    }

    if (header.bitsPerSample != 16) {
        fprintf(stderr, "Only 16-bit audio supported\n");
        fclose(fp);
        return 0;
    }

    // 计算样本数
    uint32_t totalSamples = header.subchunk2Size / (header.bitsPerSample / 8);
    uint32_t samplesPerChannel = totalSamples / header.numChannels;

    // 读取音频数据
    int16_t* stereoData = (int16_t*)malloc(header.subchunk2Size);
    if (!stereoData) {
        fprintf(stderr, "Memory allocation failed\n");
        fclose(fp);
        return 0;
    }

    if (fread(stereoData, header.subchunk2Size, 1, fp) != 1) {
        fprintf(stderr, "Failed to read audio data\n");
        free(stereoData);
        fclose(fp);
        return 0;
    }
    fclose(fp);

    // 分离左右声道
    int16_t* leftData = (int16_t*)malloc(samplesPerChannel * sizeof(int16_t));
    int16_t* rightData = (int16_t*)malloc(samplesPerChannel * sizeof(int16_t));
    
    if (!leftData || !rightData) {
        fprintf(stderr, "Memory allocation failed\n");
        free(stereoData);
        if (leftData) free(leftData);
        if (rightData) free(rightData);
        return 0;
    }

    for (uint32_t i = 0; i < samplesPerChannel; i++) {
        leftData[i] = stereoData[2 * i];     // 左声道
        rightData[i] = stereoData[2 * i + 1]; // 右声道
    }

    // 创建单声道WAV头
    WAVHeader monoHeader = header;
    monoHeader.numChannels = 1;
    monoHeader.subchunk2Size = samplesPerChannel * (monoHeader.bitsPerSample / 8);
    monoHeader.byteRate = monoHeader.sampleRate * monoHeader.numChannels * (monoHeader.bitsPerSample / 8);
    monoHeader.blockAlign = monoHeader.numChannels * (monoHeader.bitsPerSample / 8);
    monoHeader.chunkSize = 36 + monoHeader.subchunk2Size;

    // 写入左声道文件
    FILE* leftFp = fopen(leftFile, "wb");
    if (!leftFp) {
        fprintf(stderr, "Cannot create left channel file: %s\n", leftFile);
        free(stereoData);
        free(leftData);
        free(rightData);
        return 0;
    }

    fwrite(&monoHeader, sizeof(WAVHeader), 1, leftFp);
    fwrite(leftData, monoHeader.subchunk2Size, 1, leftFp);
    fclose(leftFp);

    // 写入右声道文件
    FILE* rightFp = fopen(rightFile, "wb");
    if (!rightFp) {
        fprintf(stderr, "Cannot create right channel file: %s\n", rightFile);
        free(stereoData);
        free(leftData);
        free(rightData);
        return 0;
    }

    fwrite(&monoHeader, sizeof(WAVHeader), 1, rightFp);
    fwrite(rightData, monoHeader.subchunk2Size, 1, rightFp);
    fclose(rightFp);

    // 清理内存
    free(stereoData);
    free(leftData);
    free(rightData);

    printf("Successfully split stereo file:\n");
    printf("  Left channel: %s (%u samples)\n", leftFile, samplesPerChannel);
    printf("  Right channel: %s (%u samples)\n", rightFile, samplesPerChannel);
    printf("  Sample rate: %u Hz\n", header.sampleRate);

    return 1;
}

// Simple main function for standalone execution
int main(int argc, char* argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s input_stereo.wav left_output.wav right_output.wav\n", argv[0]);
        return 1;
    }
    return split_stereo(argv[1], argv[2], argv[3]);
}

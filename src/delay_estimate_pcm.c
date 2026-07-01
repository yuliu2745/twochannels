/**
 * @file delay_estimate_pcm.c
 * @brief 裸 PCM 文件延迟估计工具
 * @details 读取两个 16-bit 单声道小端 PCM 文件，
 *          使用时域互相关或 FFT-PHAT 算法估计通道间延迟，
 *          并可选择性地生成对齐后的 WAV 输出文件。
 *
 * Usage:
 *   delay_estimate_pcm.exe pcm_file1 pcm_file2 sample_rate [output.wav]
 *
 * Examples:
 *   delay_estimate_pcm.exe mic1.pcm mic2.pcm 16000
 *   delay_estimate_pcm.exe mic1.pcm mic2.pcm 16000 aligned_output.wav
 *   echo 2 | delay_estimate_pcm.exe mic1.pcm mic2.pcm 48000 aligned.wav
 */

#include "../include/setting.h"
#include "../include/read_pcm.h"
#include "../include/delay_and_sum.h"
#include "../include/fft_path.h"
#include "../include/file_utils.h"
#include "../include/readwav.h"

int main(int argc, char* argv[]) {
    if (argc < 4) {
        fprintf(stderr, "Usage: %s pcm_file1 pcm_file2 sample_rate [output.wav]\n", argv[0]);
        fprintf(stderr, "  pcm_file1, pcm_file2 : 16-bit mono little-endian PCM files\n");
        fprintf(stderr, "  sample_rate           : sampling rate in Hz (e.g. 16000, 44100, 48000)\n");
        fprintf(stderr, "  output.wav (optional) : output beamformed WAV file\n");
        fprintf(stderr, "\nExamples:\n");
        fprintf(stderr, "  %s mic1.pcm mic2.pcm 16000\n", argv[0]);
        fprintf(stderr, "  echo 2 | %s mic1.pcm mic2.pcm 48000 aligned.wav\n", argv[0]);
        return 1;
    }

    const char* file1 = argv[1];
    const char* file2 = argv[2];
    uint32_t sampleRate = (uint32_t)atoi(argv[3]);

    if (sampleRate == 0) {
        fprintf(stderr, "ERROR: Invalid sample rate '%s'. Must be a positive integer.\n", argv[3]);
        return 1;
    }

    printf("========================================\n");
    printf("  PCM Delay Estimator\n");
    printf("========================================\n");
    printf("Input 1: %s\n", file1);
    printf("Input 2: %s\n", file2);
    printf("Sample rate: %u Hz\n", sampleRate);
    printf("Format: 16-bit mono little-endian PCM\n");
    printf("----------------------------------------\n");

    // 读取两个 PCM 文件
    uint32_t len1, len2;
    int16_t* data1 = read_pcm(file1, &len1);
    if (!data1) return 1;

    int16_t* data2 = read_pcm(file2, &len2);
    if (!data2) {
        free(data1);
        return 1;
    }

    double dur1 = (double)len1 / sampleRate;
    double dur2 = (double)len2 / sampleRate;
    printf("  Duration 1: %.3f s\n", dur1);
    printf("  Duration 2: %.3f s\n", dur2);
    printf("----------------------------------------\n");

    // 选择延迟估计算法
    printf("Delay estimation method:\n");
    printf("  1. Time-domain cross-correlation (fast, good for clean signals)\n");
    printf("  2. FFT-PHAT (robust for noisy signals)\n");
    printf("Enter choice (1 or 2): ");

    int choice;
    if (scanf("%d", &choice) != 1) {
        choice = 1;
    }
    // 清空输入缓冲区（处理管道输入）
    int c;
    while ((c = getchar()) != '\n' && c != EOF);

    int estimated_delay = 0;
    float confidence = 0.0f;
    int method_used = 0; // 0=time, 1=fft-phat
    int peak_num = 0;    // 用于 FFT-PHAT 的结果计数

    if (choice == 2) {
        printf("\nUsing FFT-PHAT delay estimation...\n");

        // int16_t -> float 转换
        float* float_data1 = (float*)malloc(len1 * sizeof(float));
        float* float_data2 = (float*)malloc(len2 * sizeof(float));
        if (!float_data1 || !float_data2) {
            fprintf(stderr, "Memory allocation failed for float conversion\n");
            free(data1); free(data2);
            free(float_data1); free(float_data2);
            return 1;
        }

        for (uint32_t i = 0; i < len1; i++) {
            float_data1[i] = (float)data1[i] / 32768.0f;
        }
        for (uint32_t i = 0; i < len2; i++) {
            float_data2[i] = (float)data2[i] / 32768.0f;
        }

        // FFT 参数
        int fft_size = 512;
        int window = (int)((len1 < (uint32_t)fft_size) ? len1 : (uint32_t)fft_size);
        int margin = 200;
        peak_num = 1;

        int* delays = (int*)malloc(peak_num * sizeof(int));
        float* peak_values = (float*)malloc(peak_num * sizeof(float));
        if (!delays || !peak_values) {
            fprintf(stderr, "Memory allocation failed for FFT results\n");
            free(data1); free(data2);
            free(float_data1); free(float_data2);
            free(delays); free(peak_values);
            return 1;
        }

        FFT_Real_Gcc_Path(delays, peak_values, &peak_num,
                                       float_data2, float_data1,
                                       margin, window, fft_size);

        if (peak_num > 0) {
            estimated_delay = delays[0];
            confidence = peak_values[0];
            method_used = 1;
            printf("  Estimated delay: %d samples (confidence: %.6f)\n",
                   estimated_delay, confidence);
        } else {
            printf("  FFT-PHAT failed, falling back to time-domain method\n");
            estimated_delay = estimate_delay(data1, len1, data2, len2, 5000);
        }

        free(float_data1);
        free(float_data2);
        free(delays);
        free(peak_values);
    }

    if (choice != 2 || (choice == 2 && peak_num <= 0)) {
        printf("\nUsing time-domain delay estimation...\n");
        int max_delay = 5000;

        // 限制搜索范围不超过信号长度
        if ((int)len1 - 1 < max_delay) max_delay = (int)len1 - 1;
        if ((int)len2 - 1 < max_delay) max_delay = (int)len2 - 1;
        if (max_delay < 0) max_delay = 0;

        estimated_delay = estimate_delay(data1, len1, data2, len2, max_delay);
        method_used = 0;
        printf("  Estimated delay: %d samples\n", estimated_delay);
    }

    // 计算结果
    float delay_ms = (float)estimated_delay / sampleRate * 1000.0f;

    printf("\n========================================\n");
    printf("  ** DELAY ESTIMATION RESULT **\n");
    printf("========================================\n");
    printf("  Method : %s\n", method_used ? "FFT-PHAT" : "Time-domain cross-correlation");
    printf("  Delay  : %+d samples\n", estimated_delay);
    printf("  Delay  : %+.4f ms (at %u Hz)\n", delay_ms, sampleRate);
    if (method_used) {
        printf("  Confidence : %.6f\n", confidence);
    }
    printf("----------------------------------------\n");
    if (estimated_delay >= 0) {
        printf("  --> file2 LAGS file1 by %d samples (%.4f ms)\n",
               estimated_delay, delay_ms);
    } else {
        printf("  --> file2 LEADS file1 by %d samples (%.4f ms)\n",
               -estimated_delay, -delay_ms);
    }
    printf("========================================\n");

    // 可选：生成对齐后的 WAV 输出文件
    if (argc >= 5) {
        const char* outFile = argv[4];
        printf("\nGenerating beamformed output WAV: %s\n", outFile);

        ensure_audio_dir();
        char* outPath = build_audio_path(outFile);

        // 分配两个通道的延迟
        int delay1, delay2;
        if (estimated_delay >= 0) {
            delay1 = 0;
            delay2 = estimated_delay;
        } else {
            delay1 = -estimated_delay;
            delay2 = 0;
        }

        // 执行延迟求和
        uint32_t outLen;
        int16_t* outData = delay_sum(data1, len1, delay1,
                                    data2, len2, delay2, &outLen);
        if (!outData) {
            fprintf(stderr, "Error: delay_sum failed\n");
            free(data1); free(data2); free(outPath);
            return 1;
        }

        // 构建 WAV 头并写入
        WAVHeader h;
        memset(&h, 0, sizeof(WAVHeader));
        memcpy(h.chunkID, "RIFF", 4);
        memcpy(h.format, "WAVE", 4);
        memcpy(h.subchunk1ID, "fmt ", 4);
        h.subchunk1Size = 16;
        h.audioFormat = 1;
        h.numChannels = 1;
        h.sampleRate = sampleRate;
        h.bitsPerSample = 16;
        h.byteRate = sampleRate * 1 * (h.bitsPerSample / 8);
        h.blockAlign = 1 * (h.bitsPerSample / 8);
        h.subchunk2Size = outLen * (h.bitsPerSample / 8);
        h.chunkSize = 36 + h.subchunk2Size;
        memcpy(h.subchunk2ID, "data", 4);

        if (write_wav(outPath, &h, outData, outLen)) {
            printf("  Output: %s\n", outPath);
            printf("  Duration: %.3f seconds\n", (float)outLen / sampleRate);
        } else {
            fprintf(stderr, "Error: failed to write WAV file\n");
        }

        free(outData);
        free(outPath);
    }

    // 清理
    free(data1);
    free(data2);

    printf("\nDone.\n");
    return 0;
}

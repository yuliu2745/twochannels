#include "../include/setting.h"
#include "../include/readwav.h"
#include "../include/delay_and_sum.h"
#include "../include/file_utils.h"
#include "../include/fft_path.h"
#include <stdio.h>
#include <stdlib.h>

int do_work_fft(int argc, char* argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s input1.wav input2.wav output.wav\n", argv[0]);
        fprintf(stderr, "       Uses FFT-PHAT delay estimation\n");
        return 1;
    }

    const char* file1 = argv[1];
    const char* file2 = argv[2];
    const char* outFile = argv[3];
    
    printf("=== FFT-PHAT Beamforming ===\n");
    
    // 确保音频目录存在
    ensure_audio_dir();
    
    // 构建输出文件的完整路径
    char* outPath = build_audio_path(outFile);
    
    // 读取第一个文件
    WAVHeader h1;
    uint32_t len1;
    int16_t* data1 = read_wav(file1, &h1, &len1);
    if (!data1) return 1;

    // 读取第二个文件
    WAVHeader h2;
    uint32_t len2;
    int16_t* data2 = read_wav(file2, &h2, &len2);
    if (!data2) {
        free(data1);
        return 1;
    }

    // 验证格式兼容
    if (h1.sampleRate != h2.sampleRate ||
        h1.bitsPerSample != h2.bitsPerSample ||
        h1.numChannels != h2.numChannels) {
        fprintf(stderr, "WAV files format incompatible\n");
        free(data1);
        free(data2);
        return 1;
    }

    printf("Audio files loaded:\n");
    printf("  Sample rate: %d Hz\n", h1.sampleRate);
    printf("  Duration: %.2f seconds\n", (float)len1 / h1.sampleRate);
    printf("  Samples: %u\n", len1);

    // 使用FFT-PHAT延迟估计
    printf("\nUsing FFT-PHAT delay estimation...\n");
    
    // 转换为float格式
    float* float_data1 = (float*)malloc(len1 * sizeof(float));
    float* float_data2 = (float*)malloc(len2 * sizeof(float));
    
    if (!float_data1 || !float_data2) {
        printf("ERROR: Memory allocation failed for float arrays\n");
        free(data1);
        free(data2);
        return 1;
    }
    
    for (uint32_t i = 0; i < len1; i++) {
        float_data1[i] = (float)data1[i] / 32768.0f;
        float_data2[i] = (float)data2[i] / 32768.0f;
    }
    
    // FFT参数设置
    int fft_size = 512;
    int window = (len1 < fft_size) ? len1 : fft_size;
    int margin = 200;
    int peak_num = 1;
    
    printf("FFT parameters:\n");
    printf("  FFT size: %d\n", fft_size);
    printf("  Window size: %d\n", window);
    printf("  Search range: +/- %d samples\n", margin);
    printf("  Peak count: %d\n", peak_num);
    
    // 分配结果数组
    int* delays = (int*)malloc(peak_num * sizeof(int));
    float* peak_values = (float*)malloc(peak_num * sizeof(float));
    
    if (!delays || !peak_values) {
        printf("ERROR: Memory allocation failed for results\n");
        free(data1);
        free(data2);
        free(float_data1);
        free(float_data2);
        return 1;
    }
    
    printf("Starting FFT-PHAT analysis...\n");
    
    // 执行FFT-PHAT
    float result = FFT_Real_Gcc_Path(delays, peak_values, &peak_num,
                                   float_data2, float_data1,
                                   margin, window, fft_size);
    
    printf("FFT-PHAT analysis completed\n");
    
    if (peak_num > 0 && result > 0) {
        int estimated_delay = delays[0];
        float delay_ms = (float)estimated_delay / h1.sampleRate * 1000.0f;
        
        printf("\nFFT-PHAT Results:\n");
        printf("  Estimated delay: %+d samples (%.2f ms)\n", estimated_delay, delay_ms);
        printf("  Confidence: %.6f\n", peak_values[0]);
        
        // 设置延迟值
        int delay1, delay2;
        if (estimated_delay >= 0) {
            delay1 = 0;
            delay2 = estimated_delay;
        } else {
            delay1 = -estimated_delay;
            delay2 = 0;
        }
        
        printf("  Applied delays: delay1=%d, delay2=%d\n", delay1, delay2);

        // 延迟求和
        printf("Performing delay-and-sum...\n");
        uint32_t outLen;
        int16_t* outData = delay_sum(data1, len1, (float)delay1, data2, len2, (float)delay2, &outLen);
        if (!outData) {
            fprintf(stderr, "Error in delay_sum operation\n");
            free(data1);
            free(data2);
            free(float_data1);
            free(float_data2);
            free(delays);
            free(peak_values);
            free(outPath);
            return 1;
        }

        // 写入输出文件
        printf("Writing output file...\n");
        if (!write_wav(outPath, &h1, outData, outLen)) {
            fprintf(stderr, "Error writing output file\n");
            free(data1);
            free(data2);
            free(outData);
            free(float_data1);
            free(float_data2);
            free(delays);
            free(peak_values);
            free(outPath);
            return 1;
        }

        printf("\nSuccess!\n");
        printf("Output file: %s\n", outPath);
        printf("Output samples: %u\n", outLen);
        printf("Output duration: %.2f seconds\n", (float)outLen / h1.sampleRate);

        free(outData);
    } else {
        fprintf(stderr, "FFT-PHAT failed to estimate delay!\n");
        fprintf(stderr, "Result: %.6f, Peak count: %d\n", result, peak_num);
    }

    // 清理内存
    free(data1);
    free(data2);
    free(float_data1);
    free(float_data2);
    free(delays);
    free(peak_values);
    free(outPath);
    
    printf("=== Process Complete ===\n");
    return 0;
}

int main(int argc, char* argv[]) {
    return do_work_fft(argc, argv);
}

#include "../include/merge_audio.h"
#include "../include/readwav.h"
#include "../include/file_utils.h"
#include <math.h>

// 重采样函数（简单的线性插值）
int16_t* resample(const int16_t* input, uint32_t inputLength, uint32_t inputRate, uint32_t outputRate, uint32_t* outputLength) {
    if (inputRate == outputRate) {
        // 不需要重采样
        *outputLength = inputLength;
        int16_t* output = (int16_t*)malloc(inputLength * sizeof(int16_t));
        if (output) {
            memcpy(output, input, inputLength * sizeof(int16_t));
        }
        return output;
    }

    *outputLength = (uint32_t)((double)inputLength * outputRate / inputRate);
    int16_t* output = (int16_t*)malloc(*outputLength * sizeof(int16_t));
    if (!output) return NULL;

    for (uint32_t i = 0; i < *outputLength; i++) {
        double srcPos = (double)i * inputRate / outputRate;
        uint32_t srcIndex = (uint32_t)srcPos;
        double frac = srcPos - srcIndex;

        if (srcIndex >= inputLength - 1) {
            output[i] = input[inputLength - 1];
        } else {
            // 线性插值
            double val = input[srcIndex] * (1.0 - frac) + input[srcIndex + 1] * frac;
            output[i] = (int16_t)val;
        }
    }

    return output;
}

// 应用延迟到信号
int16_t* apply_delay(const int16_t* input, uint32_t length, int delay_samples, uint32_t* outputLength) {
    *outputLength = length;
    int16_t* output = (int16_t*)calloc(length, sizeof(int16_t));
    if (!output) return NULL;

    if (delay_samples >= 0) {
        // 正延迟：开头补零
        for (int i = 0; i < length - delay_samples; i++) {
            if (i + delay_samples < length) {
                output[i + delay_samples] = input[i];
            }
        }
    } else {
        // 负延迟：开头截断
        int abs_delay = -delay_samples;
        for (int i = 0; i < length - abs_delay; i++) {
            output[i] = input[i + abs_delay];
        }
    }

    return output;
}

// 计算TDOA延迟样本数
int calculate_tdoa_delay(float mic_distance, float sound_speed, float angle_degrees, uint32_t sample_rate) {
    // 将角度转换为弧度
    double angle_rad = angle_degrees * M_PI / 180.0;
    
    // 计算时间差：tdoa = d * sin(θ) / c
    double time_diff = mic_distance * sin(angle_rad) / sound_speed;
    
    // 转换为样本数
    int delay_samples = (int)(time_diff * sample_rate);
    
    return delay_samples;
}

// 合并两个音频信号
int16_t* merge_signals(const int16_t* signal1, const int16_t* signal2, uint32_t length) {
    int16_t* merged = (int16_t*)malloc(length * sizeof(int16_t));
    if (!merged) return NULL;

    for (uint32_t i = 0; i < length; i++) {
        int32_t sum = (int32_t)signal1[i] + (int32_t)signal2[i];
        
        // 防止溢出
        if (sum > 32767) sum = 32767;
        if (sum < -32768) sum = -32768;
        
        merged[i] = (int16_t)sum;
    }

    return merged;
}

// 生成测试用的双通道音频
int generate_stereo_with_tdoa(const char* voiceFile, const char* sineFile,
                              const char* outputFile, float mic_distance,
                              float sound_speed, float angle_voice, float angle_noise) {
    
    // 读取人声文件
    WAVHeader voiceHeader;
    uint32_t voiceLength;
    int16_t* voiceData = read_wav(voiceFile, &voiceHeader, &voiceLength);
    if (!voiceData) {
        fprintf(stderr, "Failed to read voice file: %s\n", voiceFile);
        return 0;
    }

    // 读取正弦波文件
    WAVHeader sineHeader;
    uint32_t sineLength;
    int16_t* sineData = read_wav(sineFile, &sineHeader, &sineLength);
    if (!sineData) {
        fprintf(stderr, "Failed to read sine wave file: %s\n", sineFile);
        free(voiceData);
        return 0;
    }

    // 如果是立体声，转换为单声道（取平均值）
    int16_t* originalSineData = sineData;
    uint32_t originalSineLength = sineLength;
    
    if (sineHeader.numChannels == 2) {
        printf("Converting stereo to mono (averaging channels)\n");
        uint32_t samplesPerChannel = sineLength / 2;
        sineLength = samplesPerChannel;
        sineData = (int16_t*)malloc(sineLength * sizeof(int16_t));
        if (!sineData) {
            fprintf(stderr, "Memory allocation failed for mono conversion\n");
            free(originalSineData);
            free(voiceData);
            return 0;
        }
        
        for (uint32_t i = 0; i < sineLength; i++) {
            int32_t left = originalSineData[2 * i];
            int32_t right = originalSineData[2 * i + 1];
            int32_t mono = (left + right) / 2;
            
            // 防止溢出
            if (mono > 32767) mono = 32767;
            if (mono < -32768) mono = -32768;
            
            sineData[i] = (int16_t)mono;
        }
        
        free(originalSineData);
        sineHeader.numChannels = 1;
        sineHeader.subchunk2Size = sineLength * (sineHeader.bitsPerSample / 8);
        sineHeader.byteRate = sineHeader.sampleRate * sineHeader.numChannels * (sineHeader.bitsPerSample / 8);
        sineHeader.blockAlign = sineHeader.numChannels * (sineHeader.bitsPerSample / 8);
    }

    // 检查采样率是否一致
    uint32_t targetSampleRate = voiceHeader.sampleRate;
    int16_t* resampledSine = NULL;
    uint32_t resampledSineLength = sineLength;

    if (voiceHeader.sampleRate != sineHeader.sampleRate) {
        printf("Resampling sine wave from %u Hz to %u Hz\n", 
               sineHeader.sampleRate, voiceHeader.sampleRate);
        resampledSine = resample(sineData, sineLength, sineHeader.sampleRate, 
                                voiceHeader.sampleRate, &resampledSineLength);
        if (!resampledSine) {
            fprintf(stderr, "Failed to resample sine wave\n");
            free(voiceData);
            free(sineData);
            return 0;
        }
    } else {
        resampledSine = sineData;
    }

    // 使用较长的文件长度，循环使用较短的文件
    uint32_t maxLength = (voiceLength > resampledSineLength) ? voiceLength : resampledSineLength;
    printf("Using %u samples for processing\n", maxLength);
    
    // 如果文件长度不同，循环使用较短的文件
    int16_t* extendedVoice = voiceData;
    int16_t* extendedSine = resampledSine;
    uint32_t voiceLengthExtended = voiceLength;
    uint32_t sineLengthExtended = resampledSineLength;
    
    if (voiceLength < resampledSineLength) {
        // 扩展voice文件
        extendedVoice = (int16_t*)malloc(maxLength * sizeof(int16_t));
        if (!extendedVoice) {
            fprintf(stderr, "Failed to allocate memory for voice extension\n");
            free(voiceData);
            if (resampledSine != sineData) free(resampledSine);
            else free(sineData);
            return 0;
        }
        
        for (uint32_t i = 0; i < maxLength; i++) {
            extendedVoice[i] = voiceData[i % voiceLength];
        }
        voiceLengthExtended = maxLength;
        printf("Extended voice from %u to %u samples (looping)\n", voiceLength, maxLength);
    } else if (resampledSineLength < voiceLength) {
        // 扩展sine文件
        extendedSine = (int16_t*)malloc(maxLength * sizeof(int16_t));
        if (!extendedSine) {
            fprintf(stderr, "Failed to allocate memory for sine extension\n");
            free(voiceData);
            if (resampledSine != sineData) free(resampledSine);
            else free(sineData);
            return 0;
        }
        
        for (uint32_t i = 0; i < maxLength; i++) {
            extendedSine[i] = resampledSine[i % resampledSineLength];
        }
        sineLengthExtended = maxLength;
        printf("Extended sine from %u to %u samples (looping)\n", resampledSineLength, maxLength);
    }

    // 计算TDOA延迟
    int voiceDelay = calculate_tdoa_delay(mic_distance, sound_speed, angle_voice, targetSampleRate);
    int noiseDelay = calculate_tdoa_delay(mic_distance, sound_speed, angle_noise, targetSampleRate);
    
    printf("TDOA delays - Voice: %d samples, Noise: %d samples\n", voiceDelay, noiseDelay);
    printf("Angles - Voice: %.1f°, Noise: %.1f°\n", angle_voice, angle_noise);

    // 应用延迟到人声和正弦波
    uint32_t delayedVoiceLength, delayedNoiseLength;
    int16_t* delayedVoice = apply_delay(extendedVoice, maxLength, voiceDelay, &delayedVoiceLength);
    int16_t* delayedNoise = apply_delay(extendedSine, maxLength, noiseDelay, &delayedNoiseLength);
    
    if (!delayedVoice || !delayedNoise) {
        fprintf(stderr, "Failed to apply delays\n");
        free(voiceData);
        if (resampledSine != sineData) free(resampledSine);
        else free(sineData);
        if (delayedVoice) free(delayedVoice);
        if (delayedNoise) free(delayedNoise);
        return 0;
    }

    // 麦克风0：直接叠加（参考通道）
    int16_t* mic0 = merge_signals(extendedVoice, extendedSine, maxLength);
    if (!mic0) {
        fprintf(stderr, "Failed to create mic0 signal\n");
        free(delayedVoice);
        free(delayedNoise);
        free(voiceData);
        if (resampledSine != sineData) free(resampledSine);
        else free(sineData);
        if (extendedVoice != voiceData) free(extendedVoice);
        if (extendedSine != resampledSine) free(extendedSine);
        return 0;
    }

    // 麦克风1：延迟后叠加
    int16_t* mic1 = merge_signals(delayedVoice, delayedNoise, maxLength);
    if (!mic1) {
        fprintf(stderr, "Failed to create mic1 signal\n");
        free(mic0);
        free(delayedVoice);
        free(delayedNoise);
        free(voiceData);
        if (resampledSine != sineData) free(resampledSine);
        else free(sineData);
        if (extendedVoice != voiceData) free(extendedVoice);
        if (extendedSine != resampledSine) free(extendedSine);
        return 0;
    }

    // 创建立体声数据
    uint32_t stereoLength = maxLength;
    int16_t* stereoData = (int16_t*)malloc(stereoLength * 2 * sizeof(int16_t));
    if (!stereoData) {
        fprintf(stderr, "Failed to allocate stereo buffer\n");
        free(mic0);
        free(mic1);
        free(delayedVoice);
        free(delayedNoise);
        free(voiceData);
        if (resampledSine != sineData) free(resampledSine);
        else free(sineData);
        return 0;
    }

    // 交错左右声道
    for (uint32_t i = 0; i < stereoLength; i++) {
        stereoData[2 * i] = mic0[i];     // 左声道
        stereoData[2 * i + 1] = mic1[i];   // 右声道
    }

    // 创建立体声WAV头
    WAVHeader stereoHeader = voiceHeader;
    stereoHeader.numChannels = 2;
    stereoHeader.subchunk2Size = stereoLength * 2 * (stereoHeader.bitsPerSample / 8);
    stereoHeader.byteRate = stereoHeader.sampleRate * stereoHeader.numChannels * (stereoHeader.bitsPerSample / 8);
    stereoHeader.blockAlign = stereoHeader.numChannels * (stereoHeader.bitsPerSample / 8);
    stereoHeader.chunkSize = 36 + stereoHeader.subchunk2Size;

    // 写入立体声文件
    if (!write_wav(outputFile, &stereoHeader, stereoData, stereoLength * 2)) {
        fprintf(stderr, "Failed to write stereo file: %s\n", outputFile);
        free(stereoData);
        free(mic0);
        free(mic1);
        free(delayedVoice);
        free(delayedNoise);
        free(voiceData);
        if (resampledSine != sineData) free(resampledSine);
        else free(sineData);
        return 0;
    }

    printf("Successfully generated stereo file: %s\n", outputFile);
    printf("Duration: %.2f seconds\n", (double)stereoLength / stereoHeader.sampleRate);
    printf("Sample rate: %u Hz\n", stereoHeader.sampleRate);
    printf("Channels: %d\n", stereoHeader.numChannels);

    // 清理内存
    free(stereoData);
    free(mic0);
    free(mic1);
    free(delayedVoice);
    free(delayedNoise);
    free(voiceData);
    if (resampledSine != sineData) free(resampledSine);
    else free(sineData);
    
    // 清理扩展的音频数据
    if (extendedVoice != voiceData) free(extendedVoice);
    if (extendedSine != resampledSine) free(extendedSine);

    return 1;
}

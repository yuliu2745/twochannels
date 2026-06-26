/**
 * @file main.c
 * @brief 双通道波束成形主程序
 * @details 提供音频延迟估计和波束形成功能，支持时域和频域两种算法
 *          同时提供API接口供其他模块调用
 */

#include "include/setting.h"      // 全局设置和数据结构定义
#include "include/readwav.h"      // WAV文件读写功能
#include "include/delay_and_sum.h" // 延迟估计和求和算法
#include "include/file_utils.h"    // 文件路径处理工具
#include "include/fft_path.h"     // FFT-PHAT频域延迟估计算法

/**
 * @brief 波束成形主工作函数
 * @param argc 命令行参数个数
 * @param argv 命令行参数数组
 * @return 程序执行状态，0表示成功，非0表示失败
 * 
 * 使用方法：
 * - 自动延迟估计：program input1.wav input2.wav output.wav
 * - 手动指定延迟：program input1.wav input2.wav output.wav delay1 delay2
 * 
 * 功能说明：
 * 1. 读取两个WAV音频文件
 * 2. 自动或手动计算延迟
 * 3. 应用延迟求和算法进行波束形成
 * 4. 输出增强后的音频文件
 */
int do_work(int argc, char* argv[]) {
    // 检查命令行参数合法性
    // 支持4个参数（自动延迟）或6个参数（手动延迟）
    if (argc != 4 && argc != 6) {
        fprintf(stderr, "Usage: %s input1.wav input2.wav output.wav [delay1 delay2]\n", argv[0]);
        fprintf(stderr, "       Auto-calculate delay if not specified\n");
        return 1;
    }

    // 解析命令行参数
    const char* file1 = argv[1];  // 第一个输入音频文件
    const char* file2 = argv[2];  // 第二个输入音频文件
    const char* outFile = argv[3]; // 输出音频文件
    
    // 确保音频输出目录存在，如果不存在则创建
    ensure_audio_dir();
    
    // 构建输出文件的完整路径（包含audio_files目录）
    char* outPath = build_audio_path(outFile);
    int delay1, delay2;           // 两个通道的延迟值（样本数）
    int auto_delay = 0;           // 标志位：1表示自动计算延迟，0表示手动指定

    // 检查是否手动指定延迟参数
    // 如果有6个参数，则最后两个参数为手动指定的延迟值
    if (argc == 6) {
        delay1 = atoi(argv[4]); // 第一个通道的延迟（样本数）
        delay2 = atoi(argv[5]); // 第二个通道的延迟（样本数）
        printf("Using manual delay: delay1=%d, delay2=%d\n", delay1, delay2);
    } else {
        auto_delay = 1; // 设置标志位，表示需要自动计算延迟
        printf("Auto-calculating delay...\n");
    }

    // 读取第一个音频文件
    WAVHeader h1;      // WAV文件头信息（采样率、声道数、位深度等）
    uint32_t len1;     // 音频数据长度（样本数）
    int16_t* data1 = read_wav(file1, &h1, &len1); // 读取音频数据
    if (!data1) return 1; // 读取失败，直接退出

    // 读取第二个音频文件
    WAVHeader h2;      // 第二个文件的WAV头信息
    uint32_t len2;     // 第二个文件的音频数据长度
    int16_t* data2 = read_wav(file2, &h2, &len2); // 读取音频数据
    if (!data2) {
        free(data1);   // 释放第一个文件的内存
        return 1;      // 读取失败，退出程序
    }

    // 验证两个音频文件的格式兼容性
    // 波束形成要求两个文件的采样率、位深度、声道数完全一致
    if (h1.sampleRate != h2.sampleRate ||
        h1.bitsPerSample != h2.bitsPerSample ||
        h1.numChannels != h2.numChannels) {
        fprintf(stderr, "WAV files format incompatible (sample rate/bit depth/channels must match)\n");
        free(data1);   // 释放内存
        free(data2);
        return 1;      // 格式不兼容，退出程序
    }

    // 自动计算延迟（如果不是手动指定延迟）
    if (auto_delay) {
        // 显示延迟估计算法选择菜单
        printf("Choosing delay estimation method:\n");
        printf("1. Time-domain cross-correlation (current)\n");
        printf("2. FFT-PHAT (more accurate for noisy signals)\n");
        printf("Enter choice (1 or 2): ");
        
        // 获取用户选择的延迟估计算法
        int choice;
        if (scanf("%d", &choice) != 1) {
            choice = 1; // 输入无效时，默认使用时域方法
        }
        
        int estimated_delay; // 存储估计的延迟值
        if (choice == 2) {
            // 使用FFT-PHAT频域延迟估计算法
            printf("Using FFT-PHAT delay estimation...\n");
            
            // 将int16_t格式的音频数据转换为float格式
            // FFT-PHAT算法需要浮点数输入以提高精度
            float* float_data1 = (float*)malloc(len1 * sizeof(float));
            float* float_data2 = (float*)malloc(len2 * sizeof(float));
            
            // 归一化：将int16_t范围[-32768, 32767]转换为float范围[-1.0, 1.0]
            for (uint32_t i = 0; i < len1; i++) {
                float_data1[i] = (float)data1[i] / 32768.0f;
                float_data2[i] = (float)data2[i] / 32768.0f;
            }
            
            // FFT-PHAT算法参数设置
            int fft_size = 512;  // FFT窗口大小，影响频率分辨率
            int window = (len1 < fft_size) ? len1 : fft_size; // 实际使用的窗口大小
            int margin = 200;    // 延迟搜索范围：±200个样本
            int peak_num = 1;    // 只需要找到最强的峰值
            
            // 分配存储延迟估计结果的数组
            int* delays = (int*)malloc(peak_num * sizeof(int));       // 存储延迟值
            float* peak_values = (float*)malloc(peak_num * sizeof(float)); // 存储峰值置信度
            
            // 调用FFT-PHAT算法进行延迟估计
            // 参数：延迟数组、峰值数组、峰值数量、信号2、信号1、搜索范围、窗口大小、FFT大小
            float result = FFT_Real_Gcc_Path(delays, peak_values, &peak_num,
                                           float_data2, float_data1,
                                           margin, window, fft_size);
            
            // 检查FFT-PHAT算法是否成功找到延迟
            if (peak_num > 0) {
                estimated_delay = delays[0]; // 获取第一个（最强）峰值对应的延迟
                printf("FFT-PHAT estimated delay: %d samples (confidence: %.6f)\n", 
                       estimated_delay, peak_values[0]); // 显示延迟值和置信度
            } else {
                // FFT-PHAT算法失败，回退到时域方法
                printf("FFT-PHAT failed, falling back to time-domain method\n");
                estimated_delay = estimate_delay(data1, len1, data2, len2, 5000);
            }
            
            // 释放FFT-PHAT算法使用的临时内存
            free(float_data1);  // 释放浮点音频数据1
            free(float_data2);  // 释放浮点音频数据2
            free(delays);       // 释放延迟数组
            free(peak_values);  // 释放峰值数组
        } else {
            // 使用时域互相关延迟估计算法
            printf("Using time-domain delay estimation...\n");
            int max_delay = 5000; // 最大延迟搜索范围：5000个样本
            estimated_delay = estimate_delay(data1, len1, data2, len2, max_delay);
            printf("Time-domain estimated delay: %d samples\n", estimated_delay);
        }
        
        // 显示估计的相对延迟
        // 正值表示第二个信号滞后于第一个信号
        // 负值表示第二个信号超前于第一个信号
        printf("Estimated relative delay: %d samples (positive = second signal lags)\n", estimated_delay);
        
        // 根据估计的延迟值设置两个通道的延迟
        // 原则：保持一个通道延迟为0，另一个通道应用相对延迟
        if (estimated_delay >= 0) {
            delay1 = 0;               // 第一个通道不延迟
            delay2 = estimated_delay;   // 第二个通道延迟估计的样本数
        } else {
            delay1 = -estimated_delay;  // 第一个通道延迟绝对值
            delay2 = 0;               // 第二个通道不延迟
        }
        printf("Applied delay: delay1=%d, delay2=%d\n", delay1, delay2);
    }

    // 执行延迟求和波束形成算法
    // 将两个信号按照计算的延迟值对齐，然后求和平均
    uint32_t outLen;    // 输出信号的长度
    int16_t* outData = delay_sum(data1, len1, delay1, data2, len2, delay2, &outLen);
    if (!outData) {
        free(data1);    // 释放输入数据内存
        free(data2);
        return 1;       // 波束形成失败
    }

    // 将波束形成后的音频数据写入输出文件
    // 使用第一个文件的WAV头作为模板（保持采样率、位深度等参数）
    if (!write_wav(outPath, &h1, outData, outLen)) {
        free(data1);    // 释放所有分配的内存
        free(data2);
        free(outData);
        free(outPath);
        return 1;       // 写入失败
    }

    // 显示处理结果信息
    printf("Successfully generated enhanced WAV file: %s\n", outPath);
    printf("Output samples: %u\n", outLen);

    // 释放所有分配的内存，防止内存泄漏
    free(data1);    // 第一个输入音频数据
    free(data2);    // 第二个输入音频数据
    free(outData);  // 输出音频数据
    free(outPath);  // 输出文件路径字符串
    return 0;       // 程序成功执行
}

/**
 * @brief API函数：分离立体声文件为左右两个单声道文件
 * @param stereoFile 输入的立体声音频文件路径
 * @param leftFile 输出的左声道文件名（不包含路径）
 * @param rightFile 输出的右声道文件名（不包含路径）
 * @return 成功返回1，失败返回0
 * 
 * 功能说明：
 * 1. 自动创建audio_files输出目录
 * 2. 构建完整的输出文件路径
 * 3. 调用底层分离函数执行实际操作
 * 4. 释放路径内存并返回结果
 */
int split_stereo_api(const char* stereoFile, const char* leftFile, const char* rightFile) {
    // Ensure audio output directory exists
    ensure_audio_dir();
    
    // Build full paths for output files
    char* leftPath = build_audio_path(leftFile);
    char* rightPath = build_audio_path(rightFile);

    // For now, return error since split_stereo is not included in main compilation
    // User should use split_stereo.exe separately
    fprintf(stderr, "split_stereo_api not available in main program. Use split_stereo.exe instead.\n");
    
    free(leftPath);
    free(rightPath);
    
    return 0; // Return failure
}

/**
 * @brief API函数：生成带TDOA（到达时间差）的立体声音频
 * @param voiceFile 人声音频文件路径
 * @param sineFile 正弦波（噪声）音频文件路径
 * @param outputFile 输出立体声文件名（不包含路径）
 * @param mic_distance 麦克风间距（米）
 * @param sound_speed 声速（米/秒）
 * @param angle_voice 声源角度（度，0°=正前方，负值=左侧，正值=右侧）
 * @param angle_noise 噪声源角度（度）
 * @return 成功返回1，失败返回0
 * 
 * 功能说明：
 * 1. 模拟双麦克风阵列接收不同角度声源的场景
 * 2. 根据TDOA理论计算不同角度声源到达两个麦克风的时间差
 * 3. 生成包含人声和噪声的立体声测试信号
 * 4. 用于测试波束形成算法的性能
 */
int generate_stereo_with_tdoa_api(const char* voiceFile, const char* sineFile,
                                  const char* outputFile, float mic_distance,
                                  float sound_speed, float angle_voice, float angle_noise) {
    // 确保音频输出目录存在
    ensure_audio_dir();
    
    // 构建输出文件的完整路径
    char* outPath = build_audio_path(outputFile);

    // 显示参数信息，便于调试和验证
    printf("Generating stereo audio with TDOA:\n");
    printf("Voice file: %s\n", voiceFile);
    printf("Sine file: %s\n", sineFile);
    printf("Output file: %s\n", outPath);
    printf("Mic distance: %.3f m\n", mic_distance);
    printf("Sound speed: %.1f m/s\n", sound_speed);
    printf("Voice angle: %.1f°\n", angle_voice);
    printf("Noise angle: %.1f°\n", angle_noise);

    // For now, return error since generate_stereo_with_tdoa is not included in main compilation
    // User should use a separate program for TDOA generation
    fprintf(stderr, "generate_stereo_with_tdoa_api not available in main program.\n");
    
    int result = 0; // Return failure
    
    // 释放输出路径内存
    free(outPath);
    return result; // 返回操作结果
}

/**
 * @brief 程序主入口函数
 * @param argc 命令行参数个数
 * @param argv 命令行参数数组
 * @return 程序执行状态
 * 
 * 主函数直接调用do_work函数处理具体的波束形成逻辑
 * 这种设计便于单元测试和代码重用
 */
int main(int argc, char* argv[]) {
    return do_work(argc, argv);
}

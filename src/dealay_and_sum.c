#include "../include/setting.h"
#include "../include/fft_path.h"
#include "../include/gcc_phat_delay.h"

// 计算互相关，估计第二个信号相对于第一个信号的延迟（样本数）
// 返回值：延迟 d，满足 y[n] ≈ x[n-d] (d>0表示y滞后于x)
int estimate_delay(const int16_t* x, uint32_t len_x,
                   const int16_t* y, uint32_t len_y,
                   int max_delay) {
    // 限制最大搜索范围不超过信号长度
    int max_lag = max_delay;
    if (max_lag > (int)len_x - 1) max_lag = len_x - 1;
    if (max_lag > (int)len_y - 1) max_lag = len_y - 1;
    if (max_lag < 0) max_lag = 0;

    int best_delay = 0;
    double max_corr = -1e9;

    // 时域互相关：对于每个可能的延迟 k，计算 x[0..N-1] 与 y[k..k+N-1] 的点积
    // 为了简单，只使用两个信号重叠部分的最大长度
    int N = (len_x < len_y) ? len_x : len_y;
    // 实际使用重叠长度取决于延迟，但为了效率，我们取固定长度，忽略边界效应
    // 更精确的做法是对于每个 k，取有效重叠部分
    // 这里简化：使用所有重叠样本（动态计算）
    for (int k = -max_lag; k <= max_lag; k++) {
        double corr = 0.0;
        int start_x = (k < 0) ? -k : 0;
        int start_y = (k > 0) ? k : 0;
        int overlap_len = len_x - start_x;
        if (len_y - start_y < overlap_len) overlap_len = len_y - start_y;
        if (overlap_len <= 0) continue;

        for (int i = 0; i < overlap_len; i++) {
            corr += (double)x[start_x + i] * y[start_y + i];
        }
        if (corr > max_corr) {
            max_corr = corr;
            best_delay = k;
        }
    }
    return best_delay;
}




// 延迟求和,两个信号，各自延迟delay1和delay2样本，求和并平均
int16_t* delay_sum(const int16_t* data1, uint32_t len1, int delay1,
                   const int16_t* data2, uint32_t len2, int delay2,
                   uint32_t* outLen) 
{
    // 找到最小延迟，将所有延迟调整为非负
    int minDelay = (delay1 < delay2) ? delay1 : delay2;
    int adjDelay1 = delay1 - minDelay;
    int adjDelay2 = delay2 - minDelay;

    // 输出长度 = max(len1 + adjDelay1, len2 + adjDelay2)
    uint32_t outLen1 = len1 + adjDelay1;
    uint32_t outLen2 = len2 + adjDelay2;
    *outLen = (outLen1 > outLen2) ? outLen1 : outLen2;

    int16_t* outData = (int16_t*)calloc(*outLen, sizeof(int16_t));
    if (!outData) {
        fprintf(stderr, "Memory allocation failed\n");
        return NULL;
    }

    // 求和
    for (uint32_t i = 0; i < *outLen; i++) 
    {
        int32_t sum = 0;
        // 通道1
        if (i >= adjDelay1 && i - adjDelay1 < len1) {
            sum += data1[i - adjDelay1];
        }
        // 通道2
        if (i >= adjDelay2 && i - adjDelay2 < len2) {
            sum += data2[i - adjDelay2];
        }
        // 平均并转换为int16_t
        int16_t val = (int16_t)(sum / 2);
        outData[i] = val;
    }

    return outData;
}

/**
 * @brief 交互式延迟估计
 *
 * 显示算法选择菜单，调用用户指定的方法估算两路信号间延迟。
 * 将 int16_t 转换为 float 的逻辑封装在函数内部，main 无需关心。
 */
int estimate_delay_interactive(const int16_t* data1, uint32_t len1,
                                const int16_t* data2, uint32_t len2,
                                int sample_rate)
{
    printf("Choosing delay estimation method:\n");
    printf("  1. Time-domain cross-correlation\n");
    printf("  2. FFT-PHAT (legacy implementation)\n");
    printf("  3. GCC-PHAT (phase transform, robust)\n");
    printf("Enter choice (1, 2 or 3): ");

    int choice;
    if (scanf("%d", &choice) != 1) {
        choice = 1;
    }
    /* 清空输入缓冲区 */
    int c;
    while ((c = getchar()) != '\n' && c != EOF);

    int estimated_delay = 0;
    uint32_t min_len = (len1 < len2) ? len1 : len2;

    if (choice == 3) {
        /* ---------- GCC-PHAT ---------- */
        printf("Using GCC-PHAT delay estimation...\n");

        float* f1 = (float*)malloc(len1 * sizeof(float));
        float* f2 = (float*)malloc(len2 * sizeof(float));
        if (!f1 || !f2) {
            fprintf(stderr, "Memory allocation failed\n");
            free(f1); free(f2);
            return 0;
        }
        for (uint32_t i = 0; i < len1; i++) f1[i] = (float)data1[i] / 32768.0f;
        for (uint32_t i = 0; i < len2; i++) f2[i] = (float)data2[i] / 32768.0f;

        int fft_size = 1024;
        int max_delay = 50;
        int sig_len = (int)min_len;

        if (sig_len < fft_size) {
            fft_size = 1;
            while (fft_size * 2 <= sig_len) fft_size *= 2;
            printf("  Adjusting fft_size to %d (signal too short)\n", fft_size);
        }
        if (max_delay >= sig_len)
            max_delay = sig_len / 4;

        printf("  Parameters: fft_size=%d, max_delay=%d\n", fft_size, max_delay);
        estimated_delay = compute_gcc_phat_delay(f1, f2, sig_len, max_delay, fft_size);
        printf("GCC-PHAT estimated delay: %d samples\n", estimated_delay);

        free(f1); free(f2);

    } else if (choice == 2) {
        /* ---------- FFT-PHAT (legacy) ---------- */
        printf("Using FFT-PHAT delay estimation...\n");

        float* f1 = (float*)malloc(len1 * sizeof(float));
        float* f2 = (float*)malloc(len2 * sizeof(float));
        if (!f1 || !f2) {
            fprintf(stderr, "Memory allocation failed\n");
            free(f1); free(f2);
            return 0;
        }
        for (uint32_t i = 0; i < len1; i++) f1[i] = (float)data1[i] / 32768.0f;
        for (uint32_t i = 0; i < len2; i++) f2[i] = (float)data2[i] / 32768.0f;

        int fft_size = 512;
        int window = (min_len < (uint32_t)fft_size) ? min_len : (uint32_t)fft_size;
        int margin = 200;
        int peak_num = 1;
        int* delays = (int*)malloc(sizeof(int));
        float* peak_values = (float*)malloc(sizeof(float));
        if (!delays || !peak_values) {
            fprintf(stderr, "Memory allocation failed\n");
            free(f1); free(f2); free(delays); free(peak_values);
            return 0;
        }

        FFT_Real_Gcc_Path(delays, peak_values, &peak_num,
                          f2, f1, margin, window, fft_size);

        if (peak_num > 0) {
            estimated_delay = delays[0];
            printf("FFT-PHAT estimated delay: %d samples (confidence: %.6f)\n",
                   estimated_delay, peak_values[0]);
        } else {
            printf("FFT-PHAT failed, falling back to time-domain method\n");
            estimated_delay = estimate_delay(data1, len1, data2, len2, 5000);
        }

        free(f1); free(f2); free(delays); free(peak_values);

    } else {
        /* ---------- Time-domain ---------- */
        printf("Using time-domain delay estimation...\n");
        int max_delay = 5000;
        if ((int)len1 - 1 < max_delay) max_delay = (int)len1 - 1;
        if ((int)len2 - 1 < max_delay) max_delay = (int)len2 - 1;
        if (max_delay < 0) max_delay = 0;
        estimated_delay = estimate_delay(data1, len1, data2, len2, max_delay);
        printf("Time-domain estimated delay: %d samples\n", estimated_delay);
    }

    printf("Estimated relative delay: %d samples (positive = second signal lags)\n", estimated_delay);
    return estimated_delay;
}
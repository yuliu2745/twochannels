#include "../include/setting.h"

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
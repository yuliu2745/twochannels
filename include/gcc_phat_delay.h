#ifndef GCC_PHAT_DELAY_H
#define GCC_PHAT_DELAY_H

#include <fftw3.h>

/**
 * @brief GCC-PHAT 上下文
 *
 * 所有缓冲区、窗系数、FFTW plan 在 init 时一次性创建并复用，
 * 避免每次 compute 重复分配/释放/重建 plan。
 */
typedef struct {
    int fft_size;           /* FFT 窗口大小 */
    int input_len;          /* 零填充后传给 FFTW 的实数输入长度 = 2 * fft_size */

    float *in1, *in2;              /* 实数输入缓冲区 (input_len) */
    int    complex_len;            /* r2c 输出的复数 bin 数 = input_len/2 + 1 = fft_size + 1 */
    fftwf_complex *out1, *out2;    /* 前向 FFT 输出 (complex_len) */
    fftwf_complex *cs;             /* 归一化互功率谱   (complex_len) */
    float *xcorr;                  /* IFFT 结果        (input_len) */
    float *window_coeffs;          /* 预计算 Hamming 窗 (fft_size) */

    fftwf_plan plan_fwd1, plan_fwd2;   /* 两路 r2c 前向 plan */
    fftwf_plan plan_inv;               /* c2r 逆 plan */
} GccPhatContext;

/**
 * @brief 初始化 GCC-PHAT 上下文
 *
 * 分配所有缓冲区，预计算 Hamming 窗系数，
 * 以 FFTW_MEASURE 模式创建 plan（首次调用较慢但会深度优化）。
 *
 * @param ctx      上下文指针（未初始化时调用）
 * @param fft_size  FFT 窗口大小（内部零填充至 2*fft_size 做互相关）
 * @return 0 成功，-1 失败
 */
int gcc_phat_init(GccPhatContext *ctx, int fft_size);

/**
 * @brief 执行 GCC-PHAT 时延估计
 *
 * 使用 ctx 中预分配的资源和 plan 完成计算，不分配内存、不创建 plan。
 *
 * @param ctx       已初始化的上下文
 * @param sig1      参考信号（float 数组，长度 len）
 * @param sig2      延迟信号（float 数组，长度 len）
 * @param len       两路信号长度
 * @param max_delay 最大期望时延（样点，搜索范围为 ±max_delay）
 * @return 时延估计值（样点，正= sig2 滞后于 sig1，负=超前，0.0f=失败）
 *         含亚样点精度（三点抛物线插值）
 */
float gcc_phat_compute(GccPhatContext *ctx,
                       const float *sig1, const float *sig2,
                       int len, int max_delay,
                       int sample_rate);

/**
 * @brief 销毁 GCC-PHAT 上下文，释放所有资源
 */
void gcc_phat_destroy(GccPhatContext *ctx);

#endif /* GCC_PHAT_DELAY_H */

/**
 * @file gcc_phat_delay.c
 * @brief GCC-PHAT time delay estimation — context-based API.
 *
 * Algorithm steps (inside gcc_phat_compute):
 *   1. Apply Hamming window to both signals (precomputed coefficients)
 *   2. Zero-pad to 2 * fft_size
 *   3. Real FFT both signals -> frequency domain (pre-created plans)
 *   4. Cross-power spectrum: G(f) = X(f) * conj(Y(f))
 *   5. PHAT normalisation: G(f) / |G(f)|  (phase only)
 *   6. Inverse FFT -> time-domain cross-correlation
 *   7. Peak search within [-max_delay, +max_delay]
 *   8. Return delay at maximum peak + sub-sample refinement
 */

#include "../include/gcc_phat_delay.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ------------------------------------------------------------------ */
/*  Init                                                               */
/* ------------------------------------------------------------------ */

int gcc_phat_init(GccPhatContext *ctx, int fft_size)
{
    if (!ctx || fft_size < 8)
        return -1;

    memset(ctx, 0, sizeof(*ctx));

    ctx->fft_size      = fft_size;
    ctx->input_len     = 2 * fft_size;
    ctx->complex_len   = fft_size + 1;   /* r2c 输出复数 bin 数 = input_len/2 + 1 */

    int nbins = ctx->complex_len;        /* 分配的复数缓冲区大小 */

    /* ---------- 分配缓冲区 ---------- */
    ctx->in1          = (float *)fftwf_malloc(sizeof(float)       * ctx->input_len);
    ctx->in2          = (float *)fftwf_malloc(sizeof(float)       * ctx->input_len);
    ctx->out1         = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * nbins);
    ctx->out2         = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * nbins);
    ctx->cs           = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * nbins);
    ctx->xcorr        = (float *)fftwf_malloc(sizeof(float)       * ctx->input_len);
    ctx->window_coeffs = (float *)fftwf_malloc(sizeof(float)       * fft_size);

    if (!ctx->in1 || !ctx->in2 || !ctx->out1 || !ctx->out2 ||
        !ctx->cs || !ctx->xcorr || !ctx->window_coeffs)
        goto fail;

    /* ---------- 预计算 Hamming 窗 ---------- */
    for (int i = 0; i < fft_size; i++) {
        ctx->window_coeffs[i] = 0.54f - 0.46f * cosf(2.0f * (float)M_PI * i / (fft_size - 1));
    }

    /* ---------- 创建 FFTW plan (FFTW_MEASURE 做深度优化) ---------- */
    ctx->plan_fwd1 = fftwf_plan_dft_r2c_1d(ctx->input_len,
                                           ctx->in1, ctx->out1, FFTW_MEASURE);
    ctx->plan_fwd2 = fftwf_plan_dft_r2c_1d(ctx->input_len,
                                           ctx->in2, ctx->out2, FFTW_MEASURE);
    ctx->plan_inv  = fftwf_plan_dft_c2r_1d(ctx->input_len,
                                           ctx->cs, ctx->xcorr, FFTW_MEASURE);

    if (!ctx->plan_fwd1 || !ctx->plan_fwd2 || !ctx->plan_inv)
        goto fail;

    return 0;

fail:
    gcc_phat_destroy(ctx);
    return -1;
}

/* ------------------------------------------------------------------ */
/*  Compute                                                            */
/* ------------------------------------------------------------------ */

float gcc_phat_compute(GccPhatContext *ctx,
                       const float *sig1, const float *sig2,
                       int len, int max_delay,
                       int sample_rate)
{
    if (!ctx || !sig1 || !sig2 || len < 2 || max_delay < 1)
        return 0.0f;

    int fft_size      = ctx->fft_size;
    int input_len = ctx->input_len;
    int window_len = (len < fft_size) ? len : fft_size;

    /*DC去除，避免0延迟假峰*/
    double sum1 = 0.0, sum2 = 0.0;
    for (int i = 0; i < window_len; i++) {
            sum1 += (double)sig1[i];
        sum2 += (double)sig2[i];
    }
    double mean1 = sum1 / window_len;
    double mean2 = sum2 / window_len;

    /* ---------- 加窗 + 零填充 ---------- */
    for (int i = 0; i < window_len; i++) {
        float w = ctx->window_coeffs[i];
        ctx->in1[i] = (sig1[i] - (float)mean1) * w;
        ctx->in2[i] = (sig2[i] - (float)mean2) * w;
    }
    for (int i = window_len; i < input_len; i++) {
        ctx->in1[i] = 0.0f;
        ctx->in2[i] = 0.0f;
    }

    /* ---------- FFT ---------- */
    fftwf_execute(ctx->plan_fwd1);
    fftwf_execute(ctx->plan_fwd2);

    /* ---------- GCC-PHAT 核心 ---------- */
    int nbins = ctx->complex_len;
    for (int i = 0; i < nbins; i++) {
        float re_x = ctx->out1[i][0];
        float im_x = ctx->out1[i][1];
        float re_y = ctx->out2[i][0];
        float im_y = ctx->out2[i][1];

        /* 互功率谱 X · conj(Y) */
        float re_cs = re_x * re_y + im_x * im_y;
        float im_cs = im_x * re_y - re_x * im_y;

        /* PHAT 归一化 */
        float mag = sqrtf(re_cs * re_cs + im_cs * im_cs);
        if (mag > 1e-12f) {
            ctx->cs[i][0] = re_cs / mag;
            ctx->cs[i][1] = im_cs / mag;
        } else {
            ctx->cs[i][0] = 0.0f;
            ctx->cs[i][1] = 0.0f;
        }
    }

    // ====================== 层1：300~3400Hz 带通掩码（优化时延估计精度） ======================
    // 使用调用方传入的 sample_rate，而非硬编码值
    float bin_hz = (float)sample_rate / (float)ctx->input_len; // input_len = 2*fft_size，R2C频点分辨率
    int bin_low  = (int)ceilf(300.0f / bin_hz);
    int bin_high = (int)floorf(3400.0f / bin_hz);
    for(int k = 0; k < nbins; k++)
    {
        if(k < bin_low || k > bin_high)
        {
            ctx->cs[k][0] = 0.0f;
            ctx->cs[k][1] = 0.0f;
        }
    }
    
    /* ---------- 逆 FFT ---------- */
    fftwf_execute(ctx->plan_inv);

    /* ---------- 峰值搜索 ---------- */
    int best_lag = 0;
    float best_val = ctx->xcorr[0];   /* FFTW 未缩放，但比较时除 N 不影响峰值位置 */

    for (int i = 1; i <= max_delay && i < input_len; i++) {
        if (ctx->xcorr[i] > best_val) {
            best_val = ctx->xcorr[i];
            best_lag = i;
        }
    }
    for (int i = input_len - max_delay; i < input_len; i++) {
        if (i < 0) continue;
        if (ctx->xcorr[i] > best_val) {
            best_val = ctx->xcorr[i];
            best_lag = i - input_len;   /* 负延迟区间，回绕 */
        }
    }

    /* ---------- 亚样点抛物线插值 ---------- */
    float sub_delta = 0.0f;
    if (best_lag > -max_delay && best_lag < max_delay) {
        int idx = (best_lag >= 0) ? best_lag : (best_lag + input_len);

        /*
         * 循环互相关：idx=0 的左邻域是 idx=input_len-1，
         * idx=input_len-1 的右邻域是 idx=0。
         */
        int idx_m = (idx == 0) ? (input_len - 1) : (idx - 1);
        int idx_p = (idx == input_len - 1) ? 0 : (idx + 1);

        float y0 = ctx->xcorr[idx];
        float y1 = ctx->xcorr[idx_p];
        float ym = ctx->xcorr[idx_m];

        float denom = 2.0f * y0 - y1 - ym;
        if (fabsf(denom) > 1e-12f) {
            float delta = 0.5f * (y1 - ym) / denom;
            if (delta > -0.5f && delta < 0.5f) {
                sub_delta = delta;
            }
        }
    }

    return (float)best_lag + sub_delta;
}

/* ------------------------------------------------------------------ */
/*  Destroy                                                            */
/* ------------------------------------------------------------------ */

void gcc_phat_destroy(GccPhatContext *ctx)
{
    if (!ctx) return;

    if (ctx->plan_fwd1) fftwf_destroy_plan(ctx->plan_fwd1);
    if (ctx->plan_fwd2) fftwf_destroy_plan(ctx->plan_fwd2);
    if (ctx->plan_inv)  fftwf_destroy_plan(ctx->plan_inv);

    if (ctx->in1)          fftwf_free(ctx->in1);
    if (ctx->in2)          fftwf_free(ctx->in2);
    if (ctx->out1)         fftwf_free(ctx->out1);
    if (ctx->out2)         fftwf_free(ctx->out2);
    if (ctx->cs)           fftwf_free(ctx->cs);
    if (ctx->xcorr)        fftwf_free(ctx->xcorr);
    if (ctx->window_coeffs) fftwf_free(ctx->window_coeffs);

    memset(ctx, 0, sizeof(*ctx));
}

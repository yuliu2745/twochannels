#include "../include/signal_restore.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define SR_EPSILON 1e-10f

/* ------------------------------------------------------------------ */
/*  Init                                                               */
/* ------------------------------------------------------------------ */

SignalRestoreContext* restore_init(int complex_len, float lambda, float eta, float eta_low,
                                   float beta_low, float beta_high, float snr_thresh_db)
{
    SignalRestoreContext *ctx = (SignalRestoreContext *)calloc(1, sizeof(SignalRestoreContext));
    if (!ctx) return NULL;

    ctx->complex_len    = complex_len;
    ctx->lambda         = lambda;
    ctx->eta            = eta;
    ctx->eta_low        = eta_low;
    ctx->beta_low       = beta_low;
    ctx->beta_high      = beta_high;
    ctx->snr_thresh_db  = snr_thresh_db;
    /* initialized 由 calloc 置 0 */

    return ctx;
}

/* ------------------------------------------------------------------ */
/*  帧级余弦相似度:  Σ Re(X* · Y) / (||X|| · ||Y||)                   */
/* ------------------------------------------------------------------ */

static float cos_sim(const fftwf_complex *X, const fftwf_complex *Y, int nbins)
{
    double dot = 0.0, nrm_x = 0.0, nrm_y = 0.0;

    for (int k = 0; k < nbins; k++) {
        float re_x = X[k][0], im_x = X[k][1];
        float re_y = Y[k][0], im_y = Y[k][1];
        dot   += (double)(re_x * re_y + im_x * im_y);
        nrm_x += (double)(re_x * re_x + im_x * im_x);
        nrm_y += (double)(re_y * re_y + im_y * im_y);
    }

    float denom = (float)(sqrt(nrm_x) * sqrt(nrm_y));
    if (denom < SR_EPSILON) return 0.0f;
    return (float)(dot / denom);
}

/* ------------------------------------------------------------------ */
/*  Process one frame                                                  */
/* ------------------------------------------------------------------ */

void restore_process_frame(SignalRestoreContext *ctx,
                           const fftwf_complex *Y_FBF,
                           const fftwf_complex *Y_U,
                           const fftwf_complex *X1,
                           const fftwf_complex *X2,
                           fftwf_complex *Y_out,
                           int nbins)
{
    if (!ctx || !Y_FBF || !Y_U || !X1 || !X2 || !Y_out)
        return;

    int n = (nbins < ctx->complex_len) ? nbins : ctx->complex_len;
    float lam = ctx->lambda;

    /* ---- Step 1: 帧级余弦相似度 ---- */
    float c_x1f = cos_sim(X1, Y_FBF, n);
    float c_x2f = cos_sim(X2, Y_FBF, n);
    float c_x1u = cos_sim(X1, Y_U,   n);
    float c_x2u = cos_sim(X2, Y_U,   n);

    /* 时间平滑 (lambda=0.88~0.90) */
    if (!ctx->initialized) {
        ctx->S_X1_Y_FBF = c_x1f;
        ctx->S_X2_Y_FBF = c_x2f;
        ctx->S_X1_Y_U   = c_x1u;
        ctx->S_X2_Y_U   = c_x2u;
        ctx->initialized = 1;
    } else {
        ctx->S_X1_Y_FBF = lam * ctx->S_X1_Y_FBF + (1.0f - lam) * c_x1f;
        ctx->S_X2_Y_FBF = lam * ctx->S_X2_Y_FBF + (1.0f - lam) * c_x2f;
        ctx->S_X1_Y_U   = lam * ctx->S_X1_Y_U   + (1.0f - lam) * c_x1u;
        ctx->S_X2_Y_U   = lam * ctx->S_X2_Y_U   + (1.0f - lam) * c_x2u;
    }

    /* ---- Step 2: 帧级 SNR 计算 → 自适应 beta ---- */
    double p_fbf_sum = 0.0, p_yu_sum = 0.0;
    for (int k = 0; k < n; k++) {
        p_fbf_sum += (double)(Y_FBF[k][0]*Y_FBF[k][0] + Y_FBF[k][1]*Y_FBF[k][1]);
        p_yu_sum  += (double)(Y_U[k][0]*Y_U[k][0]   + Y_U[k][1]*Y_U[k][1]);
    }
    float snr_db = 10.0f * log10f((float)(p_fbf_sum + SR_EPSILON) / (float)(p_yu_sum + SR_EPSILON));
    float beta = (snr_db > ctx->snr_thresh_db) ? ctx->beta_high : ctx->beta_low;

    /* ---- Step 3: 三段式决策 ---- */
    float d12_fbf = ctx->S_X1_Y_FBF - ctx->S_X2_Y_FBF;
    float SD = fabsf(d12_fbf);

    int sel = 0;
    if (SD > ctx->eta) {
        /* 相似度差值大 → 选匹配度高的单麦 */
        sel = (d12_fbf > 0.0f) ? 1 : 2;
    } else if (SD > ctx->eta_low) {
        /* 中性区间 → 强制 (X1+X2)/2，绝不闭锁 */
        sel = 3;
    } else {
        /* 强干扰 → 关闭补偿，仅保留 FBF */
        sel = 0;
    }

    /* ---- Step 4: 逐 bin 增益 + 恢复叠加 ---- */
    for (int k = 0; k < n; k++) {
        /* 选择 X_selected */
        float re_xs = 0.0f, im_xs = 0.0f;
        switch (sel) {
        case 1:
            re_xs = X1[k][0]; im_xs = X1[k][1];
            break;
        case 2:
            re_xs = X2[k][0]; im_xs = X2[k][1];
            break;
        case 3:
            re_xs = (X1[k][0] + X2[k][0]) * 0.5f;
            im_xs = (X1[k][1] + X2[k][1]) * 0.5f;
            break;
        default: /* 不恢复 */
            break;
        }

        /* G_X = sqrt(|Y_FBF|² / (β·|Y_U|² + |Y_FBF|² + ε)) */
        float p_fbf = Y_FBF[k][0]*Y_FBF[k][0] + Y_FBF[k][1]*Y_FBF[k][1];
        float p_yu  = Y_U[k][0]*Y_U[k][0]   + Y_U[k][1]*Y_U[k][1];

        float gx = sqrtf(p_fbf / (beta * p_yu + p_fbf + SR_EPSILON));

        /* Ȳ_FBF = Y_FBF + G_X · X_selected */
        Y_out[k][0] = Y_FBF[k][0] + gx * re_xs;
        Y_out[k][1] = Y_FBF[k][1] + gx * im_xs;
    }
}

/* ------------------------------------------------------------------ */
/*  Destroy                                                            */
/* ------------------------------------------------------------------ */

void restore_destroy(SignalRestoreContext *ctx)
{
    if (!ctx) return;
    memset(ctx, 0, sizeof(SignalRestoreContext));
    free(ctx);
}

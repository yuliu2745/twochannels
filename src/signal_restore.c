#include "../include/signal_restore.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define SR_EPSILON 1e-10f

/* ------------------------------------------------------------------ */
/*  Init                                                               */
/* ------------------------------------------------------------------ */

SignalRestoreContext* restore_init(int complex_len, float lambda, float eta, float beta)
{
    SignalRestoreContext *ctx = (SignalRestoreContext *)calloc(1, sizeof(SignalRestoreContext));
    if (!ctx) return NULL;

    ctx->complex_len = complex_len;
    ctx->lambda = lambda;
    ctx->eta    = eta;
    ctx->beta   = beta;
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

/* ---- sign function ---- */
static inline float sgn(float x)
{
    if (x > 0.0f) return 1.0f;
    if (x < 0.0f) return -1.0f;
    return 0.0f;
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

    /* 时间平滑 */
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

    /* ---- Step 2: 决策差值 ---- */
    float d12_fbf = ctx->S_X1_Y_FBF - ctx->S_X2_Y_FBF;   /* SD_X12_Y_FBF */
    float d21_u   = ctx->S_X2_Y_U   - ctx->S_X1_Y_U;     /* SD_X21_Y_U   */

    float eta = ctx->eta;

    /* 选择: 0=不恢复, 1=X1, 2=X2, 3=(X1+X2)/2 */
    int sel = 0;

    if (d12_fbf > eta && d21_u >= 0.0f) {
        sel = 1;                            /* X1 */
    } else if (d12_fbf <= -eta && d21_u <= 0.0f) {
        sel = 2;                            /* X2 */
    } else {
        float prod = d12_fbf * sgn(d21_u);
        if (prod > 0.0f && prod <= eta) {
            sel = 3;                        /* (X1+X2)/2 */
        }
    }

    /* ---- Step 3: 逐 bin 增益 + 恢复叠加 ---- */
    for (int k = 0; k < n; k++) {
        /* 选择 X_selected */
        float re_xs, im_xs;
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
            re_xs = 0.0f; im_xs = 0.0f;
            break;
        }

        /* G_X = sqrt(|Y_FBF|² / (β·|Y_U|² + |Y_FBF|²)) */
        float p_fbf = Y_FBF[k][0]*Y_FBF[k][0] + Y_FBF[k][1]*Y_FBF[k][1];
        float p_yu  = Y_U[k][0]*Y_U[k][0]   + Y_U[k][1]*Y_U[k][1];

        float gx = sqrtf(p_fbf / (ctx->beta * p_yu + p_fbf + SR_EPSILON));

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

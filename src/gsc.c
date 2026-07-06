/**
 * @file gsc.c
 * @brief GSC (Generalized Sidelobe Canceller) — frequency-domain implementation
 *
 * Algorithm per frequency bin (each frame):
 *   Speech ref:  S = (X1 + X2_rotated) / 2
 *   Noise ref:   N = (X1 - X2_rotated) / 2
 *   GSC output:  Y = S - W * N
 *
 * Adaptive filter W updated via leaky NLMS (bin-wise VAD gated):
 *   Pn = alpha * Pn + (1-alpha) * |N|^2
 *   W += mu * conj(N) * Y / (Pn + epsilon)
 *
 * Per-bin VAD: when |S|^2 > vad_thresh * |N|^2, freeze W update
 * to avoid cancelling the desired speech.
 */

#include "../include/gsc.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define GSC_EPSILON 1e-10f

/* ------------------------------------------------------------------ */
/*  Init                                                               */
/* ------------------------------------------------------------------ */

GscContext* gsc_init(int complex_len, float mu, float alpha, float vad_thresh)
{
    if (complex_len < 2)
        return NULL;

    GscContext *gsc = (GscContext *)calloc(1, sizeof(GscContext));
    if (!gsc)
        return NULL;

    gsc->W = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * complex_len);
    gsc->Pn = (float *)fftwf_malloc(sizeof(float) * complex_len);

    if (!gsc->W || !gsc->Pn) {
        gsc_destroy(gsc);
        return NULL;
    }

    /* 权重初始化为 0 (calloc 已经清空) */
    memset(gsc->W,  0, sizeof(fftwf_complex) * complex_len);
    memset(gsc->Pn, 0, sizeof(float)          * complex_len);

    gsc->complex_len  = complex_len;
    gsc->mu           = mu;
    gsc->alpha        = alpha;
    gsc->vad_thresh   = vad_thresh;

    return gsc;
}

/* ------------------------------------------------------------------ */
/*  Process one frame                                                  */
/* ------------------------------------------------------------------ */

void gsc_process_frame(GscContext *gsc,
                       const fftwf_complex *X1,
                       const fftwf_complex *X2,
                       fftwf_complex *Y,
                       int nbins)
{
    if (!gsc || !X1 || !X2 || !Y)
        return;

    int n = (nbins < gsc->complex_len) ? nbins : gsc->complex_len;

    float mu    = gsc->mu;
    float alpha = gsc->alpha;
    float beta  = gsc->vad_thresh;

    for (int k = 0; k < n; k++)
    {
        /* ---- 读取输入频谱 ---- */
        float re_x1 = X1[k][0];
        float im_x1 = X1[k][1];
        float re_x2 = X2[k][0];
        float im_x2 = X2[k][1];

        /* ---- Speech reference: S = (X1 + X2) / 2 ---- */
        float re_s = (re_x1 + re_x2) * 0.5f;
        float im_s = (im_x1 + im_x2) * 0.5f;

        /* ---- Noise reference:  N = (X1 - X2) / 2 ---- */
        float re_n = (re_x1 - re_x2) * 0.5f;
        float im_n = (im_x1 - im_x2) * 0.5f;

        /* ---- 读取自适应权重 W[k] ---- */
        float re_w = gsc->W[k][0];
        float im_w = gsc->W[k][1];

        /* ---- GSC 输出: Y = S - W * N ---- */
        /* W * N = (Wr + jWi)*(Nr + jNi) = (Wr*Nr - Wi*Ni) + j(Wr*Ni + Wi*Nr) */
        float re_wn = re_w * re_n - im_w * im_n;
        float im_wn = re_w * im_n + im_w * re_n;

        float re_y = re_s - re_wn;
        float im_y = im_s - im_wn;

        /* ---- 频点级简易 VAD ---- */
        float pow_s = re_s * re_s + im_s * im_s;
        float pow_n = re_n * re_n + im_n * im_n;

        int speech_present = (pow_s > beta * pow_n + GSC_EPSILON) ? 1 : 0;

        /* ---- 泄漏 NLMS 更新（仅无语音时更新 W） ---- */
        if (!speech_present)
        {
            /* 平滑噪声参考功率: Pn = α * Pn + (1-α) * |N|^2 */
            gsc->Pn[k] = alpha * gsc->Pn[k] + (1.0f - alpha) * pow_n;

            /* NLMS: ΔW = μ * conj(N) * Y / (Pn + ε) */
            /* conj(N) * Y = (Nr - jNi) * (Yr + jYi)
             *            = (Nr*Yr + Ni*Yi) + j(Nr*Yi - Ni*Yr) */
            float denom = gsc->Pn[k] + GSC_EPSILON;
            float scale = mu / denom;

            float re_dw = (re_n * re_y + im_n * im_y) * scale;
            float im_dw = (re_n * im_y - im_n * re_y) * scale;

            gsc->W[k][0] += re_dw;
            gsc->W[k][1] += im_dw;
        }

        /* ---- 写入输出 Y（乘以 2 保持与原始 DSB 幅度一致性） ---- */
        Y[k][0] = re_y;
        Y[k][1] = im_y;
    }
}

/* ------------------------------------------------------------------ */
/*  Destroy                                                            */
/* ------------------------------------------------------------------ */

void gsc_destroy(GscContext *gsc)
{
    if (!gsc) return;

    if (gsc->W)  fftwf_free(gsc->W);
    if (gsc->Pn) fftwf_free(gsc->Pn);

    memset(gsc, 0, sizeof(GscContext));
    free(gsc);
}

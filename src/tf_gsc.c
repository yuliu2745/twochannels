/**
 * @file tf_gsc.c
 * @brief TF_GSC (Transfer Function GSC) — frequency-domain implementation
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
 *
 * BM leakage suppression via per-bin Wiener gain G_U.
 */

#include "../include/tf_gsc.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define TF_GSC_EPSILON 1e-10f

/* ------------------------------------------------------------------ */
/*  Init                                                               */
/* ------------------------------------------------------------------ */

TfGscContext* tf_gsc_init(int complex_len, float mu, float alpha, float vad_thresh, float leak, float smooth_factor)
{
    if (complex_len < 2)
        return NULL;

    TfGscContext *tf_gsc = (TfGscContext *)calloc(1, sizeof(TfGscContext));
    if (!tf_gsc)
        return NULL;

    tf_gsc->W = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * complex_len);
    tf_gsc->Pn = (float *)fftwf_malloc(sizeof(float) * complex_len);
    tf_gsc->G_smooth = (float *)fftwf_malloc(sizeof(float) * complex_len);
    tf_gsc->S_floor = (float *)fftwf_malloc(sizeof(float) * complex_len);

    if (!tf_gsc->W || !tf_gsc->Pn || !tf_gsc->G_smooth || !tf_gsc->S_floor) {
        tf_gsc_destroy(tf_gsc);
        return NULL;
    }

    /* 权重初始化为 0 (calloc 已清空) */
    memset(tf_gsc->W,  0, sizeof(fftwf_complex) * complex_len);
    memset(tf_gsc->Pn, 0, sizeof(float)          * complex_len);
    memset(tf_gsc->G_smooth, 0, sizeof(float)    * complex_len);
    memset(tf_gsc->S_floor, 0, sizeof(float)     * complex_len);

    tf_gsc->complex_len    = complex_len;
    tf_gsc->mu             = mu;
    tf_gsc->alpha          = alpha;
    tf_gsc->vad_thresh     = vad_thresh;
    tf_gsc->leak           = leak;
    tf_gsc->smooth_factor  = smooth_factor;

    return tf_gsc;
}

/* ------------------------------------------------------------------ */
/*  Process one frame                                                  */
/* ------------------------------------------------------------------ */

void tf_gsc_process_frame(TfGscContext *tf_gsc,
                          const fftwf_complex *X1,
                          const fftwf_complex *X2,
                          fftwf_complex *Y_FBF,
                          fftwf_complex *Y_U,
                          fftwf_complex *Y,
                          int nbins)
{
    if (!tf_gsc || !X1 || !X2 || !Y)
        return;

    int n = (nbins < tf_gsc->complex_len) ? nbins : tf_gsc->complex_len;

    float mu    = tf_gsc->mu;
    float alpha = tf_gsc->alpha;
    float beta  = tf_gsc->vad_thresh;
    float leak  = tf_gsc->leak;
    float sf    = tf_gsc->smooth_factor;

    for (int k = 0; k < n; k++)
    {
        /* ---- 读取输入频谱 ---- */
        float re_x1 = X1[k][0];
        float im_x1 = X1[k][1];
        float re_x2 = X2[k][0];
        float im_x2 = X2[k][1];

        /* ---- FBF: Speech reference S = (X1 + X2) / 2 ---- */
        float re_s = (re_x1 + re_x2) * 0.5f;
        float im_s = (im_x1 + im_x2) * 0.5f;

        /* 若需要，输出 S */
        if (Y_FBF) {
            Y_FBF[k][0] = re_s;
            Y_FBF[k][1] = im_s;
        }

        /* ---- BM: Noise reference  N = (X1 - X2) / 2 ---- */
        float re_n = (re_x1 - re_x2) * 0.5f;
        float im_n = (im_x1 - im_x2) * 0.5f;

        /* ---- BM 泄漏抑制（S噪声底归一化，净化噪声参考） ---- */
        float pow_s = re_s * re_s + im_s * im_s;
        float raw_pow_n = re_n * re_n + im_n * im_n;  // 保存原始 N 能量供 VAD 使用

        // 追踪 S 的噪声底：最小值追踪器
        //   S_floor 作为 |S|² 中"噪声分量"的估计
        //   |S|² - S_floor ≈ 纯语音能量（去除了噪声基底）
        if (tf_gsc->S_floor[k] < 1e-15f) {
            tf_gsc->S_floor[k] = pow_s;
        } else if (pow_s < tf_gsc->S_floor[k]) {
            tf_gsc->S_floor[k] = pow_s;
        } else {
            tf_gsc->S_floor[k] += 1e-3f * tf_gsc->S_floor[k] + 1e-12f;
        }
        float s_floor = tf_gsc->S_floor[k];

        // G_U = |N|² / (max(|S|² - S_floor, 0) + |N|² + ε)
        //   |S|² - S_floor = 去除噪声基底后的纯语音能量
        //     语音帧：|S|² >> s_floor → 分母大 → G_U ≈ 0（抑制BM泄漏）
        //     噪声帧：|S|² ≈ s_floor → 分母→|N|² → G_U ≈ 1（保留噪声供ANC）
        float speech_energy = pow_s - s_floor;
        if (speech_energy < 0.0f) speech_energy = 0.0f;
        float G_U = raw_pow_n / (speech_energy + raw_pow_n + TF_GSC_EPSILON);

        // 时间平滑，防止音乐噪声
        tf_gsc->G_smooth[k] = sf * tf_gsc->G_smooth[k] + (1.0f - sf) * G_U;
        float g = tf_gsc->G_smooth[k];

        // 净化后的噪声参考 Ñ = G_U · N
        re_n *= g;
        im_n *= g;
        float cleaned_pow_n = re_n * re_n + im_n * im_n;

        // 若需要，输出 Ñ（在前两句增益应用之后、自适应滤波之前）
        if (Y_U) {
            Y_U[k][0] = re_n;
            Y_U[k][1] = im_n;
        }

        /* ---- 读取自适应权重 W[k] ---- */
        float re_w = tf_gsc->W[k][0];
        float im_w = tf_gsc->W[k][1];

        /* ---- 首帧用净化后 |Ñ|² 初始化 Pn ---- */
        if (!tf_gsc->initialized)
        {
            tf_gsc->Pn[k] = cleaned_pow_n;
            Y[k][0] = re_s;   /* W=0 → Y = S */
            Y[k][1] = im_s;
            if (k == n - 1)
                tf_gsc->initialized = 1;
            continue;
        }

        /* ---- 频点级 VAD（用原始 N 能量判定，不被 G_U 抑制影响） ---- */
        int speech_present = (pow_s > beta * raw_pow_n + TF_GSC_EPSILON) ? 1 : 0;

        /* ---- GSC 输出: Y = S - W * Ñ  (使用净化噪声参考) ---- */
        float re_wn = re_w * re_n - im_w * im_n;
        float im_wn = re_w * im_n + im_w * re_n;

        float re_y = re_s - re_wn;
        float im_y = im_s - im_wn;

        /* ---- NLMS 自适应更新（仅无语音时，使用净化后 Ñ） ---- */
        if (!speech_present)
        {
            tf_gsc->Pn[k] = alpha * tf_gsc->Pn[k] + (1.0f - alpha) * cleaned_pow_n;

            float denom = tf_gsc->Pn[k] + TF_GSC_EPSILON;
            float scale = mu / denom;

            float re_dw = (re_n * re_y + im_n * im_y) * scale;
            float im_dw = (re_n * im_y - im_n * re_y) * scale;

            tf_gsc->W[k][0] += re_dw;
            tf_gsc->W[k][1] += im_dw;
        }
        else
        {
            /* 语音帧：W 缓慢泄漏回零，防止累积误差持续抵消语音 */
            tf_gsc->W[k][0] *= leak;
            tf_gsc->W[k][1] *= leak;
        }

        Y[k][0] = re_y;
        Y[k][1] = im_y;
    }
}

/* ------------------------------------------------------------------ */
/*  Destroy                                                            */
/* ------------------------------------------------------------------ */

void tf_gsc_destroy(TfGscContext *tf_gsc)
{
    if (!tf_gsc) return;

    if (tf_gsc->W)        fftwf_free(tf_gsc->W);
    if (tf_gsc->Pn)       fftwf_free(tf_gsc->Pn);
    if (tf_gsc->G_smooth) fftwf_free(tf_gsc->G_smooth);
    if (tf_gsc->S_floor)  fftwf_free(tf_gsc->S_floor);

    memset(tf_gsc, 0, sizeof(TfGscContext));
    free(tf_gsc);
}

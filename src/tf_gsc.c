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
 *
 * AMC (Adaptive Mode Controller):
 *   Spatial SNR = 10*log10(FBF_power / BM_power) — 双麦空间信噪比
 *   Coherence   = 帧级宽带相干系数
 *   Silence (空间 SNR<6dB + 相干性>0.85): μ=0, 冻结 W
 *   Strong speech (空间 SNR>12dB):         μ=基准步长
 *   Weak / reverb / multi (else):          μ=0.5×基准步长
 *   Weight clamp: |W[k]| < W_max 防止发散
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
    tf_gsc->onset_thresh   = 3.0f;   /* 能量跳变 3× 触发语音起振保护 */

    /* AMC 默认值 */
    tf_gsc->amc_mu         = mu;     /* 初始 = 基准步长 */
    tf_gsc->W_max          = 10.0f;  /* 权重幅值钳位 ±10（约 20dB） */
    tf_gsc->coh_smooth     = 0.0f;

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

    /* ============================================================ */
    /*  Phase 1: 语音起振检测 + AMC 场景判别                            */
    /* ============================================================ */
    int   in_onset = 0;
    float fbf_energy = 0.0f, bm_energy = 0.0f;
    float sum_cross_mag = 0.0f, sum_pow1 = 0.0f, sum_pow2 = 0.0f;

    for (int k = 0; k < n; k++) {
        float re_x1 = X1[k][0], im_x1 = X1[k][1];
        float re_x2 = X2[k][0], im_x2 = X2[k][1];

        float re_s = (re_x1 + re_x2) * 0.5f;
        float im_s = (im_x1 + im_x2) * 0.5f;
        float re_n = (re_x1 - re_x2) * 0.5f;
        float im_n = (im_x1 - im_x2) * 0.5f;

        float pow_s = re_s*re_s + im_s*im_s;
        float pow_n = re_n*re_n + im_n*im_n;
        fbf_energy += pow_s;
        bm_energy  += pow_n;

        /* 帧级相干系数：Σ|X₁·X₂*| / √(Σ|X₁|² · Σ|X₂|²) */
        sum_cross_mag += sqrtf((re_x1*re_x2 + im_x1*im_x2)*(re_x1*re_x2 + im_x1*im_x2)
                              + (im_x1*re_x2 - re_x1*im_x2)*(im_x1*re_x2 - re_x1*im_x2));
        float p1 = re_x1*re_x1 + im_x1*im_x1;
        float p2 = re_x2*re_x2 + im_x2*im_x2;
        sum_pow1 += p1;
        sum_pow2 += p2;
    }

    /* ---- 起振检测（延续现有 onset 逻辑） ---- */
    if (tf_gsc->initialized) {
        fbf_energy /= (float)n;
        if (tf_gsc->prev_fbf_energy > 1e-10f &&
            fbf_energy > tf_gsc->onset_thresh * tf_gsc->prev_fbf_energy) {
            tf_gsc->onset_frames = 8;
        }
        in_onset = (tf_gsc->onset_frames > 0) ? 1 : 0;
        tf_gsc->prev_fbf_energy = 0.8f * tf_gsc->prev_fbf_energy + 0.2f * fbf_energy;
    }

    /* ---- AMC 场景判别（仅在 initialized 后生效） ---- */
    float effective_mu = mu;  /* 默认用基准步长 */

    if (tf_gsc->initialized) {
        float spatial_snr_db = 10.0f * log10f((fbf_energy*bm_energy > 0.0f)
                                               ? (fbf_energy + 1e-10f) / (bm_energy + 1e-10f)
                                               : 1.0f);
        float coh_frame = (sum_pow1 * sum_pow2 > 1e-20f)
                          ? sum_cross_mag / (sqrtf(sum_pow1 * sum_pow2) + 1e-10f)
                          : 0.0f;
        tf_gsc->coh_smooth = 0.8f * tf_gsc->coh_smooth + 0.2f * coh_frame;
        float coh = tf_gsc->coh_smooth;

        if (spatial_snr_db < 6.0f && coh > 0.85f) {
            /* 纯静音 / 扩散噪声场 → 完全冻结 NLMS */
            effective_mu = 0.0f;
        } else if (spatial_snr_db > 12.0f) {
            /* 强单目标语音 → 正常步长 */
            effective_mu = mu;
        } else {
            /* 弱语音 / 混响 / 多声源 → 缩小步长，慢速收敛 */
            effective_mu = mu * 0.5f;
        }
    }

    tf_gsc->amc_mu = effective_mu;  /* 暴露给外部（调试用） */

    /* ============================================================ */
    /*  Phase 2: 逐 bin GSC 处理（使用 effective_mu）                   */
    /* ============================================================ */

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

        if (Y_FBF) {
            Y_FBF[k][0] = re_s;
            Y_FBF[k][1] = im_s;
        }

        /* ---- BM: Noise reference  N = (X1 - X2) / 2 ---- */
        float re_n = (re_x1 - re_x2) * 0.5f;
        float im_n = (im_x1 - im_x2) * 0.5f;

        /* ---- BM 泄漏抑制 ---- */
        float pow_s = re_s * re_s + im_s * im_s;
        float raw_pow_n = re_n * re_n + im_n * im_n;

        if (tf_gsc->S_floor[k] < 1e-15f) {
            tf_gsc->S_floor[k] = pow_s;
        } else if (pow_s < tf_gsc->S_floor[k]) {
            tf_gsc->S_floor[k] = pow_s;
        } else {
            tf_gsc->S_floor[k] += 1e-3f * tf_gsc->S_floor[k] + 1e-12f;
        }
        float s_floor = tf_gsc->S_floor[k];

        float speech_energy = pow_s - s_floor;
        if (speech_energy < 0.0f) speech_energy = 0.0f;
        float G_U = raw_pow_n / (speech_energy + raw_pow_n + TF_GSC_EPSILON);

        tf_gsc->G_smooth[k] = sf * tf_gsc->G_smooth[k] + (1.0f - sf) * G_U;
        float g = tf_gsc->G_smooth[k];

        re_n *= g;
        im_n *= g;
        float cleaned_pow_n = re_n * re_n + im_n * im_n;

        if (Y_U) {
            Y_U[k][0] = re_n;
            Y_U[k][1] = im_n;
        }

        /* ---- 读取自适应权重 W[k] ---- */
        float re_w = tf_gsc->W[k][0];
        float im_w = tf_gsc->W[k][1];

        /* ---- 首帧初始化 ---- */
        if (!tf_gsc->initialized)
        {
            tf_gsc->Pn[k] = cleaned_pow_n;
            Y[k][0] = re_s;
            Y[k][1] = im_s;
            if (k == n - 1)
                tf_gsc->initialized = 1;
            continue;
        }

        /* ---- 频点级 VAD ---- */
        int speech_present = (pow_s > beta * raw_pow_n + TF_GSC_EPSILON) ? 1 : 0;

        /* ---- GSC 输出: Y = S - W * Ñ ---- */
        float re_wn = re_w * re_n - im_w * im_n;
        float im_wn = re_w * im_n + im_w * re_n;

        float re_y = re_s - re_wn;
        float im_y = im_s - im_wn;

        /* ---- NLMS 自适应更新（AMC 控制有效步长） ---- */
        if (!speech_present)
        {
            /* Pn 始终更新（噪声功率跟踪不中断） */
            tf_gsc->Pn[k] = alpha * tf_gsc->Pn[k] + (1.0f - alpha) * cleaned_pow_n;

            if (effective_mu > 0.0f) {
                /* AMC 允许更新 → 使用有效步长 */
                float denom = tf_gsc->Pn[k] + TF_GSC_EPSILON;
                float scale = effective_mu / denom;

                float re_dw = (re_n * re_y + im_n * im_y) * scale;
                float im_dw = (re_n * im_y - im_n * re_y) * scale;

                tf_gsc->W[k][0] += re_dw;
                tf_gsc->W[k][1] += im_dw;
            } else {
                /* μ=0：完全冻结 W，仅 leak 防止噪作漂移 */
                tf_gsc->W[k][0] *= leak;
                tf_gsc->W[k][1] *= leak;
            }
        }
        else
        {
            if (in_onset) {
                tf_gsc->W[k][0] *= 0.85f;
                tf_gsc->W[k][1] *= 0.85f;
            } else {
                tf_gsc->W[k][0] *= leak;
                tf_gsc->W[k][1] *= leak;
            }
        }

        Y[k][0] = re_y;
        Y[k][1] = im_y;
    }

    /* ---- 权重幅值钳位：|W[k]| < W_max ---- */
    float W_max_sq = tf_gsc->W_max * tf_gsc->W_max;
    for (int k = 0; k < n; k++) {
        float mag_sq = tf_gsc->W[k][0]*tf_gsc->W[k][0] + tf_gsc->W[k][1]*tf_gsc->W[k][1];
        if (mag_sq > W_max_sq) {
            float scale = tf_gsc->W_max / sqrtf(mag_sq);
            tf_gsc->W[k][0] *= scale;
            tf_gsc->W[k][1] *= scale;
        }
    }

    /* ---- 语音起振保护递减 ---- */
    if (tf_gsc->onset_frames > 0)
        tf_gsc->onset_frames--;
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

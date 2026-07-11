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
    tf_gsc->Pn_bm = (float *)fftwf_malloc(sizeof(float) * complex_len);
    tf_gsc->S_floor = (float *)fftwf_malloc(sizeof(float) * complex_len);

    /* 模块 A: 扩散度掩膜阵列 */
    tf_gsc->diff_P1     = (float *)fftwf_malloc(sizeof(float) * complex_len);
    tf_gsc->diff_P2     = (float *)fftwf_malloc(sizeof(float) * complex_len);
    tf_gsc->diff_C12_re = (float *)fftwf_malloc(sizeof(float) * complex_len);
    tf_gsc->diff_C12_im = (float *)fftwf_malloc(sizeof(float) * complex_len);

    if (!tf_gsc->W || !tf_gsc->Pn || !tf_gsc->G_smooth || !tf_gsc->Pn_bm || !tf_gsc->S_floor
        || !tf_gsc->diff_P1 || !tf_gsc->diff_P2 || !tf_gsc->diff_C12_re || !tf_gsc->diff_C12_im) {
        tf_gsc_destroy(tf_gsc);
        return NULL;
    }

    /* 权重初始化为 0 (calloc 已清空) */
    memset(tf_gsc->W,  0, sizeof(fftwf_complex) * complex_len);
    memset(tf_gsc->Pn, 0, sizeof(float)          * complex_len);
    memset(tf_gsc->G_smooth, 0, sizeof(float)    * complex_len);
    memset(tf_gsc->Pn_bm,    0, sizeof(float)    * complex_len);
    memset(tf_gsc->S_floor, 0, sizeof(float)     * complex_len);
    memset(tf_gsc->diff_P1, 0, sizeof(float)     * complex_len);
    memset(tf_gsc->diff_P2, 0, sizeof(float)     * complex_len);
    memset(tf_gsc->diff_C12_re, 0, sizeof(float) * complex_len);
    memset(tf_gsc->diff_C12_im, 0, sizeof(float) * complex_len);

    tf_gsc->complex_len    = complex_len;
    tf_gsc->mu             = mu;
    tf_gsc->alpha          = alpha;
    tf_gsc->vad_thresh     = vad_thresh;
    tf_gsc->leak           = leak;
    tf_gsc->smooth_factor  = smooth_factor;
    tf_gsc->onset_thresh   = 3.0f;   /* 能量跳变 3× 触发语音起振保护 */

    /* AMC 默认值 */
    tf_gsc->amc_mu         = mu;     /* 初始 = 基准步长 */
    tf_gsc->W_max          = 0.15f;  /* 权重幅值钳位 ±0.15，限制抵消强度 */
    tf_gsc->W_min          = 0.0f;   /* 权重幅值下界 0.0=禁用，由上层按需配置 */
    tf_gsc->coh_smooth     = 0.0f;

    /* BM 泄漏保护默认：未配置时不生效 */
    tf_gsc->gU_low_bin  = 0;
    tf_gsc->gU_high_bin = 0;
    tf_gsc->gU_min      = 0.0f;

    /* 模块 A: 扩散度掩膜默认 */
    tf_gsc->diff_alpha    = 0.85f;
    tf_gsc->diff_thresh   = 0.6f;
    tf_gsc->diff_suppress = 0.25f;

    /* 模块 B: FBF 前置增强默认关闭 */
    tf_gsc->gfb_enabled = 0;

    /* 模块 C: ANC 输出钳位默认关闭 */
    tf_gsc->clamp_enabled   = 0;
    tf_gsc->clamp_min_ratio = 0.35f;
    tf_gsc->fbf_peak_db     = -100.0f;
    tf_gsc->clamp_vad       = 0;

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

        /* 模块 A: 逐频扩散度统计量平滑 (auto & cross power) */
        float cr = re_x1*re_x2 + im_x1*im_x2;
        float ci = im_x1*re_x2 - re_x1*im_x2;
        float da = tf_gsc->diff_alpha;
        tf_gsc->diff_P1[k]     = da * tf_gsc->diff_P1[k]     + (1.0f - da) * p1;
        tf_gsc->diff_P2[k]     = da * tf_gsc->diff_P2[k]     + (1.0f - da) * p2;
        tf_gsc->diff_C12_re[k] = da * tf_gsc->diff_C12_re[k] + (1.0f - da) * cr;
        tf_gsc->diff_C12_im[k] = da * tf_gsc->diff_C12_im[k] + (1.0f - da) * ci;
    }

    /* ---- 帧级 γ_avg = Σ|S|² / Σ|N|²，用于分级泄漏与 ANC 钳位门控 ---- */
    float frame_gamma = fbf_energy / (bm_energy + 1e-10f);
    float frame_leak;
    if (frame_gamma > 1.5f) {
        frame_leak = 0.995f;       /* 强语音：几乎不衰减，保收敛精度 */
    } else if (frame_gamma > 1.0f) {
        frame_leak = 0.99f;        /* 弱语音：轻微泄放，防过拟合 */
    } else {
        frame_leak = 0.98f;        /* 静音：缓慢泄放 W → 纯 FBF 输出 */
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
        } else if (spatial_snr_db < 6.0f) {
            effective_mu = 0.002f;      /* 低 SNR 非纯静音，较小步长 */
        } else if (spatial_snr_db < 12.0f) {
            float t = (spatial_snr_db - 6.0f) / 6.0f;
            effective_mu = mu * t;      /* 6~12dB 线性渐升至 μ */
        } else {
            effective_mu = mu;          /* >12dB 全速更新 */
        }
    }

    tf_gsc->amc_mu = effective_mu;  /* 暴露给外部（调试用） */

    /* ---- 模块 C: 内部 VAD 跟踪（用于 ANC 输出钳位门控） ---- */
    if (tf_gsc->clamp_enabled && tf_gsc->initialized) {
        float fbf_db = 10.0f * log10f(fbf_energy / (float)n + 1e-30f);
        if (fbf_db > tf_gsc->fbf_peak_db)
            tf_gsc->fbf_peak_db = fbf_db;
        else
            tf_gsc->fbf_peak_db -= 0.15f;  /* 6dB/s @ 50fps */
        tf_gsc->clamp_vad = (fbf_db > tf_gsc->fbf_peak_db - 25.0f) ? 1 : 0;
    } else if (!tf_gsc->initialized) {
        tf_gsc->clamp_vad = 0;
    }

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

        /* 模块 B: 前置直达语音增强 G_FB — 从源头抬高直达人声基底 */
        if (tf_gsc->gfb_enabled) {
            float denom_fb = re_x1*re_x1 + im_x1*im_x1 + re_x2*re_x2 + im_x2*im_x2 + TF_GSC_EPSILON;
            float num_fb = (re_x1+re_x2)*(re_x1+re_x2) + (im_x1+im_x2)*(im_x1+im_x2);
            float g_fb = sqrtf(num_fb / denom_fb);
            re_s *= g_fb;
            im_s *= g_fb;
        }

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

        /* BM 噪声功率帧间平滑，稳定维纳增益 */
        tf_gsc->Pn_bm[k] = 0.85f * tf_gsc->Pn_bm[k] + 0.15f * raw_pow_n;
        float smooth_pow_n = tf_gsc->Pn_bm[k];
        float G_U = smooth_pow_n / (speech_energy + smooth_pow_n + TF_GSC_EPSILON);

        tf_gsc->G_smooth[k] = sf * tf_gsc->G_smooth[k] + (1.0f - sf) * G_U;
        float g = tf_gsc->G_smooth[k];

        /* BM 泄漏抑制：语音基频段 (80-800Hz) 钳位增益 ≥ gU_min，保弱人声 */
        if (k >= tf_gsc->gU_low_bin && k <= tf_gsc->gU_high_bin && g < tf_gsc->gU_min)
            g = tf_gsc->gU_min;

        re_n *= g;
        im_n *= g;

        /* 模块 A: 时频扩散度掩膜 — 高扩散频点额外压制 BM 混响泄露 */
        {
            float P1 = tf_gsc->diff_P1[k], P2 = tf_gsc->diff_P2[k];
            float C12_re = tf_gsc->diff_C12_re[k], C12_im = tf_gsc->diff_C12_im[k];
            float denom_d = P1 * P2 + 1e-10f;
            float coh = (C12_re*C12_re + C12_im*C12_im) / denom_d;
            if (coh > 1.0f) coh = 1.0f;
            float diff = 1.0f - coh;
            float diff_g = 1.0f;
            if (diff > tf_gsc->diff_thresh) {
                diff_g = tf_gsc->diff_suppress;    /* 混响多径 → ×0.4 */
            } else if (diff > 0.3f) {
                float t = (diff - 0.3f) / (tf_gsc->diff_thresh - 0.3f);
                diff_g = 1.0f - t * (1.0f - tf_gsc->diff_suppress);
            }                                       /* diff ≤ 0.3: 直达 → 1.0 */
            re_n *= diff_g;
            im_n *= diff_g;
        }

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

        /* 模块 C: ANC 输出人声下限钳位 — VAD 门控 + 静音帧放开 */
        if (tf_gsc->clamp_enabled && tf_gsc->clamp_vad && frame_gamma > 1.0f) {
            float mag_y = sqrtf(re_y*re_y + im_y*im_y);
            float mag_s = sqrtf(re_s*re_s + im_s*im_s);
            float min_mag = tf_gsc->clamp_min_ratio * mag_s;
            if (mag_y < min_mag && mag_y > 1e-10f) {
                float scale = min_mag / mag_y;
                re_y *= scale;
                im_y *= scale;
            }
        }

        /* ---- NLMS 自适应更新（AMC 控制有效步长） ---- */
        if (!speech_present)
        {
            /* Pn 始终更新（噪声功率跟踪不中断） */
            tf_gsc->Pn[k] = alpha * tf_gsc->Pn[k] + (1.0f - alpha) * cleaned_pow_n;

            if (effective_mu > 0.0f) {
                /* 统一 NLMS 步长（不分段），噪声快速收敛，防止 pd_n 持续抬升 */
                float bin_mu = effective_mu;

                float denom = tf_gsc->Pn[k] + TF_GSC_EPSILON;
                float scale = bin_mu / denom;

                float re_dw = (re_n * re_y + im_n * im_y) * scale;
                float im_dw = (re_n * im_y - im_n * re_y) * scale;

                tf_gsc->W[k][0] += re_dw;
                tf_gsc->W[k][1] += im_dw;
            } else {
                /* μ=0：完全冻结 W，仅分级 leak 泄放权重 */
                tf_gsc->W[k][0] *= frame_leak;
                tf_gsc->W[k][1] *= frame_leak;
            }
        }
        else
        {
            if (in_onset) {
                /* 每帧 ×0.975，8 帧累计 ≈ ×0.817，从 W 末态平滑过渡无跳变 */
                tf_gsc->W[k][0] *= 0.975f;
                tf_gsc->W[k][1] *= 0.975f;
            } else {
                tf_gsc->W[k][0] *= frame_leak;
                tf_gsc->W[k][1] *= frame_leak;
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

    /* ---- 权重幅值钳位下界：|W[k]| ≥ W_min，防止过度衰减弱语音 ---- */
    float W_min_val = tf_gsc->W_min;
    if (W_min_val > 0.0f) {
        float W_min_sq = W_min_val * W_min_val;
        for (int k = 0; k < n; k++) {
            float mag_sq = tf_gsc->W[k][0]*tf_gsc->W[k][0] + tf_gsc->W[k][1]*tf_gsc->W[k][1];
            if (mag_sq < W_min_sq) {
                float scale = W_min_val / (sqrtf(mag_sq) + TF_GSC_EPSILON);
                tf_gsc->W[k][0] *= scale;
                tf_gsc->W[k][1] *= scale;
            }
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
    if (tf_gsc->Pn_bm)    fftwf_free(tf_gsc->Pn_bm);
    if (tf_gsc->S_floor)  fftwf_free(tf_gsc->S_floor);
    if (tf_gsc->diff_P1)     fftwf_free(tf_gsc->diff_P1);
    if (tf_gsc->diff_P2)     fftwf_free(tf_gsc->diff_P2);
    if (tf_gsc->diff_C12_re) fftwf_free(tf_gsc->diff_C12_re);
    if (tf_gsc->diff_C12_im) fftwf_free(tf_gsc->diff_C12_im);

    memset(tf_gsc, 0, sizeof(TfGscContext));
    free(tf_gsc);
}

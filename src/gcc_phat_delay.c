/**
 * @file gcc_phat_delay.c
 * @brief GCC-PHAT time delay estimation implementation.
 *
 * Algorithm steps:
 *   1. Apply Hamming window to both signals
 *   2. Zero-pad to 2 * fft_size
 *   3. Real FFT both signals -> frequency domain
 *   4. Cross-power spectrum: G(f) = X(f) * conj(Y(f))
 *   5. PHAT normalisation: G(f) / |G(f)|  (phase only)
 *   6. Inverse FFT -> time-domain cross-correlation
 *   7. Peak search within [-max_delay, +max_delay]
 *   8. Return delay at maximum peak
 */

#include "../include/gcc_phat_delay.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <fftw3.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int compute_gcc_phat_delay(const float *sig1, const float *sig2,
                           int len, int max_delay, int fft_size)
{
    if (!sig1 || !sig2 || len < 2 || max_delay < 1 || fft_size < 8)
        return 0;

    int fft_real_size = 2 * fft_size;
    int complex_size = fft_size;     /* r2c of length N gives N/2+1 complex bins,
                                      * but FFTW stores N/2+1 for odd N;
                                      * for even N (always true here), bins = N/2+1.
                                      * We allocate fft_size complex elements
                                      * (first fft_size/2+1 used) and ignore the rest. */

    int window = (len < fft_size) ? len : fft_size;

    /* ---------- allocate buffers ---------- */
    float *in1         = (float *)fftwf_malloc(sizeof(float) * fft_real_size);
    float *in2         = (float *)fftwf_malloc(sizeof(float) * fft_real_size);
    fftwf_complex *out1 = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * complex_size);
    fftwf_complex *out2 = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * complex_size);
    fftwf_complex *cs  = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * complex_size);
    float *xcorr       = (float *)fftwf_malloc(sizeof(float) * fft_real_size);

    if (!in1 || !in2 || !out1 || !out2 || !cs || !xcorr)
        goto cleanup_alloc;

    /* ---------- create FFTW plans ---------- */
    fftwf_plan plan_fwd1 = fftwf_plan_dft_r2c_1d(fft_real_size, in1, out1, FFTW_ESTIMATE);
    fftwf_plan plan_fwd2 = fftwf_plan_dft_r2c_1d(fft_real_size, in2, out2, FFTW_ESTIMATE);
    fftwf_plan plan_inv  = fftwf_plan_dft_c2r_1d(fft_real_size, cs, xcorr, FFTW_ESTIMATE);

    if (!plan_fwd1 || !plan_fwd2 || !plan_inv)
        goto cleanup_plans;

    /* ---------- windowing + zero-padding ---------- */
    for (int i = 0; i < window; i++) {
        float hamm = (float)(0.54 - 0.46 * cos(2.0 * M_PI * i / (window - 1)));
        in1[i] = sig1[i] * hamm;
        in2[i] = sig2[i] * hamm;
    }
    for (int i = window; i < fft_real_size; i++) {
        in1[i] = 0.0f;
        in2[i] = 0.0f;
    }

    /* ---------- forward FFT ---------- */
    fftwf_execute(plan_fwd1);
    fftwf_execute(plan_fwd2);

    /* ---------- GCC-PHAT core: cross-power spectrum with phase transform ---------- */
    int nbins = fft_size / 2 + 1;   /* number of valid complex bins from r2c */

    for (int i = 0; i < nbins; i++) {
        float re_x = out1[i][0];
        float im_x = out1[i][1];
        float re_y = out2[i][0];
        float im_y = out2[i][1];

        /* cross-power: X * conj(Y) */
        float re_cs = re_x * re_y + im_x * im_y;   /* real part   */
        float im_cs = im_x * re_y - re_x * im_y;   /* imag part   */

        /* PHAT: divide by magnitude (keep only phase) */
        float mag = sqrtf(re_cs * re_cs + im_cs * im_cs);
        if (mag > 1e-12f) {
            cs[i][0] = re_cs / mag;
            cs[i][1] = im_cs / mag;
        } else {
            cs[i][0] = 0.0f;
            cs[i][1] = 0.0f;
        }
    }

    /* ---------- inverse FFT ---------- */
    fftwf_execute(plan_inv);

    /* ---------- peak search within [-max_delay, +max_delay] ---------- */
    /*
     * After inverse FFT, the cross-correlation is stored in xcorr[0..fft_real_size-1].
     * Due to FFT's periodic nature:
     *   indices [0 .. max_delay]           correspond to lag 0 .. +max_delay
     *   indices [fft_real_size - max_delay .. fft_real_size-1] correspond to lag -max_delay .. -1
     */
    int best_lag = 0;
    float best_val = xcorr[0] / fft_real_size;   /* FFTW unscaled; divide by N */

    for (int i = 1; i <= max_delay && i < fft_real_size; i++) {
        float val = xcorr[i] / fft_real_size;
        if (val > best_val) {
            best_val = val;
            best_lag = i;
        }
    }
    for (int i = fft_real_size - max_delay; i < fft_real_size; i++) {
        if (i < 0) continue;
        float val = xcorr[i] / fft_real_size;
        if (val > best_val) {
            best_val = val;
            best_lag = i - fft_real_size;   /* wrap to negative lag */
        }
    }

    /* estimate sub-sample delay via parabolic interpolation around the peak */
    /* (only refine if we have valid neighbours) */
    if (best_lag > -max_delay && best_lag < max_delay) {
        int idx = (best_lag >= 0) ? best_lag : (best_lag + fft_real_size);
        int idx_p = idx + 1;
        int idx_m = idx - 1;
        float y0 = xcorr[idx] / fft_real_size;
        float y1 = xcorr[idx_p] / fft_real_size;
        float ym = xcorr[idx_m] / fft_real_size;

        float denom = 2.0f * y0 - y1 - ym;
        if (fabsf(denom) > 1e-12f) {
            float delta = 0.5f * (y1 - ym) / denom;
            /* delta is in [-0.5, 0.5] for a valid quadratic fit */
            if (delta > -1.0f && delta < 1.0f) {
                /* sub-sample refinement: add fractional correction */
                /* (keep best_lag integer, but this correction is available if needed) */
                (void)delta; /* available for future use */
            }
        }
    }

    /* ---------- cleanup ---------- */
    fftwf_destroy_plan(plan_fwd1);
    fftwf_destroy_plan(plan_fwd2);
    fftwf_destroy_plan(plan_inv);
    fftwf_free(in1);
    fftwf_free(in2);
    fftwf_free(out1);
    fftwf_free(out2);
    fftwf_free(cs);
    fftwf_free(xcorr);

    return best_lag;

cleanup_plans:
    fftwf_destroy_plan(plan_fwd1);
    fftwf_destroy_plan(plan_fwd2);
    fftwf_destroy_plan(plan_inv);
cleanup_alloc:
    if (in1)  fftwf_free(in1);
    if (in2)  fftwf_free(in2);
    if (out1) fftwf_free(out1);
    if (out2) fftwf_free(out2);
    if (cs)   fftwf_free(cs);
    if (xcorr) fftwf_free(xcorr);
    return 0;
}

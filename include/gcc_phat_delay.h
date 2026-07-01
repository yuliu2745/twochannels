#ifndef GCC_PHAT_DELAY_H
#define GCC_PHAT_DELAY_H

/**
 * @brief Compute time delay between two signals using GCC-PHAT.
 *
 * GCC-PHAT (Generalized Cross-Correlation with Phase Transform)
 * normalises the cross-power spectrum to retain only phase information,
 * making it robust against reverberation and noise.
 *
 * @param sig1       Reference signal (float array, length len)
 * @param sig2       Delayed signal   (float array, length len)
 * @param len        Length of both signals
 * @param max_delay  Maximum expected delay in samples (± search range)
 * @param fft_size   FFT window size (zero-padded to 2*fft_size)
 * @return Estimated delay in samples (positive  -> sig2 LAGS sig1,
 *                                     negative -> sig2 LEADS sig1,
 *                                     0         on failure)
 */
int compute_gcc_phat_delay(const float *sig1, const float *sig2,
                           int len, int max_delay, int fft_size);

#endif /* GCC_PHAT_DELAY_H */

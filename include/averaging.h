/**
 * @file averaging.h
 * @brief Signal averaging — coherent, exponential, and median filtering.
 *
 * ── Coherent (Synchronous) Averaging ─────────────────────────────
 *
 *   Trial 1:  ────/\──/\──────  (signal + noise)
 *   Trial 2:  ────/\──/\──────
 *   Trial 3:  ────/\──/\──────
 *    ...         (K trials)
 *             ─────────────────
 *   Average:  ────/\──/\──────  (noise reduced by √K)
 *
 *   SNR improvement = 10·log₁₀(K) dB
 *
 * ── Exponential Moving Average (EMA) ─────────────────────────────
 *
 *   y[n] = α·x[n] + (1-α)·y[n-1]
 *
 *   Time constant τ ≈ 1/α samples.
 *   Equivalent to IIR lowpass with pole at (1-α).
 *
 * ── Moving Median Filter ─────────────────────────────────────────
 *
 *   y[n] = median{x[n-k], ..., x[n], ..., x[n+k]}
 *
 *   Non-linear: preserves edges, rejects impulse noise.
 */

#ifndef AVERAGING_H
#define AVERAGING_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Coherent averaging of K aligned trials.
 *
 * @param trials   Array of K signal pointers, each of length n.
 * @param K        Number of trials.
 * @param n        Length of each trial.
 * @param out      Output averaged signal (length n).
 */
void coherent_average(const double **trials, int K, int n, double *out);

/**
 * @brief Exponential moving average (1st-order IIR smoother).
 *
 * @param x      Input signal (length n).
 * @param n      Signal length.
 * @param alpha  Smoothing factor (0 < alpha <= 1). Small = more smoothing.
 * @param y      Output smoothed signal (length n).
 */
void ema_filter(const double *x, int n, double alpha, double *y);

/**
 * @brief Moving average (FIR boxcar) filter.
 *
 * @param x      Input signal (length n).
 * @param n      Signal length.
 * @param M      Window size (odd recommended).
 * @param y      Output smoothed signal (length n).
 */
void moving_average(const double *x, int n, int M, double *y);

/**
 * @brief Moving median filter (non-linear, edge-preserving).
 *
 * @param x      Input signal (length n).
 * @param n      Signal length.
 * @param M      Window size (odd, typically 3-11).
 * @param y      Output filtered signal (length n).
 */
void median_filter(const double *x, int n, int M, double *y);

/**
 * @brief Compute SNR improvement from K-trial averaging.
 *
 * @param x_noisy   Noisy signal (length n).
 * @param x_clean   Clean reference (length n).
 * @param x_avg     Averaged result (length n).
 * @param n         Signal length.
 * @param snr_before  Output: SNR before averaging (dB).
 * @param snr_after   Output: SNR after averaging (dB).
 */
void compute_snr_improvement(const double *x_noisy, const double *x_clean,
                             const double *x_avg, int n,
                             double *snr_before, double *snr_after);

#ifdef __cplusplus
}
#endif

#endif /* AVERAGING_H */

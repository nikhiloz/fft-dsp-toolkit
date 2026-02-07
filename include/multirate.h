/**
 * @file multirate.h
 * @brief Multirate DSP — decimation, interpolation, and polyphase filters.
 *
 * ── Decimation by M ──────────────────────────────────────────────
 *
 *   x[n] ──► [Anti-alias LPF (fc=π/M)] ──► ↓M ──► y[n]
 *
 *   Keep every M-th sample. LPF prevents aliasing.
 *
 * ── Interpolation by L ───────────────────────────────────────────
 *
 *   x[n] ──► ↑L ──► [Anti-image LPF (fc=π/L, gain=L)] ──► y[n]
 *
 *   Insert L-1 zeros between samples, then LPF removes images.
 *
 * ── Rational Rate Conversion (L/M) ───────────────────────────────
 *
 *   x[n] ──► ↑L ──► [LPF (fc=π/max(L,M))] ──► ↓M ──► y[n]
 *
 *   Always interpolate first to avoid information loss.
 *
 * ── Polyphase Decomposition ──────────────────────────────────────
 *
 *        h[n] = h₀[n] + z⁻¹·h₁[n] + ... + z⁻⁽ᴹ⁻¹⁾·h_{M-1}[n]
 *
 *   Each sub-filter hₖ operates at the lower rate → M× savings.
 */

#ifndef MULTIRATE_H
#define MULTIRATE_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Decimate a signal by factor M with anti-alias filtering.
 *
 * Applies an internal FIR lowpass (cutoff = 0.5/M) before down-sampling.
 *
 * @param x      Input signal (length n).
 * @param n      Length of input.
 * @param M      Decimation factor (>= 1).
 * @param y      Output buffer (length >= n/M).
 * @return       Number of output samples.
 */
int decimate(const double *x, int n, int M, double *y);

/**
 * @brief Interpolate a signal by factor L with anti-image filtering.
 *
 * Inserts L-1 zeros, then applies an FIR lowpass (cutoff = 0.5/L, gain = L).
 *
 * @param x      Input signal (length n).
 * @param n      Length of input.
 * @param L      Interpolation factor (>= 1).
 * @param y      Output buffer (length >= n*L).
 * @return       Number of output samples (n * L).
 */
int interpolate(const double *x, int n, int L, double *y);

/**
 * @brief Rational sample-rate conversion by factor L/M.
 *
 * Interpolates by L, then decimates by M.
 *
 * @param x      Input signal (length n).
 * @param n      Length of input.
 * @param L      Interpolation factor.
 * @param M      Decimation factor.
 * @param y      Output buffer (length >= (n * L) / M).
 * @return       Number of output samples.
 */
int resample(const double *x, int n, int L, int M, double *y);

/**
 * @brief Polyphase decimation — efficient M-fold downsampling.
 *
 * Decomposes filter h into M polyphase branches and accumulates
 * at the lower rate. Complexity: (taps/M) multiplies per output sample.
 *
 * @param x      Input signal (length n).
 * @param n      Length of input.
 * @param h      FIR filter coefficients (length taps).
 * @param taps   Number of filter taps.
 * @param M      Decimation factor.
 * @param y      Output buffer (length >= n/M).
 * @return       Number of output samples.
 */
int polyphase_decimate(const double *x, int n,
                       const double *h, int taps, int M, double *y);

#ifdef __cplusplus
}
#endif

#endif /* MULTIRATE_H */

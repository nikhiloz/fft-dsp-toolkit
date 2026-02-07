/**
 * @file fft.h
 * @brief Cooley-Tukey Radix-2 FFT — the heart of this toolkit.
 *
 * Provides forward and inverse FFT for complex and real-valued signals.
 * All sizes must be powers of 2.
 *
 * See chapters/02-fft-fundamentals.md for the full theory walkthrough.
 */

#ifndef FFT_H
#define FFT_H

#include "dsp_utils.h"  /* Complex type */

/* ── Forward FFT ─────────────────────────────────────────────────── */

/**
 * Compute in-place complex FFT using Cooley-Tukey radix-2 DIT.
 * @param x  Array of N complex samples (modified in-place)
 * @param n  Transform size (MUST be a power of 2)
 *
 * After calling: x[k] contains the k-th frequency bin (0 ≤ k < N).
 *   x[0]   = DC component (sum of all samples)
 *   x[N/2] = Nyquist component
 */
void fft(Complex *x, int n);

/**
 * Compute FFT of a real-valued signal.
 * @param in   Real input array of length n
 * @param out  Complex output array of length n  (caller allocates)
 * @param n    Transform size (power of 2)
 */
void fft_real(const double *in, Complex *out, int n);

/* ── Inverse FFT ─────────────────────────────────────────────────── */

/**
 * Compute in-place inverse FFT (frequency → time domain).
 * @param x  Array of N complex frequency bins (modified in-place)
 * @param n  Transform size (power of 2)
 *
 * After calling: x[i] contains the i-th time-domain sample.
 */
void ifft(Complex *x, int n);

/* ── Feature extraction ──────────────────────────────────────────── */

/**
 * Extract magnitude spectrum: mag[k] = |X[k]| for k = 0..n-1
 */
void fft_magnitude(const Complex *x, double *mag, int n);

/**
 * Extract phase spectrum: phase[k] = atan2(im, re) for k = 0..n-1
 */
void fft_phase(const Complex *x, double *phase, int n);

#endif /* FFT_H */

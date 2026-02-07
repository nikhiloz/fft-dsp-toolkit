/**
 * @file hilbert.h
 * @brief Hilbert transform — analytic signal, envelope, instantaneous frequency.
 *
 * ── Analytic Signal ──────────────────────────────────────────────
 *
 *   x[n] (real)
 *     │
 *     ├──────────────────────────────► x[n]         (real part)
 *     │
 *     └──► [Hilbert FIR] ──► x̂[n]     (imaginary part)
 *
 *   z[n] = x[n] + j·x̂[n]    (analytic signal)
 *
 * ── Envelope & Instantaneous Phase ───────────────────────────────
 *
 *   Envelope:     A[n] = |z[n]| = √(x[n]² + x̂[n]²)
 *   Inst. phase:  φ[n] = atan2(x̂[n], x[n])
 *   Inst. freq:   f[n] = (1/2π) · dφ/dn
 *
 * ── Hilbert FIR Impulse Response ─────────────────────────────────
 *
 *   h[n] = { 2/(π·n)   for n odd
 *          { 0          for n even
 *
 *   (windowed to finite length and shifted for causality)
 */

#ifndef HILBERT_H
#define HILBERT_H

#include "dsp_utils.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Design a Hilbert transform FIR filter (Type III, odd length).
 *
 * @param h      Output coefficient array (length taps, must be odd).
 * @param taps   Number of filter taps (must be odd).
 */
void hilbert_design(double *h, int taps);

/**
 * @brief Compute the analytic signal using an FIR Hilbert transformer.
 *
 * @param x      Real input signal (length n).
 * @param n      Signal length.
 * @param z      Output analytic signal (Complex array, length n).
 * @param taps   Hilbert FIR order (odd, typically 31-127).
 */
void analytic_signal(const double *x, int n, Complex *z, int taps);

/**
 * @brief Compute the analytic signal via FFT (zero negative frequencies).
 *
 * Faster for long signals. Uses N-point FFT where N = next power of 2.
 *
 * @param x      Real input signal (length n).
 * @param n      Signal length.
 * @param z      Output analytic signal (Complex array, length n).
 */
void analytic_signal_fft(const double *x, int n, Complex *z);

/**
 * @brief Extract the envelope (amplitude modulation) of a signal.
 *
 * @param x      Real input signal (length n).
 * @param n      Signal length.
 * @param env    Output envelope (length n).
 * @param taps   Hilbert FIR order (odd, 0 = use FFT method).
 */
void envelope(const double *x, int n, double *env, int taps);

/**
 * @brief Compute instantaneous frequency (normalised, 0 to 0.5).
 *
 * @param x      Real input signal (length n).
 * @param n      Signal length.
 * @param freq   Output instantaneous frequency (length n).
 * @param taps   Hilbert FIR order (odd, 0 = use FFT method).
 */
void inst_frequency(const double *x, int n, double *freq, int taps);

#ifdef __cplusplus
}
#endif

#endif /* HILBERT_H */

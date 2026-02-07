/**
 * @file filter.h
 * @brief FIR (Finite Impulse Response) digital filter.
 *
 * A FIR filter computes each output sample as a weighted sum of
 * current and past input samples:
 *
 *     y[n] = h[0]·x[n] + h[1]·x[n-1] + ... + h[M-1]·x[n-M+1]
 *
 * where h[] are the filter coefficients (the "impulse response")
 * and M is the filter order.
 *
 * See chapters/04-digital-filters.md for theory, design, and examples.
 */

#ifndef FILTER_H
#define FILTER_H

/* ── Direct FIR filtering (whole-buffer) ─────────────────────────── */

/**
 * Apply FIR filter to an input signal.
 * @param in     Input signal array, length n
 * @param out    Output signal array, length n  (caller allocates)
 * @param n      Number of samples
 * @param h      Filter coefficients (impulse response), length order
 * @param order  Number of filter taps (coefficients)
 *
 * The first (order-1) output samples use zero-padding for
 * unavailable past samples.
 */
void fir_filter(const double *in, double *out, int n,
                const double *h, int order);

/* ── Simple coefficient generators ───────────────────────────────── */

/**
 * Generate a simple moving-average (boxcar) lowpass filter.
 * @param h      Output coefficient array (caller allocates, length taps)
 * @param taps   Number of taps (larger = smoother but more latency)
 *
 * Each coefficient = 1.0 / taps → uniform averaging.
 */
void fir_moving_average(double *h, int taps);

/**
 * Generate a windowed-sinc lowpass filter.
 * @param h       Output coefficient array (caller allocates, length taps)
 * @param taps    Number of taps (must be odd for symmetric filter)
 * @param cutoff  Normalized cutoff frequency (0.0 to 0.5, where 0.5 = Nyquist)
 *
 * Uses a Hamming window to control side-lobe levels.
 * See chapters/04-digital-filters.md § "Windowed-Sinc Design".
 */
void fir_lowpass(double *h, int taps, double cutoff);

#endif /* FILTER_H */

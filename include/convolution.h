/**
 * @file convolution.h
 * @brief Discrete convolution, correlation, and LTI system utilities.
 *
 * Provides linear convolution, cross-correlation, and related operations
 * fundamental to LTI system analysis.
 *
 * See chapters/04-lti-systems.md for theory and walkthrough.
 */

#ifndef CONVOLUTION_H
#define CONVOLUTION_H

#include <stddef.h>

/* ── Linear convolution ─────────────────────────────────────────── */

/**
 * Linear (direct) convolution: y[n] = Σ x[k] * h[n-k].
 *
 * Output length = x_len + h_len - 1 (caller must allocate).
 *
 * @param x      Input signal, length x_len
 * @param x_len  Length of input signal
 * @param h      Impulse response (filter), length h_len
 * @param h_len  Length of impulse response
 * @param y      Output buffer, must be at least (x_len + h_len - 1) long
 * @return       Output length (x_len + h_len - 1)
 */
int convolve(const double *x, int x_len,
             const double *h, int h_len,
             double *y);

/**
 * Causal convolution with truncation: output length equals x_len.
 * Equivalent to the first x_len samples of full convolution.
 *
 * @param x      Input signal
 * @param x_len  Length of input signal
 * @param h      Impulse response
 * @param h_len  Length of impulse response
 * @param y      Output buffer, length x_len
 */
void convolve_causal(const double *x, int x_len,
                     const double *h, int h_len,
                     double *y);

/* ── Cross-correlation ──────────────────────────────────────────── */

/**
 * Cross-correlation: r_xy[n] = Σ x[k] * y[k+n].
 *
 * Output length = x_len + y_len - 1 (same as convolution).
 * Correlation of x and y equals convolution of x with time-reversed y.
 *
 * @param x      First signal
 * @param x_len  Length of first signal
 * @param y      Second signal
 * @param y_len  Length of second signal
 * @param r      Output buffer (x_len + y_len - 1)
 * @return       Output length
 */
int cross_correlate(const double *x, int x_len,
                    const double *y, int y_len,
                    double *r);

/**
 * Auto-correlation: r_xx[n] = Σ x[k] * x[k+n].
 *
 * Output length = 2*x_len - 1, symmetric about the centre.
 *
 * @param x      Signal
 * @param x_len  Length of signal
 * @param r      Output buffer (2*x_len - 1)
 * @return       Output length
 */
int auto_correlate(const double *x, int x_len, double *r);

/* ── LTI system properties ──────────────────────────────────────── */

/**
 * Check if an impulse response represents a BIBO-stable system.
 * A DT LTI system is BIBO stable iff Σ|h[n]| < ∞.
 *
 * @param h      Impulse response
 * @param h_len  Length of impulse response
 * @return       1 if stable (finite sum), 0 otherwise
 */
int is_bibo_stable(const double *h, int h_len);

/**
 * Compute the energy of a discrete signal: E = Σ |x[n]|².
 */
double signal_energy(const double *x, int x_len);

/**
 * Compute the power of a discrete signal: P = (1/N) Σ |x[n]|².
 */
double signal_power(const double *x, int x_len);

#endif /* CONVOLUTION_H */

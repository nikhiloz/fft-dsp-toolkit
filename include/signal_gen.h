/**
 * @file signal_gen.h
 * @brief Discrete-time signal generation — impulses, steps, sinusoids, chirps, noise.
 *
 * Provides functions to generate the fundamental building-block signals
 * used throughout DSP theory and practice.
 *
 * See chapters/01-signals-and-sequences.md for theory and walkthrough.
 */

#ifndef SIGNAL_GEN_H
#define SIGNAL_GEN_H

#include <stddef.h>
#include "dsp_utils.h"  /* Complex type */

/* ── Elementary signals ──────────────────────────────────────────── */

/**
 * Unit impulse (Kronecker delta): x[n] = 1 if n == delay, else 0.
 * @param out    Output buffer of length n (caller allocates)
 * @param n      Number of samples
 * @param delay  Sample index where the impulse occurs (0-based)
 */
void gen_impulse(double *out, int n, int delay);

/**
 * Unit step: x[n] = 1 if n >= start, else 0.
 * @param out    Output buffer of length n
 * @param n      Number of samples
 * @param start  Sample index where the step begins
 */
void gen_step(double *out, int n, int start);

/**
 * Real exponential: x[n] = amplitude * base^n.
 * @param out        Output buffer of length n
 * @param n          Number of samples
 * @param amplitude  Starting amplitude
 * @param base       Decay factor (|base| < 1 for decay, > 1 for growth)
 */
void gen_exponential(double *out, int n, double amplitude, double base);

/* ── Sinusoidal signals ──────────────────────────────────────────── */

/**
 * Cosine wave: x[n] = amplitude * cos(2π * freq * n / sample_rate + phase).
 * @param out          Output buffer of length n
 * @param n            Number of samples
 * @param amplitude    Peak amplitude
 * @param freq_hz      Frequency in Hz
 * @param sample_rate  Samples per second
 * @param phase_rad    Initial phase in radians
 */
void gen_cosine(double *out, int n, double amplitude,
                double freq_hz, double sample_rate, double phase_rad);

/**
 * Sine wave: x[n] = amplitude * sin(2π * freq * n / sample_rate + phase).
 */
void gen_sine(double *out, int n, double amplitude,
              double freq_hz, double sample_rate, double phase_rad);

/**
 * Complex exponential: x[n] = amplitude * e^(j * (2π * freq * n / sample_rate + phase)).
 * @param out  Complex output buffer of length n (caller allocates)
 */
void gen_complex_exp(Complex *out, int n, double amplitude,
                     double freq_hz, double sample_rate, double phase_rad);

/* ── Composite signals ───────────────────────────────────────────── */

/**
 * Linear chirp (swept sine): frequency sweeps linearly from f0 to f1.
 * @param out          Output buffer of length n
 * @param n            Number of samples
 * @param amplitude    Peak amplitude
 * @param f0_hz        Start frequency in Hz
 * @param f1_hz        End frequency in Hz
 * @param sample_rate  Samples per second
 */
void gen_chirp(double *out, int n, double amplitude,
               double f0_hz, double f1_hz, double sample_rate);

/**
 * Sum of sinusoids: x[n] = Σ amplitudes[i] * sin(2π * freqs[i] * n / sr).
 * @param out          Output buffer of length n
 * @param n            Number of samples
 * @param freqs_hz     Array of num_tones frequencies
 * @param amplitudes   Array of num_tones amplitudes
 * @param num_tones    Number of sinusoidal components
 * @param sample_rate  Samples per second
 */
void gen_multi_tone(double *out, int n,
                    const double *freqs_hz, const double *amplitudes,
                    int num_tones, double sample_rate);

/* ── Noise ───────────────────────────────────────────────────────── */

/**
 * White noise (uniform distribution, range [-amplitude, +amplitude]).
 * @param out        Output buffer of length n
 * @param n          Number of samples
 * @param amplitude  Peak amplitude
 * @param seed       RNG seed (0 for time-based)
 */
void gen_white_noise(double *out, int n, double amplitude, unsigned int seed);

/**
 * Gaussian (normal) noise using Box-Muller transform.
 * @param out     Output buffer of length n
 * @param n       Number of samples
 * @param mean    Mean of distribution
 * @param stddev  Standard deviation
 * @param seed    RNG seed (0 for time-based)
 */
void gen_gaussian_noise(double *out, int n, double mean, double stddev,
                        unsigned int seed);

/* ── Utility ─────────────────────────────────────────────────────── */

/**
 * Add signal b into signal a: a[i] += b[i] for i = 0..n-1.
 */
void signal_add(double *a, const double *b, int n);

/**
 * Scale signal in-place: x[i] *= scale for i = 0..n-1.
 */
void signal_scale(double *x, int n, double scale);

#endif /* SIGNAL_GEN_H */

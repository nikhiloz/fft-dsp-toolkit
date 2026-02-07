/**
 * @file signal_gen.c
 * @brief Implementation of discrete-time signal generators.
 *
 * All generators write into caller-allocated buffers.
 * See chapters/01-signals-and-sequences.md for theory.
 *
 * References:
 *   - Oppenheim & Willsky, "Signals and Systems", Ch 1-2
 *   - Oppenheim & Schafer, "Discrete-Time Signal Processing", Ch 2
 *   - Lyons, "Understanding DSP", Ch 1
 */

#define _GNU_SOURCE
#include "signal_gen.h"
#include "dsp_utils.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ── Elementary signals ──────────────────────────────────────────── */

void gen_impulse(double *out, int n, int delay)
{
    memset(out, 0, (size_t)n * sizeof(double));
    if (delay >= 0 && delay < n) {
        out[delay] = 1.0;
    }
}

void gen_step(double *out, int n, int start)
{
    for (int i = 0; i < n; i++) {
        out[i] = (i >= start) ? 1.0 : 0.0;
    }
}

void gen_exponential(double *out, int n, double amplitude, double base)
{
    double val = amplitude;
    for (int i = 0; i < n; i++) {
        out[i] = val;
        val *= base;
    }
}

/* ── Sinusoidal signals ──────────────────────────────────────────── */

void gen_cosine(double *out, int n, double amplitude,
                double freq_hz, double sample_rate, double phase_rad)
{
    double omega = 2.0 * M_PI * freq_hz / sample_rate;
    for (int i = 0; i < n; i++) {
        out[i] = amplitude * cos(omega * i + phase_rad);
    }
}

void gen_sine(double *out, int n, double amplitude,
              double freq_hz, double sample_rate, double phase_rad)
{
    double omega = 2.0 * M_PI * freq_hz / sample_rate;
    for (int i = 0; i < n; i++) {
        out[i] = amplitude * sin(omega * i + phase_rad);
    }
}

void gen_complex_exp(Complex *out, int n, double amplitude,
                     double freq_hz, double sample_rate, double phase_rad)
{
    double omega = 2.0 * M_PI * freq_hz / sample_rate;
    for (int i = 0; i < n; i++) {
        double angle = omega * i + phase_rad;
        out[i] = complex_from_polar(amplitude, angle);
    }
}

/* ── Composite signals ───────────────────────────────────────────── */

void gen_chirp(double *out, int n, double amplitude,
               double f0_hz, double f1_hz, double sample_rate)
{
    /*
     * Linear chirp: instantaneous frequency sweeps from f0 to f1 over n samples.
     * Phase at sample i:  φ(i) = 2π * [f0 * t + (f1 - f0) * t² / (2 * T)]
     * where t = i / sample_rate, T = (n - 1) / sample_rate.
     */
    double T = (double)(n - 1) / sample_rate;
    for (int i = 0; i < n; i++) {
        double t = (double)i / sample_rate;
        double phase = 2.0 * M_PI * (f0_hz * t + (f1_hz - f0_hz) * t * t / (2.0 * T));
        out[i] = amplitude * sin(phase);
    }
}

void gen_multi_tone(double *out, int n,
                    const double *freqs_hz, const double *amplitudes,
                    int num_tones, double sample_rate)
{
    memset(out, 0, (size_t)n * sizeof(double));
    for (int t = 0; t < num_tones; t++) {
        double omega = 2.0 * M_PI * freqs_hz[t] / sample_rate;
        for (int i = 0; i < n; i++) {
            out[i] += amplitudes[t] * sin(omega * i);
        }
    }
}

/* ── Noise ───────────────────────────────────────────────────────── */

void gen_white_noise(double *out, int n, double amplitude, unsigned int seed)
{
    if (seed == 0) seed = (unsigned int)time(NULL);
    srand(seed);
    for (int i = 0; i < n; i++) {
        /* Uniform in [-amplitude, +amplitude] */
        out[i] = amplitude * (2.0 * ((double)rand() / RAND_MAX) - 1.0);
    }
}

void gen_gaussian_noise(double *out, int n, double mean, double stddev,
                        unsigned int seed)
{
    /*
     * Box-Muller transform: generates pairs of independent Gaussian samples
     * from uniform random samples.
     *
     * Given U1, U2 ~ Uniform(0,1):
     *   Z0 = sqrt(-2 ln U1) * cos(2π U2)
     *   Z1 = sqrt(-2 ln U1) * sin(2π U2)
     * Then Z0, Z1 ~ Normal(0,1).
     *
     * Reference: Proakis & Manolakis, Ch 12 (noise simulation).
     */
    if (seed == 0) seed = (unsigned int)time(NULL);
    srand(seed);

    for (int i = 0; i < n; i += 2) {
        double u1 = ((double)rand() + 1.0) / ((double)RAND_MAX + 1.0); /* (0,1] */
        double u2 = (double)rand() / (double)RAND_MAX;
        double r  = sqrt(-2.0 * log(u1));
        double z0 = r * cos(2.0 * M_PI * u2);
        double z1 = r * sin(2.0 * M_PI * u2);

        out[i] = mean + stddev * z0;
        if (i + 1 < n) {
            out[i + 1] = mean + stddev * z1;
        }
    }
}

/* ── Utility ─────────────────────────────────────────────────────── */

void signal_add(double *a, const double *b, int n)
{
    for (int i = 0; i < n; i++) {
        a[i] += b[i];
    }
}

void signal_scale(double *x, int n, double scale)
{
    for (int i = 0; i < n; i++) {
        x[i] *= scale;
    }
}

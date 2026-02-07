/**
 * @file dsp_utils.c
 * @brief Core DSP utilities — complex math, windows, helpers.
 *
 * TUTORIAL CROSS-REFERENCES:
 *   Complex arithmetic  → chapters/01-complex-numbers.md
 *   Window functions    → chapters/03-window-functions.md
 */

#define _GNU_SOURCE
#include "dsp_utils.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ════════════════════════════════════════════════════════════════════
 *  Complex arithmetic
 *  Tutorial ref: chapters/01-complex-numbers.md § "C Implementation"
 * ════════════════════════════════════════════════════════════════════ */

Complex complex_add(Complex a, Complex b) {
    return (Complex){ a.re + b.re, a.im + b.im };
}

Complex complex_sub(Complex a, Complex b) {
    return (Complex){ a.re - b.re, a.im - b.im };
}

/*
 * Multiply two complex numbers:
 *   (a + bi)(c + di) = (ac - bd) + (ad + bc)i
 *
 * This is the core operation inside every FFT butterfly.
 * See fft.c butterfly_radix2() for where this is used.
 */
Complex complex_mul(Complex a, Complex b) {
    return (Complex){
        a.re * b.re - a.im * b.im,
        a.re * b.im + a.im * b.re
    };
}

double complex_mag(Complex z) {
    return sqrt(z.re * z.re + z.im * z.im);
}

double complex_phase(Complex z) {
    return atan2(z.im, z.re);
}

Complex complex_from_polar(double mag, double phase) {
    return (Complex){ mag * cos(phase), mag * sin(phase) };
}

/* ════════════════════════════════════════════════════════════════════
 *  Window functions
 *  Tutorial ref: chapters/03-window-functions.md
 *
 *  Why windows?  When we take an FFT of a finite chunk of signal,
 *  the abrupt edges cause "spectral leakage" — energy smears across
 *  frequency bins.  Multiplying by a window tapers the edges to zero,
 *  trading frequency resolution for reduced leakage.
 * ════════════════════════════════════════════════════════════════════ */

/*
 * Hann window: w[i] = 0.5 * (1 - cos(2π·i / (N-1)))
 *   - Good general-purpose window
 *   - Side-lobe level: -31 dB
 *   - Main-lobe width: 4 bins
 */
double hann_window(int n, int i) {
    return 0.5 * (1.0 - cos(2.0 * M_PI * i / (n - 1)));
}

/*
 * Hamming window: w[i] = 0.54 - 0.46 * cos(2π·i / (N-1))
 *   - Similar to Hann but doesn't touch zero at edges
 *   - Side-lobe level: -42 dB (better than Hann)
 *   - Main-lobe width: 4 bins
 */
double hamming_window(int n, int i) {
    return 0.54 - 0.46 * cos(2.0 * M_PI * i / (n - 1));
}

/*
 * Blackman window: w[i] = 0.42 - 0.5·cos(2π·i/(N-1)) + 0.08·cos(4π·i/(N-1))
 *   - Excellent side-lobe suppression: -58 dB
 *   - Wider main lobe: 6 bins (poorer frequency resolution)
 *   - Best for detecting weak signals near strong ones
 */
double blackman_window(int n, int i) {
    double t = 2.0 * M_PI * i / (n - 1);
    return 0.42 - 0.5 * cos(t) + 0.08 * cos(2.0 * t);
}

void apply_window(double *signal, int n, window_fn w) {
    for (int i = 0; i < n; i++) {
        signal[i] *= w(n, i);
    }
}

/* ════════════════════════════════════════════════════════════════════
 *  Utility helpers
 * ════════════════════════════════════════════════════════════════════ */

/*
 * Round up to the next power of 2.
 * FFT requires power-of-2 lengths; this helper tells you how much
 * zero-padding is needed.
 */
int next_power_of_2(int n) {
    int p = 1;
    while (p < n) p <<= 1;
    return p;
}

/* Convert magnitude to decibels: dB = 20·log₁₀(mag) */
double db_from_magnitude(double mag) {
    if (mag <= 0.0) return -200.0;  /* floor for silence */
    return 20.0 * log10(mag);
}

/* Root-mean-square of a signal */
double rms(const double *signal, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += signal[i] * signal[i];
    }
    return sqrt(sum / n);
}

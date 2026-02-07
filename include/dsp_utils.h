/**
 * @file dsp_utils.h
 * @brief Core DSP utility functions — complex arithmetic, window functions, helpers.
 *
 * This is the foundational module. Every other module depends on it.
 * See chapters/01-complex-numbers.md for theory and walkthrough.
 */

#ifndef DSP_UTILS_H
#define DSP_UTILS_H

#include <stddef.h>

/* ── Complex number type ─────────────────────────────────────────── */

typedef struct {
    double re;  /* Real part      */
    double im;  /* Imaginary part */
} Complex;

/* ── Complex arithmetic ──────────────────────────────────────────── */

Complex complex_add(Complex a, Complex b);
Complex complex_sub(Complex a, Complex b);
Complex complex_mul(Complex a, Complex b);
double  complex_mag(Complex z);          /* |z| = sqrt(re² + im²) */
double  complex_phase(Complex z);        /* atan2(im, re) in radians */
Complex complex_from_polar(double mag, double phase);

/* ── Window functions ────────────────────────────────────────────── */
/* Each returns w[i] for a window of length n.
 * See chapters/03-window-functions.md for spectral leakage theory. */

double hann_window(int n, int i);
double hamming_window(int n, int i);
double blackman_window(int n, int i);

/* Apply a window in-place: signal[i] *= window(n, i) */
typedef double (*window_fn)(int n, int i);
void apply_window(double *signal, int n, window_fn w);

/* ── Utility helpers ─────────────────────────────────────────────── */

int    next_power_of_2(int n);
double db_from_magnitude(double mag);   /* 20 * log10(mag) */
double rms(const double *signal, int n);

#endif /* DSP_UTILS_H */

/**
 * @file hilbert.c
 * @brief Hilbert transform — analytic signal, envelope, instantaneous frequency.
 *
 * ── FIR Hilbert Filter Design ────────────────────────────────────
 *
 *   For an ideal Hilbert transformer:
 *     H(ω) = { -j   for 0 < ω < π
 *            {  j   for -π < ω < 0
 *
 *   Impulse response (infinite):
 *     h[n] =  2/(π·n)   for n odd
 *     h[n] =  0          for n even
 *
 *   We window this to get a causal FIR of length (2K+1),
 *   where the centre tap is at index K.
 *
 * ── FFT-based Method ─────────────────────────────────────────────
 *
 *   x[n] ──► FFT ──► zero negative freq ──► IFFT ──► z[n]
 *
 *   Specifically:
 *     X[k] for k=0:     multiply by 1
 *     X[k] for k=1..N/2-1: multiply by 2
 *     X[k] for k=N/2:   multiply by 1
 *     X[k] for k=N/2+1..N-1: multiply by 0
 *
 *   Result is the analytic signal z[n] = x[n] + j·x̂[n].
 */

#include "hilbert.h"
#include "fft.h"
#include "dsp_utils.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ── Hilbert FIR Design ───────────────────────────────────────── */

void hilbert_design(double *h, int taps)
{
    if (!h || taps < 3) return;
    /* Must be odd for Type III FIR */
    if ((taps & 1) == 0) taps--;

    int K = taps / 2;

    for (int i = 0; i < taps; i++) {
        int n = i - K;
        if (n == 0) {
            h[i] = 0.0;
        } else if ((n & 1) != 0) {
            /* Odd n: h[n] = 2/(π·n), windowed with Hann */
            double win = hann_window(taps, i);
            h[i] = (2.0 / (M_PI * (double)n)) * win;
        } else {
            h[i] = 0.0;
        }
    }
}

/* ── Analytic Signal (FIR) ────────────────────────────────────── */

void analytic_signal(const double *x, int n, Complex *z, int taps)
{
    if (!x || !z || n <= 0) return;
    if (taps < 3) taps = 31;
    if ((taps & 1) == 0) taps--;

    double *h = (double *)calloc((size_t)taps, sizeof(double));
    hilbert_design(h, taps);

    int K = taps / 2;

    for (int i = 0; i < n; i++) {
        z[i].re = x[i];

        /* Convolve with Hilbert filter */
        double acc = 0.0;
        for (int j = 0; j < taps; j++) {
            int idx = i - j + K;
            if (idx >= 0 && idx < n)
                acc += h[j] * x[idx];
        }
        z[i].im = acc;
    }

    free(h);
}

/* ── Next power of 2 ─────────────────────────────────────────── */

static int next_pow2(int n)
{
    int p = 1;
    while (p < n) p <<= 1;
    return p;
}

/* ── Analytic Signal (FFT) ────────────────────────────────────── */

void analytic_signal_fft(const double *x, int n, Complex *z)
{
    if (!x || !z || n <= 0) return;

    int N = next_pow2(n);
    Complex *X = (Complex *)calloc((size_t)N, sizeof(Complex));

    /* Copy real signal */
    for (int i = 0; i < n; i++) {
        X[i].re = x[i];
        X[i].im = 0.0;
    }

    /* Forward FFT */
    fft(X, N);

    /*
     * Create one-sided spectrum:
     *   k = 0:         keep (DC)
     *   k = 1..N/2-1:  multiply by 2 (positive frequencies)
     *   k = N/2:       keep (Nyquist)
     *   k = N/2+1..N-1: zero out (negative frequencies)
     */
    /* k=0: unchanged */
    for (int k = 1; k < N / 2; k++) {
        X[k].re *= 2.0;
        X[k].im *= 2.0;
    }
    /* k = N/2: unchanged */
    for (int k = N / 2 + 1; k < N; k++) {
        X[k].re = 0.0;
        X[k].im = 0.0;
    }

    /* Inverse FFT */
    ifft(X, N);

    /* Copy result, trim to original length */
    for (int i = 0; i < n; i++) {
        z[i].re = X[i].re;
        z[i].im = X[i].im;
    }

    free(X);
}

/* ── Envelope ─────────────────────────────────────────────────── */

void envelope(const double *x, int n, double *env, int taps)
{
    if (!x || !env || n <= 0) return;

    Complex *z = (Complex *)calloc((size_t)n, sizeof(Complex));

    if (taps > 0)
        analytic_signal(x, n, z, taps);
    else
        analytic_signal_fft(x, n, z);

    for (int i = 0; i < n; i++)
        env[i] = complex_mag(z[i]);

    free(z);
}

/* ── Instantaneous Frequency ──────────────────────────────────── */

void inst_frequency(const double *x, int n, double *freq, int taps)
{
    if (!x || !freq || n <= 0) return;

    Complex *z = (Complex *)calloc((size_t)n, sizeof(Complex));

    if (taps > 0)
        analytic_signal(x, n, z, taps);
    else
        analytic_signal_fft(x, n, z);

    /* Instantaneous phase via atan2, then differentiate */
    double prev_phase = atan2(z[0].im, z[0].re);
    freq[0] = 0.0;

    for (int i = 1; i < n; i++) {
        double phase = atan2(z[i].im, z[i].re);
        double dp = phase - prev_phase;

        /* Unwrap */
        while (dp > M_PI)  dp -= 2.0 * M_PI;
        while (dp < -M_PI) dp += 2.0 * M_PI;

        /* Normalised frequency (0 to 0.5) */
        freq[i] = dp / (2.0 * M_PI);
        prev_phase = phase;
    }

    free(z);
}

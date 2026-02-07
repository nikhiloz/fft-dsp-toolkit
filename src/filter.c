/**
 * @file filter.c
 * @brief FIR digital filter implementation.
 *
 * TUTORIAL CROSS-REFERENCES:
 *   FIR theory & design     → chapters/04-digital-filters.md
 *   Window for sinc design  → chapters/03-window-functions.md
 *
 * IMPLEMENTATION NOTES:
 *   The core FIR operation is a convolution (dot product of
 *   coefficients with delayed input samples):
 *
 *     y[n] = Σ h[k] · x[n - k],  k = 0 to order-1
 *
 *   This is an O(N·M) operation where N = signal length, M = order.
 *   For large filters, FFT-based convolution is faster — but this
 *   direct form is simpler to understand and sufficient for
 *   real-time per-sample processing.
 */

#define _GNU_SOURCE
#include "filter.h"
#include "dsp_utils.h"   /* hamming_window */
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ════════════════════════════════════════════════════════════════════
 *  FIR filter — direct convolution
 *
 *  For each output sample y[n]:
 *    y[n] = h[0]·x[n] + h[1]·x[n-1] + ... + h[M-1]·x[n-M+1]
 *
 *  When (n - k) < 0, we assume x was zero before the signal started.
 *  This is called "zero-padding" or "causal assumption".
 *
 *  Tutorial: see chapters/04-digital-filters.md § "Convolution"
 * ════════════════════════════════════════════════════════════════════ */

void fir_filter(const double *in, double *out, int n,
                const double *h, int order)
{
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int k = 0; k < order; k++) {
            int idx = i - k;
            if (idx >= 0) {
                sum += h[k] * in[idx];
            }
            /* else: x[idx] = 0 (zero-padding) */
        }
        out[i] = sum;
    }
}

/* ════════════════════════════════════════════════════════════════════
 *  Moving-average filter
 *
 *  The simplest lowpass filter: each coefficient = 1/N.
 *  Equivalent to averaging the last N samples.
 *
 *  Frequency response: sinc-shaped with nulls at multiples of fs/N.
 *  Not great side-lobe performance, but trivial to understand.
 *
 *  Tutorial: see chapters/04-digital-filters.md § "Moving Average"
 * ════════════════════════════════════════════════════════════════════ */

void fir_moving_average(double *h, int taps) {
    double coeff = 1.0 / taps;
    for (int i = 0; i < taps; i++) {
        h[i] = coeff;
    }
}

/* ════════════════════════════════════════════════════════════════════
 *  Windowed-sinc lowpass filter
 *
 *  The ideal lowpass impulse response is a sinc function:
 *    h_ideal[k] = sin(2π·fc·k) / (π·k)
 *
 *  But sinc extends to ±∞, so we truncate and apply a window.
 *  The Hamming window gives -42 dB side-lobe attenuation.
 *
 *  Steps:
 *    1. Generate sinc at the desired cutoff
 *    2. Apply Hamming window
 *    3. Normalize so DC gain = 1.0
 *
 *  Tutorial: see chapters/04-digital-filters.md § "Windowed-Sinc"
 * ════════════════════════════════════════════════════════════════════ */

void fir_lowpass(double *h, int taps, double cutoff) {
    int center = taps / 2;

    /* Step 1 + 2: Generate windowed sinc */
    double sum = 0.0;
    for (int i = 0; i < taps; i++) {
        int offset = i - center;
        if (offset == 0) {
            h[i] = 2.0 * cutoff;
        } else {
            /* sinc: sin(2π·fc·n) / (π·n) */
            h[i] = sin(2.0 * M_PI * cutoff * offset) / (M_PI * offset);
        }
        /* Apply Hamming window */
        h[i] *= hamming_window(taps, i);
        sum += h[i];
    }

    /* Step 3: Normalize for unity DC gain */
    if (sum != 0.0) {
        for (int i = 0; i < taps; i++) {
            h[i] /= sum;
        }
    }
}

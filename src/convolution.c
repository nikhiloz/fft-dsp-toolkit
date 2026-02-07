/**
 * @file convolution.c
 * @brief Implementation of discrete convolution, correlation, and LTI utilities.
 *
 * All functions are direct implementations — no FFT-based fast convolution here.
 * That belongs in Ch16 (Overlap-Add/Save) once we have the DFT chapter.
 */

#define _GNU_SOURCE
#include "convolution.h"
#include <math.h>
#include <string.h>

/* ── Linear convolution ─────────────────────────────────────────── */

int convolve(const double *x, int x_len,
             const double *h, int h_len,
             double *y)
{
    int y_len = x_len + h_len - 1;

    /* Zero output first */
    memset(y, 0, (size_t)y_len * sizeof(double));

    /* Direct sum: y[n] = Σ_{k=0}^{h_len-1} h[k] * x[n-k] */
    for (int n = 0; n < y_len; n++) {
        for (int k = 0; k < h_len; k++) {
            int idx = n - k;
            if (idx >= 0 && idx < x_len) {
                y[n] += h[k] * x[idx];
            }
        }
    }

    return y_len;
}

void convolve_causal(const double *x, int x_len,
                     const double *h, int h_len,
                     double *y)
{
    /* Only compute first x_len samples of full convolution */
    for (int n = 0; n < x_len; n++) {
        double sum = 0.0;
        for (int k = 0; k < h_len; k++) {
            int idx = n - k;
            if (idx >= 0 && idx < x_len) {
                sum += h[k] * x[idx];
            }
        }
        y[n] = sum;
    }
}

/* ── Cross-correlation ──────────────────────────────────────────── */

int cross_correlate(const double *x, int x_len,
                    const double *y, int y_len,
                    double *r)
{
    int r_len = x_len + y_len - 1;
    memset(r, 0, (size_t)r_len * sizeof(double));

    /*
     * r_xy[lag] = Σ x[k] * y[k + lag]
     * where lag ranges from -(y_len-1) to +(x_len-1).
     * We store r[0] = lag -(y_len-1), ..., r[y_len-1] = lag 0, etc.
     */
    for (int lag_idx = 0; lag_idx < r_len; lag_idx++) {
        int lag = lag_idx - (y_len - 1);
        double sum = 0.0;
        for (int k = 0; k < x_len; k++) {
            int j = k + lag;
            if (j >= 0 && j < y_len) {
                sum += x[k] * y[j];
            }
        }
        r[lag_idx] = sum;
    }

    return r_len;
}

int auto_correlate(const double *x, int x_len, double *r)
{
    return cross_correlate(x, x_len, x, x_len, r);
}

/* ── LTI system properties ──────────────────────────────────────── */

int is_bibo_stable(const double *h, int h_len)
{
    double sum = 0.0;
    for (int i = 0; i < h_len; i++) {
        sum += fabs(h[i]);
    }
    /* Finite-length impulse responses are always BIBO stable
       as long as the sum is finite (not NaN/Inf) */
    return isfinite(sum) ? 1 : 0;
}

double signal_energy(const double *x, int x_len)
{
    double e = 0.0;
    for (int i = 0; i < x_len; i++) {
        e += x[i] * x[i];
    }
    return e;
}

double signal_power(const double *x, int x_len)
{
    if (x_len <= 0) return 0.0;
    return signal_energy(x, x_len) / (double)x_len;
}

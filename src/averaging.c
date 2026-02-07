/**
 * @file averaging.c
 * @brief Signal averaging — coherent, exponential, moving, and median.
 *
 * ── Coherent Averaging ───────────────────────────────────────────
 *
 *   out[n] = (1/K) · Σ_{k=0}^{K-1} trial_k[n]
 *
 *   Noise (zero-mean, uncorrelated) averages to zero.
 *   Signal (deterministic, repeatable) reinforces.
 *   SNR improvement ≈ 10·log₁₀(K) dB.
 *
 * ── EMA (Exponential Moving Average) ─────────────────────────────
 *
 *   y[n] = α·x[n] + (1-α)·y[n-1]
 *
 *   Equivalent to 1st-order IIR:
 *     H(z) = α / (1 - (1-α)·z⁻¹)
 *     pole at z = (1-α), stable for 0 < α ≤ 1
 *
 * ── Median Filter ────────────────────────────────────────────────
 *
 *   Sorts window, picks middle value.
 *   O(M·log M) per sample (insertion sort for small M).
 *
 *   Key property: removes impulse noise while preserving edges.
 */

#include "averaging.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ── Coherent Averaging ───────────────────────────────────────── */

void coherent_average(const double **trials, int K, int n, double *out)
{
    if (!trials || !out || K <= 0 || n <= 0) return;

    memset(out, 0, (size_t)n * sizeof(double));

    for (int k = 0; k < K; k++) {
        if (!trials[k]) continue;
        for (int i = 0; i < n; i++)
            out[i] += trials[k][i];
    }

    double inv_K = 1.0 / (double)K;
    for (int i = 0; i < n; i++)
        out[i] *= inv_K;
}

/* ── Exponential Moving Average ───────────────────────────────── */

void ema_filter(const double *x, int n, double alpha, double *y)
{
    if (!x || !y || n <= 0) return;
    if (alpha <= 0.0) alpha = 0.01;
    if (alpha > 1.0) alpha = 1.0;

    y[0] = x[0];
    for (int i = 1; i < n; i++)
        y[i] = alpha * x[i] + (1.0 - alpha) * y[i - 1];
}

/* ── Moving Average (boxcar FIR) ──────────────────────────────── */

void moving_average(const double *x, int n, int M, double *y)
{
    if (!x || !y || n <= 0 || M <= 0) return;

    int half = M / 2;

    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        int count = 0;
        for (int j = i - half; j <= i + half; j++) {
            if (j >= 0 && j < n) {
                sum += x[j];
                count++;
            }
        }
        y[i] = sum / (double)count;
    }
}

/* ── Median Filter ────────────────────────────────────────────── */

/* Insertion sort for small windows */
static void insertion_sort(double *arr, int n)
{
    for (int i = 1; i < n; i++) {
        double key = arr[i];
        int j = i - 1;
        while (j >= 0 && arr[j] > key) {
            arr[j + 1] = arr[j];
            j--;
        }
        arr[j + 1] = key;
    }
}

void median_filter(const double *x, int n, int M, double *y)
{
    if (!x || !y || n <= 0 || M <= 0) return;
    if ((M & 1) == 0) M++;  /* Ensure odd */

    int half = M / 2;
    double *window = (double *)malloc((size_t)M * sizeof(double));

    for (int i = 0; i < n; i++) {
        int count = 0;
        for (int j = i - half; j <= i + half; j++) {
            if (j >= 0 && j < n)
                window[count++] = x[j];
        }
        insertion_sort(window, count);
        y[i] = window[count / 2];
    }

    free(window);
}

/* ── SNR Improvement Computation ──────────────────────────────── */

void compute_snr_improvement(const double *x_noisy, const double *x_clean,
                             const double *x_avg, int n,
                             double *snr_before, double *snr_after)
{
    if (!x_noisy || !x_clean || !x_avg || n <= 0) return;

    double sig_pow = 0.0, noise_before = 0.0, noise_after = 0.0;

    for (int i = 0; i < n; i++) {
        sig_pow += x_clean[i] * x_clean[i];
        double nb = x_noisy[i] - x_clean[i];
        double na = x_avg[i] - x_clean[i];
        noise_before += nb * nb;
        noise_after  += na * na;
    }

    if (snr_before) {
        *snr_before = (noise_before > 0.0)
            ? 10.0 * log10(sig_pow / noise_before) : 100.0;
    }
    if (snr_after) {
        *snr_after = (noise_after > 0.0)
            ? 10.0 * log10(sig_pow / noise_after) : 100.0;
    }
}

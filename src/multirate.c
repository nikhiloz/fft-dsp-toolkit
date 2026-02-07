/**
 * @file multirate.c
 * @brief Multirate DSP — decimation, interpolation, polyphase filtering.
 *
 * ── Decimation Pipeline ──────────────────────────────────────────
 *
 *   x[n] ──► [FIR LPF fc=0.5/M] ──► keep every M-th ──► y[m]
 *
 *   Anti-alias filter prevents spectral folding after down-sampling.
 *   Output length = floor(n / M).
 *
 * ── Interpolation Pipeline ───────────────────────────────────────
 *
 *   x[n] ──► insert L-1 zeros ──► [FIR LPF fc=0.5/L, gain=L] ──► y[m]
 *
 *   Anti-image filter removes spectral images from zero-insertion.
 *   Output length = n * L.
 *
 * ── Polyphase Decimator ──────────────────────────────────────────
 *
 *   h[n] split into M sub-filters:
 *
 *     x[n] ──┬──► [h₀] ──┐
 *             ├──► [h₁] ──┤
 *             ├──► [h₂] ──┼──► Σ ──► y[m]
 *             └──► [hₘ₋₁]─┘
 *
 *   Each hₖ processes at rate 1/M → M× fewer multiplies.
 */

#include "multirate.h"
#include "filter.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Internal anti-alias filter taps (auto-sized) */
#define ANTI_ALIAS_TAPS_PER_FACTOR 8
#define MIN_TAPS 15
#define MAX_TAPS 127

static int compute_filter_len(int factor)
{
    int taps = factor * ANTI_ALIAS_TAPS_PER_FACTOR + 1;
    if (taps < MIN_TAPS) taps = MIN_TAPS;
    if (taps > MAX_TAPS) taps = MAX_TAPS;
    /* Ensure odd */
    if ((taps & 1) == 0) taps++;
    return taps;
}

/* ── Decimation ───────────────────────────────────────────────── */

int decimate(const double *x, int n, int M, double *y)
{
    if (M <= 0 || n <= 0 || !x || !y) return 0;
    if (M == 1) {
        memcpy(y, x, (size_t)n * sizeof(double));
        return n;
    }

    /* Design anti-alias lowpass */
    int taps = compute_filter_len(M);
    double *h = (double *)calloc((size_t)taps, sizeof(double));
    double cutoff = 0.5 / (double)M;
    fir_lowpass(h, taps, cutoff);

    /* Filter, then downsample */
    double *filtered = (double *)calloc((size_t)n, sizeof(double));
    fir_filter(x, filtered, n, h, taps);

    int out_len = 0;
    for (int i = 0; i < n; i += M)
        y[out_len++] = filtered[i];

    free(h);
    free(filtered);
    return out_len;
}

/* ── Interpolation ────────────────────────────────────────────── */

int interpolate(const double *x, int n, int L, double *y)
{
    if (L <= 0 || n <= 0 || !x || !y) return 0;
    if (L == 1) {
        memcpy(y, x, (size_t)n * sizeof(double));
        return n;
    }

    int out_len = n * L;

    /* Zero-insert: place original samples at multiples of L */
    double *upsampled = (double *)calloc((size_t)out_len, sizeof(double));
    for (int i = 0; i < n; i++)
        upsampled[i * L] = x[i];

    /* Design anti-image lowpass with gain = L */
    int taps = compute_filter_len(L);
    double *h = (double *)calloc((size_t)taps, sizeof(double));
    double cutoff = 0.5 / (double)L;
    fir_lowpass(h, taps, cutoff);

    /* Scale by L to compensate for zero-inserted energy loss */
    for (int i = 0; i < taps; i++)
        h[i] *= (double)L;

    fir_filter(upsampled, y, out_len, h, taps);

    free(upsampled);
    free(h);
    return out_len;
}

/* ── Rational Resampling ──────────────────────────────────────── */

int resample(const double *x, int n, int L, int M, double *y)
{
    if (L <= 0 || M <= 0 || n <= 0 || !x || !y) return 0;

    /* Interpolate by L first */
    int interp_len = n * L;
    double *tmp = (double *)calloc((size_t)interp_len, sizeof(double));
    interpolate(x, n, L, tmp);

    /* Then decimate by M (no extra filter — interpolation filter
     * already limits bandwidth to min(1/L, 1/M) if L >= M) */
    int out_len = 0;
    if (L >= M) {
        /* Already bandwidth-limited by interpolation filter */
        for (int i = 0; i < interp_len; i += M)
            y[out_len++] = tmp[i];
    } else {
        /* Need additional anti-alias filter */
        out_len = decimate(tmp, interp_len, M, y);
    }

    free(tmp);
    return out_len;
}

/* ── Polyphase Decimation ─────────────────────────────────────── */

int polyphase_decimate(const double *x, int n,
                       const double *h, int taps, int M, double *y)
{
    if (M <= 0 || n <= 0 || taps <= 0 || !x || !h || !y) return 0;

    /*
     * Polyphase decomposition:
     *   sub-filter k: hₖ[m] = h[k + m·M]   for m = 0, 1, ...
     *   sub-filter length = ceil(taps / M)
     */
    int sub_len = (taps + M - 1) / M;
    int out_len = n / M;

    for (int out_idx = 0; out_idx < out_len; out_idx++) {
        double acc = 0.0;
        int base = out_idx * M;  /* Input index for this output */

        for (int k = 0; k < M; k++) {
            /* Sub-filter k, applied to x[base - k] (time-reversed) */
            for (int m = 0; m < sub_len; m++) {
                int h_idx = k + m * M;
                int x_idx = base - k - m * M;
                if (h_idx < taps && x_idx >= 0 && x_idx < n)
                    acc += h[h_idx] * x[x_idx];
            }
        }
        y[out_idx] = acc;
    }

    return out_len;
}

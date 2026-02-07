/**
 * @file spectrum.c
 * @brief Power spectral density estimation — periodogram & Welch's method.
 *
 * ── Welch's Method Pipeline ─────────────────────────────────────
 *
 *   Input x[N]
 *     │
 *     ▼
 *   ┌─────────────────────────────────────────────┐
 *   │  Segment into overlapping blocks             │
 *   │                                              │
 *   │  seg[0]: x[0 .. L-1]                        │
 *   │  seg[1]: x[D .. D+L-1]     D = L - overlap  │
 *   │  seg[2]: x[2D .. 2D+L-1]                    │
 *   │  ...                                         │
 *   └─────────────────────────────────────────────┘
 *     │
 *     ▼  (for each segment)
 *   ┌───────────┐   ┌───────────┐   ┌───────────┐
 *   │  Window   │──►│    FFT    │──►│  |X[k]|²  │
 *   │  (Hann)   │   │  (NFFT)   │   │  / (N·U)  │
 *   └───────────┘   └───────────┘   └───────────┘
 *     │
 *     ▼
 *   Average all segment periodograms → PSD estimate
 *
 *   U = (1/L) Σ w[n]²   (window power normalisation)
 *
 * ── Periodogram (single segment, no averaging) ─────────────────
 *
 *   x[N] ──► window ──► FFT ──► |X[k]|² / N  =  PSD
 *
 *   Simple but high variance (no averaging).
 */
#define _POSIX_C_SOURCE 200809L
#include "spectrum.h"
#include "fft.h"
#include "dsp_utils.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/*  Helpers                                                           */
/* ------------------------------------------------------------------ */

/** Check that n is a power of 2 */
static int is_power_of_2(int n)
{
    return n > 0 && (n & (n - 1)) == 0;
}

/**
 *  Compute |X[k]|² for k = 0 … nfft/2.
 *  buf must be of length nfft (Complex).  Modified in place.
 */
static void accumulate_power(Complex *buf, int nfft,
                             double *psd, double scale, int accumulate)
{
    int n_bins = nfft / 2 + 1;
    for (int k = 0; k < n_bins; k++) {
        double power = (buf[k].re * buf[k].re + buf[k].im * buf[k].im) * scale;
        if (accumulate)
            psd[k] += power;
        else
            psd[k] = power;
    }
}

/* ------------------------------------------------------------------ */
/*  periodogram                                                       */
/* ------------------------------------------------------------------ */

int periodogram(const double *x, int n, double *psd, int nfft)
{
    return periodogram_windowed(x, n, psd, nfft, NULL);
}

int periodogram_windowed(const double *x, int n, double *psd, int nfft,
                         window_fn win)
{
    if (!x || !psd || n <= 0 || nfft <= 0 || !is_power_of_2(nfft))
        return -1;
    if (nfft < n) return -1;

    int n_bins = nfft / 2 + 1;

    /* Allocate working buffer */
    Complex *buf = (Complex *)calloc((size_t)nfft, sizeof(Complex));
    if (!buf) return -1;

    /* Copy signal with optional window */
    double win_power = 0.0;
    for (int i = 0; i < n; i++) {
        double w = win ? win(n, i) : 1.0;
        win_power += w * w;
        buf[i].re = x[i] * w;
        buf[i].im = 0.0;
    }
    /* zero-pad (calloc already zeroed) */

    if (win_power < 1e-30) win_power = (double)n;   /* safety for rectangular */

    /* FFT */
    fft(buf, nfft);

    /* Power spectrum: |X|² / (win_power) */
    double scale = 1.0 / win_power;
    accumulate_power(buf, nfft, psd, scale, 0);

    /* One-sided: double non-DC, non-Nyquist bins */
    for (int k = 1; k < n_bins - 1; k++)
        psd[k] *= 2.0;

    free(buf);
    return n_bins;
}

/* ------------------------------------------------------------------ */
/*  Welch PSD                                                         */
/* ------------------------------------------------------------------ */

int welch_psd(const double *x, int n, double *psd, int nfft,
              int seg_len, int overlap, window_fn win)
{
    if (!x || !psd || n <= 0 || nfft <= 0 || !is_power_of_2(nfft))
        return -1;
    if (seg_len <= 0 || seg_len > nfft || overlap < 0 || overlap >= seg_len)
        return -1;

    int n_bins = nfft / 2 + 1;
    int hop    = seg_len - overlap;
    int n_segs = 0;

    /* Pre-compute window and its power */
    double *w = (double *)malloc((size_t)seg_len * sizeof(double));
    Complex *buf = (Complex *)calloc((size_t)nfft, sizeof(Complex));
    if (!w || !buf) {
        free(w);
        free(buf);
        return -1;
    }

    double win_power = 0.0;
    for (int i = 0; i < seg_len; i++) {
        w[i] = win ? win(seg_len, i) : 1.0;
        win_power += w[i] * w[i];
    }
    if (win_power < 1e-30) win_power = (double)seg_len;

    double scale = 1.0 / win_power;

    /* Zero the accumulator */
    memset(psd, 0, (size_t)n_bins * sizeof(double));

    /* Iterate over segments */
    for (int start = 0; start + seg_len <= n; start += hop) {
        /* Window the segment into buf */
        memset(buf, 0, (size_t)nfft * sizeof(Complex));
        for (int i = 0; i < seg_len; i++) {
            buf[i].re = x[start + i] * w[i];
            buf[i].im = 0.0;
        }

        /* FFT */
        fft(buf, nfft);

        /* Accumulate power */
        accumulate_power(buf, nfft, psd, scale, 1);
        n_segs++;
    }

    if (n_segs == 0) {
        free(w);
        free(buf);
        return -1;
    }

    /* Average and make one-sided */
    double inv_segs = 1.0 / (double)n_segs;
    for (int k = 0; k < n_bins; k++)
        psd[k] *= inv_segs;
    for (int k = 1; k < n_bins - 1; k++)
        psd[k] *= 2.0;

    free(w);
    free(buf);
    return n_segs;
}

/* ------------------------------------------------------------------ */
/*  Cross PSD                                                         */
/* ------------------------------------------------------------------ */

int cross_psd(const double *x, const double *y, int n,
              Complex *cpsd, int nfft,
              int seg_len, int overlap, window_fn win)
{
    if (!x || !y || !cpsd || n <= 0 || nfft <= 0 || !is_power_of_2(nfft))
        return -1;
    if (seg_len <= 0 || seg_len > nfft || overlap < 0 || overlap >= seg_len)
        return -1;

    int n_bins = nfft / 2 + 1;
    int hop    = seg_len - overlap;
    int n_segs = 0;

    double *w = (double *)malloc((size_t)seg_len * sizeof(double));
    Complex *bx = (Complex *)calloc((size_t)nfft, sizeof(Complex));
    Complex *by = (Complex *)calloc((size_t)nfft, sizeof(Complex));
    if (!w || !bx || !by) {
        free(w); free(bx); free(by);
        return -1;
    }

    double win_power = 0.0;
    for (int i = 0; i < seg_len; i++) {
        w[i] = win ? win(seg_len, i) : 1.0;
        win_power += w[i] * w[i];
    }
    if (win_power < 1e-30) win_power = (double)seg_len;

    double scale = 1.0 / win_power;

    memset(cpsd, 0, (size_t)n_bins * sizeof(Complex));

    for (int start = 0; start + seg_len <= n; start += hop) {
        memset(bx, 0, (size_t)nfft * sizeof(Complex));
        memset(by, 0, (size_t)nfft * sizeof(Complex));
        for (int i = 0; i < seg_len; i++) {
            bx[i].re = x[start + i] * w[i];
            by[i].re = y[start + i] * w[i];
        }

        fft(bx, nfft);
        fft(by, nfft);

        /* Pxy += conj(X) · Y */
        for (int k = 0; k < n_bins; k++) {
            double xr = bx[k].re, xi = bx[k].im;
            double yr = by[k].re, yi = by[k].im;
            cpsd[k].re += (xr * yr + xi * yi) * scale;   /* Re(conj(X)·Y) */
            cpsd[k].im += (xr * yi - xi * yr) * scale;   /* Im(conj(X)·Y) */
        }
        n_segs++;
    }

    if (n_segs == 0) {
        free(w); free(bx); free(by);
        return -1;
    }

    double inv = 1.0 / (double)n_segs;
    for (int k = 0; k < n_bins; k++) {
        cpsd[k].re *= inv;
        cpsd[k].im *= inv;
    }
    /* One-sided doubling */
    for (int k = 1; k < n_bins - 1; k++) {
        cpsd[k].re *= 2.0;
        cpsd[k].im *= 2.0;
    }

    free(w); free(bx); free(by);
    return n_segs;
}

/* ------------------------------------------------------------------ */
/*  Utility                                                           */
/* ------------------------------------------------------------------ */

void psd_to_db(const double *psd, double *psd_db, int n_bins, double floor_db)
{
    for (int k = 0; k < n_bins; k++) {
        double val = 10.0 * log10(psd[k] + 1e-30);
        psd_db[k] = val < floor_db ? floor_db : val;
    }
}

void psd_freq_axis(double *freq, int n_bins, double fs)
{
    int nfft = (n_bins - 1) * 2;
    for (int k = 0; k < n_bins; k++)
        freq[k] = (double)k * fs / (double)nfft;
}

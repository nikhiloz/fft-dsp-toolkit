/**
 * @file remez.c
 * @brief Optimal equiripple FIR design via iteratively reweighted least-squares.
 *
 * ── Method ───────────────────────────────────────────────────────
 *
 *   For a Type I symmetric FIR of length N = 2R + 1, the zero-phase
 *   frequency response is a cosine polynomial:
 *
 *     A(f) = Σ_{k=0}^{R} a[k] · cos(2πfk)
 *
 *   where a[0] = h[R], a[k] = 2·h[R−k] for k = 1..R.
 *
 *   We minimise the weighted maximum error:
 *
 *     min_a max_f { W(f) · |D(f) − A(f)| }
 *
 *   using IRLS: start with uniform weights, solve weighted
 *   least-squares, then increase weights where error is large.
 *   This converges to the minimax (equiripple) solution.
 *
 * ── Normal Equations ─────────────────────────────────────────────
 *
 *   Cᵀ·W·C · a = Cᵀ·W · d
 *
 *   C_mk = cos(2π·f_m·k), W = diag(weights), d = desired response.
 *   Solved via Cholesky-like Gauss elimination (R+1 × R+1 system).
 */

#include "remez.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define GRID_PTS_PER_TAP 20
#define IRLS_ITERATIONS  30

/* ── Solve Ax = b by Gauss elimination (small dense system) ──── */

static int solve_linear(double *A, double *b, int n)
{
    /* In-place: A is n×n row-major, b is n×1 */
    for (int col = 0; col < n; col++) {
        /* Partial pivot */
        int best = col;
        for (int row = col + 1; row < n; row++) {
            if (fabs(A[row * n + col]) > fabs(A[best * n + col]))
                best = row;
        }
        if (best != col) {
            for (int k = 0; k < n; k++) {
                double t = A[col * n + k];
                A[col * n + k] = A[best * n + k];
                A[best * n + k] = t;
            }
            double t = b[col]; b[col] = b[best]; b[best] = t;
        }
        if (fabs(A[col * n + col]) < 1e-30) return -1;  /* singular */

        /* Eliminate below */
        for (int row = col + 1; row < n; row++) {
            double f = A[row * n + col] / A[col * n + col];
            for (int k = col; k < n; k++)
                A[row * n + k] -= f * A[col * n + k];
            b[row] -= f * b[col];
        }
    }
    /* Back-substitution */
    for (int row = n - 1; row >= 0; row--) {
        for (int k = row + 1; k < n; k++)
            b[row] -= A[row * n + k] * b[k];
        b[row] /= A[row * n + row];
    }
    return 0;
}

/* ── Weighted Least-Squares FIR Design ────────────────────────── */

static void wls_design(const double *freq, const double *des,
                       const double *wt, int ng, int R,
                       double *a /* R+1 cosine coefficients */)
{
    int nc = R + 1;

    /* Build Cᵀ·W·C and Cᵀ·W·d */
    double *CWC = (double *)calloc((size_t)(nc * nc), sizeof(double));
    double *CWd = (double *)calloc((size_t)nc, sizeof(double));

    for (int m = 0; m < ng; m++) {
        double w2 = wt[m] * wt[m];  /* squared weight for LS */
        /* Compute cosine basis row for this grid point */
        double cos_val[256];  /* R+1 values, R < 256 */
        for (int k = 0; k < nc && k < 256; k++)
            cos_val[k] = cos(2.0 * M_PI * freq[m] * (double)k);

        for (int i = 0; i < nc; i++) {
            CWd[i] += w2 * cos_val[i] * des[m];
            for (int j = i; j < nc; j++)
                CWC[i * nc + j] += w2 * cos_val[i] * cos_val[j];
        }
    }
    /* Symmetrise */
    for (int i = 0; i < nc; i++)
        for (int j = 0; j < i; j++)
            CWC[i * nc + j] = CWC[j * nc + i];

    /* Solve */
    solve_linear(CWC, CWd, nc);

    for (int k = 0; k < nc; k++)
        a[k] = CWd[k];

    free(CWC); free(CWd);
}

/* ── Main Equiripple FIR Design ───────────────────────────────── */

int remez_fir(double *h, int taps, const RemezBand *bands,
              int n_bands, int max_iter)
{
    if (!h || taps < 3 || !bands || n_bands < 1) return -1;
    if ((taps & 1) == 0) taps--;  /* ensure odd for Type I */

    int R = (taps - 1) / 2;

    /* ── Build dense frequency grid across all bands ─────────── */
    int ng = 0;
    /* Count grid points */
    for (int b = 0; b < n_bands; b++) {
        double bw = bands[b].high - bands[b].low;
        int np = (int)(bw * taps * GRID_PTS_PER_TAP) + 4;
        ng += np;
    }

    double *freq = (double *)malloc((size_t)ng * sizeof(double));
    double *des  = (double *)malloc((size_t)ng * sizeof(double));
    double *wt   = (double *)malloc((size_t)ng * sizeof(double));
    double *irls_wt = (double *)malloc((size_t)ng * sizeof(double));
    int idx = 0;

    for (int b = 0; b < n_bands; b++) {
        double lo = bands[b].low, hi = bands[b].high;
        double bw = hi - lo;
        int np = (int)(bw * taps * GRID_PTS_PER_TAP) + 4;
        for (int i = 0; i < np; i++) {
            freq[idx] = lo + bw * (double)i / (double)(np - 1);
            des[idx]  = bands[b].desired;
            wt[idx]   = bands[b].weight;
            irls_wt[idx] = bands[b].weight;
            idx++;
        }
    }
    ng = idx;

    /* ── IRLS Iterations ─────────────────────────────────────── */
    double *a = (double *)calloc((size_t)(R + 1), sizeof(double));
    double *err_buf = (double *)malloc((size_t)ng * sizeof(double));

    if (max_iter <= 0) max_iter = IRLS_ITERATIONS;

    for (int iter = 0; iter < max_iter; iter++) {
        /* Weighted least-squares design */
        wls_design(freq, des, irls_wt, ng, R, a);

        /* Evaluate error and update IRLS weights */
        double max_err = 0.0;
        for (int m = 0; m < ng; m++) {
            double Af = 0.0;
            for (int k = 0; k <= R; k++)
                Af += a[k] * cos(2.0 * M_PI * freq[m] * (double)k);
            err_buf[m] = des[m] - Af;
            double we = fabs(wt[m] * err_buf[m]);
            if (we > max_err) max_err = we;
        }

        if (max_err < 1e-15) break;

        /* IRLS reweighting: increase weight where error is large */
        for (int m = 0; m < ng; m++) {
            double we = fabs(wt[m] * err_buf[m]);
            double ratio = we / (max_err + 1e-30);
            /* Lawson's iteration: multiply weight by |error|^p */
            irls_wt[m] = wt[m] * (1.0 + 2.0 * ratio);
        }
    }

    /* ── Convert cosine coefficients to impulse response ─────── */
    int mid = R;
    h[mid] = a[0];
    for (int k = 1; k <= R; k++) {
        h[mid - k] = a[k] / 2.0;
        h[mid + k] = a[k] / 2.0;
    }

    free(a); free(err_buf); free(freq); free(des); free(wt); free(irls_wt);
    return 0;
}

/* ── Convenience: Lowpass ─────────────────────────────────────── */

int remez_lowpass(double *h, int taps, double fpass, double fstop,
                  double wpass, double wstop)
{
    RemezBand bands[2] = {
        { 0.0,   fpass, 1.0, wpass },
        { fstop, 0.5,   0.0, wstop }
    };
    return remez_fir(h, taps, bands, 2, IRLS_ITERATIONS);
}

/* ── Convenience: Bandpass ────────────────────────────────────── */

int remez_bandpass(double *h, int taps,
                   double fstop1, double fpass1,
                   double fpass2, double fstop2)
{
    RemezBand bands[3] = {
        { 0.0,    fstop1, 0.0, 1.0 },
        { fpass1, fpass2, 1.0, 1.0 },
        { fstop2, 0.5,    0.0, 1.0 }
    };
    return remez_fir(h, taps, bands, 3, IRLS_ITERATIONS);
}

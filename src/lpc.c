/**
 * @file lpc.c
 * @brief Linear Prediction Coding — Levinson-Durbin, AR spectral envelope.
 *
 * ── Autocorrelation Method Flow ──────────────────────────────────
 *
 *   x[n] ──► Autocorrelation ──► Levinson-Durbin ──► a[1..p], E
 *                r[0..p]                  │
 *                                         ├──► Residual e[n]
 *                                         └──► AR Spectrum S(f)
 *
 * ── Levinson-Durbin Step Detail ──────────────────────────────────
 *
 *   m=1: k₁ = -r[1]/r[0]
 *         a₁[1] = k₁
 *         E₁ = (1-k₁²)·r[0]
 *
 *   m=2: k₂ = -(r[2] + a₁[1]·r[1]) / E₁
 *         a₂[1] = a₁[1] + k₂·a₁[1]
 *         a₂[2] = k₂
 *         E₂ = (1-k₂²)·E₁
 *
 *   Each stage reuses previous coefficients — O(p²) total.
 */

#include "lpc.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ================================================================== */
/*  Autocorrelation                                                    */
/* ================================================================== */

void lpc_autocorrelation(const double *x, int n, double *r, int p)
{
    for (int lag = 0; lag <= p; lag++) {
        double sum = 0.0;
        for (int i = 0; i < n - lag; i++)
            sum += x[i] * x[i + lag];
        r[lag] = sum;
    }
}

/* ================================================================== */
/*  Levinson-Durbin Recursion                                          */
/* ================================================================== */

int levinson_durbin(const double *r, int p,
                    double *a, double *k_out, double *E_out)
{
    if (r[0] <= 0.0) return -1;

    double *a_prev = (double *)calloc((size_t)(p + 1), sizeof(double));
    double *a_curr = (double *)calloc((size_t)(p + 1), sizeof(double));
    if (!a_prev || !a_curr) { free(a_prev); free(a_curr); return -1; }

    double E = r[0];

    for (int m = 1; m <= p; m++) {
        /* Compute reflection coefficient k[m] */
        double sum = r[m];
        for (int i = 1; i < m; i++)
            sum += a_prev[i] * r[m - i];
        double km = -sum / E;

        /* Check stability: |k| must be < 1 */
        if (fabs(km) >= 1.0) km = (km > 0) ? 0.9999 : -0.9999;

        /* Update coefficients */
        a_curr[m] = km;
        for (int i = 1; i < m; i++)
            a_curr[i] = a_prev[i] + km * a_prev[m - i];

        /* Update error energy */
        E *= (1.0 - km * km);

        if (k_out) k_out[m - 1] = km;

        /* Copy current to previous */
        memcpy(a_prev, a_curr, (size_t)(p + 1) * sizeof(double));
    }

    /* Output: a[0..p-1] = a_curr[1..p] */
    for (int i = 0; i < p; i++)
        a[i] = a_curr[i + 1];

    if (E_out) *E_out = E;

    free(a_prev);
    free(a_curr);
    return 0;
}

/* ================================================================== */
/*  Convenience wrappers                                               */
/* ================================================================== */

int lpc_coefficients(const double *x, int n, int p,
                     double *a, double *E_out)
{
    double *r = (double *)malloc((size_t)(p + 1) * sizeof(double));
    if (!r) return -1;

    lpc_autocorrelation(x, n, r, p);
    int ret = levinson_durbin(r, p, a, NULL, E_out);

    free(r);
    return ret;
}

void lpc_residual(const double *x, int n, const double *a, int p,
                  double *e)
{
    for (int i = 0; i < n; i++) {
        double pred = 0.0;
        for (int k = 1; k <= p && k <= i; k++)
            pred += a[k - 1] * x[i - k];
        e[i] = x[i] + pred;  /* a[k] already has sign convention */
    }
}

void lpc_synthesise(const double *e, int n, const double *a, int p,
                    double *x)
{
    for (int i = 0; i < n; i++) {
        double pred = 0.0;
        for (int k = 1; k <= p && k <= i; k++)
            pred += a[k - 1] * x[i - k];
        x[i] = e[i] - pred;
    }
}

void lpc_spectrum(const double *a, int p, double E,
                  double *spec, int nfft)
{
    int half = nfft / 2;
    for (int i = 0; i < half; i++) {
        double f = (double)i / (double)nfft;  /* normalised freq 0..0.5 */
        double w = 2.0 * M_PI * f;

        /* A(e^{jω}) = 1 + Σ a[k] · e^{-jωk} */
        double re = 1.0, im = 0.0;
        for (int k = 0; k < p; k++) {
            re += a[k] * cos(-(double)(k + 1) * w);
            im += a[k] * sin(-(double)(k + 1) * w);
        }

        /* |A(e^{jω})|² */
        double mag2 = re * re + im * im;

        /* S(f) = E / |A|² in dB */
        spec[i] = 10.0 * log10(E / (mag2 + 1e-30) + 1e-30);
    }
}

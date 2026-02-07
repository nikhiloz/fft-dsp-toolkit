/**
 * @file adaptive.c
 * @brief Adaptive filter implementations — LMS, NLMS, RLS.
 *
 * ── Algorithm Comparison ─────────────────────────────────────────
 *
 *   Algorithm   Cost/sample   Convergence     Tracking
 *   ─────────   ───────────   ───────────     ────────
 *   LMS         O(L)          Slow (~10·L)    Moderate
 *   NLMS        O(L)          Medium (~5·L)   Good
 *   RLS         O(L²)         Fast  (~2·L)    Excellent
 *
 *   L = filter length (taps)
 *
 * ── LMS Weight Update ────────────────────────────────────────────
 *
 *   ┌─────────────────────────────────────────┐
 *   │  for k = 0..L-1:                        │
 *   │    w[k] += μ · e[n] · x[n-k]           │
 *   │                                         │
 *   │  Gradient descent on MSE surface        │
 *   │  E{e²} is quadratic → bowl shape        │
 *   │  μ controls step size down the bowl     │
 *   └─────────────────────────────────────────┘
 *
 * ── RLS Matrix Recursion ─────────────────────────────────────────
 *
 *   P[n] tracks inverse of weighted autocorrelation:
 *
 *   P ← (1/λ)(P - k·xᵀ·P)
 *
 *   Woodbury identity avoids explicit matrix inversion.
 *   O(L²) storage + compute per sample.
 */

#include "adaptive.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ================================================================== */
/*  LMS — Least Mean Squares                                          */
/* ================================================================== */

int lms_init(LmsState *s, int taps, double mu)
{
    s->taps = taps;
    s->mu = mu;
    s->pos = 0;
    s->w = (double *)calloc((size_t)taps, sizeof(double));
    s->x_buf = (double *)calloc((size_t)taps, sizeof(double));
    if (!s->w || !s->x_buf) { lms_free(s); return -1; }
    return 0;
}

void lms_update(LmsState *s, double x, double d,
                double *y_out, double *e_out)
{
    int L = s->taps;

    /* Insert new sample into circular delay line */
    s->x_buf[s->pos] = x;

    /* Compute output: y = wᵀ · x_buf */
    double y = 0.0;
    for (int k = 0; k < L; k++) {
        int idx = (s->pos - k + L) % L;
        y += s->w[k] * s->x_buf[idx];
    }

    /* Error */
    double e = d - y;

    /* Weight update: w[k] += μ · e · x[n-k] */
    for (int k = 0; k < L; k++) {
        int idx = (s->pos - k + L) % L;
        s->w[k] += s->mu * e * s->x_buf[idx];
    }

    /* Advance circular buffer position */
    s->pos = (s->pos + 1) % L;

    *y_out = y;
    *e_out = e;
}

void lms_free(LmsState *s)
{
    free(s->w);     s->w = NULL;
    free(s->x_buf); s->x_buf = NULL;
}

/* ================================================================== */
/*  NLMS — Normalised LMS                                             */
/* ================================================================== */

int nlms_init(NlmsState *s, int taps, double mu, double eps)
{
    s->taps = taps;
    s->mu = mu;
    s->eps = eps;
    s->pos = 0;
    s->w = (double *)calloc((size_t)taps, sizeof(double));
    s->x_buf = (double *)calloc((size_t)taps, sizeof(double));
    if (!s->w || !s->x_buf) { nlms_free(s); return -1; }
    return 0;
}

void nlms_update(NlmsState *s, double x, double d,
                 double *y_out, double *e_out)
{
    int L = s->taps;
    s->x_buf[s->pos] = x;

    /* Compute output and input energy */
    double y = 0.0;
    double energy = 0.0;
    for (int k = 0; k < L; k++) {
        int idx = (s->pos - k + L) % L;
        double xk = s->x_buf[idx];
        y += s->w[k] * xk;
        energy += xk * xk;
    }

    double e = d - y;

    /* Normalised weight update */
    double step = s->mu / (energy + s->eps);
    for (int k = 0; k < L; k++) {
        int idx = (s->pos - k + L) % L;
        s->w[k] += step * e * s->x_buf[idx];
    }

    s->pos = (s->pos + 1) % L;
    *y_out = y;
    *e_out = e;
}

void nlms_free(NlmsState *s)
{
    free(s->w);     s->w = NULL;
    free(s->x_buf); s->x_buf = NULL;
}

/* ================================================================== */
/*  RLS — Recursive Least Squares                                     */
/* ================================================================== */

int rls_init(RlsState *s, int taps, double lambda, double delta)
{
    s->taps = taps;
    s->lambda = lambda;
    s->pos = 0;
    s->w = (double *)calloc((size_t)taps, sizeof(double));
    s->x_buf = (double *)calloc((size_t)taps, sizeof(double));
    s->k = (double *)calloc((size_t)taps, sizeof(double));
    s->P = (double *)calloc((size_t)(taps * taps), sizeof(double));
    if (!s->w || !s->x_buf || !s->k || !s->P) { rls_free(s); return -1; }

    /* Initialise P = (1/delta) · I */
    double inv_delta = 1.0 / delta;
    for (int i = 0; i < taps; i++)
        s->P[i * taps + i] = inv_delta;

    return 0;
}

void rls_update(RlsState *s, double x, double d,
                double *y_out, double *e_out)
{
    int L = s->taps;
    s->x_buf[s->pos] = x;

    /* Build input vector (from circular buffer) */
    double *xv = (double *)malloc((size_t)L * sizeof(double));
    for (int k = 0; k < L; k++)
        xv[k] = s->x_buf[(s->pos - k + L) % L];

    /* Compute P·x */
    double *Px = (double *)malloc((size_t)L * sizeof(double));
    for (int i = 0; i < L; i++) {
        Px[i] = 0.0;
        for (int j = 0; j < L; j++)
            Px[i] += s->P[i * L + j] * xv[j];
    }

    /* Compute denominator: λ + xᵀ·P·x */
    double denom = s->lambda;
    for (int i = 0; i < L; i++)
        denom += xv[i] * Px[i];

    /* Gain vector: k = P·x / denom */
    for (int i = 0; i < L; i++)
        s->k[i] = Px[i] / denom;

    /* A priori error: e = d - wᵀx */
    double y = 0.0;
    for (int i = 0; i < L; i++)
        y += s->w[i] * xv[i];
    double e = d - y;

    /* Weight update: w += k·e */
    for (int i = 0; i < L; i++)
        s->w[i] += s->k[i] * e;

    /* P update: P = (1/λ)(P - k·xᵀ·P) */
    /* P_new[i][j] = (1/λ)(P[i][j] - k[i] · (xᵀ·P)[j]) */
    /* Note: (xᵀ·P)[j] = Σ_m x[m]·P[m][j] = Px transposed if P symmetric */
    double inv_lam = 1.0 / s->lambda;
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            /* xᵀ·P column j = Σ_m x[m]·P[m][j] */
            double xPj = 0.0;
            for (int m = 0; m < L; m++)
                xPj += xv[m] * s->P[m * L + j];
            s->P[i * L + j] = inv_lam * (s->P[i * L + j] - s->k[i] * xPj);
        }
    }

    s->pos = (s->pos + 1) % L;
    *y_out = y;
    *e_out = e;

    free(xv);
    free(Px);
}

void rls_free(RlsState *s)
{
    free(s->w);     s->w = NULL;
    free(s->x_buf); s->x_buf = NULL;
    free(s->P);     s->P = NULL;
    free(s->k);     s->k = NULL;
}

/* ================================================================== */
/*  Batch convenience wrappers                                        */
/* ================================================================== */

void lms_filter(const double *x, const double *d, int n,
                int taps, double mu,
                double *y, double *e, double *w_final)
{
    LmsState s;
    lms_init(&s, taps, mu);
    for (int i = 0; i < n; i++)
        lms_update(&s, x[i], d[i], &y[i], &e[i]);
    if (w_final)
        memcpy(w_final, s.w, (size_t)taps * sizeof(double));
    lms_free(&s);
}

void nlms_filter(const double *x, const double *d, int n,
                 int taps, double mu, double eps,
                 double *y, double *e, double *w_final)
{
    NlmsState s;
    nlms_init(&s, taps, mu, eps);
    for (int i = 0; i < n; i++)
        nlms_update(&s, x[i], d[i], &y[i], &e[i]);
    if (w_final)
        memcpy(w_final, s.w, (size_t)taps * sizeof(double));
    nlms_free(&s);
}

void rls_filter(const double *x, const double *d, int n,
                int taps, double lambda, double delta,
                double *y, double *e, double *w_final)
{
    RlsState s;
    rls_init(&s, taps, lambda, delta);
    for (int i = 0; i < n; i++)
        rls_update(&s, x[i], d[i], &y[i], &e[i]);
    if (w_final)
        memcpy(w_final, s.w, (size_t)taps * sizeof(double));
    rls_free(&s);
}

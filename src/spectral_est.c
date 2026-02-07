/**
 * @file spectral_est.c
 * @brief Parametric spectral estimation — MUSIC, Capon, Jacobi eigensolver.
 *
 * ── Jacobi Eigendecomposition ────────────────────────────────────
 *
 *   For real symmetric A:
 *     1. Find largest off-diagonal |A[p][q]|
 *     2. Compute rotation angle θ to zero A[p][q]
 *     3. Apply Givens rotation: A ← GᵀAG
 *     4. Repeat until off-diagonal < tolerance
 *
 *   Convergence: quadratic, typically 5-10 sweeps for small matrices.
 *
 * ── MUSIC Pseudospectrum ─────────────────────────────────────────
 *
 *       1                  1
 *   ─────────── = ──────────────────────
 *   P_MUSIC(f)    Σ |eᴴ(f) · vₖ|²
 *                k∈noise
 *
 *   Sharp peaks at signal frequencies (not a true PSD).
 */

#include "spectral_est.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ================================================================== */
/*  Autocorrelation Matrix                                             */
/* ================================================================== */

void correlation_matrix(const double *x, int n, int p, double *R)
{
    /* R[i][j] = (1/(n-p+1)) * Σ_{m=0}^{n-p} x[m+i]*x[m+j] */
    int L = n - p + 1;
    if (L < 1) L = 1;
    for (int i = 0; i < p; i++) {
        for (int j = i; j < p; j++) {
            double sum = 0.0;
            for (int m = 0; m < L; m++)
                sum += x[m + i] * x[m + j];
            R[i * p + j] = sum / (double)L;
            R[j * p + i] = R[i * p + j];  /* symmetric */
        }
    }
}

/* ================================================================== */
/*  Jacobi Eigendecomposition (real symmetric)                         */
/* ================================================================== */

int eigen_symmetric(double *A, int p,
                    double *evals, double *evecs, int max_iter)
{
    /* Initialise V = Identity */
    for (int i = 0; i < p; i++)
        for (int j = 0; j < p; j++)
            evecs[i * p + j] = (i == j) ? 1.0 : 0.0;

    for (int iter = 0; iter < max_iter; iter++) {
        /* Find largest off-diagonal element */
        double max_off = 0.0;
        int ip = 0, iq = 1;
        for (int i = 0; i < p; i++) {
            for (int j = i + 1; j < p; j++) {
                if (fabs(A[i * p + j]) > max_off) {
                    max_off = fabs(A[i * p + j]);
                    ip = i; iq = j;
                }
            }
        }

        if (max_off < 1e-12) break;  /* converged */

        /* Compute rotation to zero A[ip][iq] */
        double aii = A[ip * p + ip];
        double ajj = A[iq * p + iq];
        double aij = A[ip * p + iq];

        double theta;
        if (fabs(aii - ajj) < 1e-30)
            theta = M_PI / 4.0;
        else
            theta = 0.5 * atan2(2.0 * aij, aii - ajj);

        double c = cos(theta);
        double s = sin(theta);

        /* Apply Givens rotation: A ← GᵀAG */
        /* Update rows/cols ip and iq */
        double *tmp_i = (double *)malloc((size_t)p * sizeof(double));
        double *tmp_j = (double *)malloc((size_t)p * sizeof(double));

        for (int k = 0; k < p; k++) {
            tmp_i[k] = c * A[ip * p + k] + s * A[iq * p + k];
            tmp_j[k] = -s * A[ip * p + k] + c * A[iq * p + k];
        }
        for (int k = 0; k < p; k++) {
            A[ip * p + k] = tmp_i[k];
            A[iq * p + k] = tmp_j[k];
        }

        for (int k = 0; k < p; k++) {
            tmp_i[k] = c * A[k * p + ip] + s * A[k * p + iq];
            tmp_j[k] = -s * A[k * p + ip] + c * A[k * p + iq];
        }
        for (int k = 0; k < p; k++) {
            A[k * p + ip] = tmp_i[k];
            A[k * p + iq] = tmp_j[k];
        }

        /* Update eigenvector matrix: V ← V·G */
        for (int k = 0; k < p; k++) {
            double vi = evecs[k * p + ip];
            double vj = evecs[k * p + iq];
            evecs[k * p + ip] = c * vi + s * vj;
            evecs[k * p + iq] = -s * vi + c * vj;
        }

        free(tmp_i);
        free(tmp_j);
    }

    /* Extract eigenvalues from diagonal */
    for (int i = 0; i < p; i++)
        evals[i] = A[i * p + i];

    /* Sort by descending eigenvalue (bubble sort, p is small) */
    for (int i = 0; i < p - 1; i++) {
        for (int j = i + 1; j < p; j++) {
            if (evals[j] > evals[i]) {
                /* Swap eigenvalues */
                double tmp = evals[i];
                evals[i] = evals[j];
                evals[j] = tmp;
                /* Swap corresponding eigenvector columns */
                for (int k = 0; k < p; k++) {
                    tmp = evecs[k * p + i];
                    evecs[k * p + i] = evecs[k * p + j];
                    evecs[k * p + j] = tmp;
                }
            }
        }
    }

    return 0;
}

/* ================================================================== */
/*  MUSIC Pseudospectrum                                               */
/* ================================================================== */

void music_spectrum(const double *x, int n, int p, int n_sigs,
                    double *spec, int nfft)
{
    int half = nfft / 2;

    /* Build correlation matrix */
    double *R = (double *)malloc((size_t)(p * p) * sizeof(double));
    correlation_matrix(x, n, p, R);

    /* Eigendecompose */
    double *evals = (double *)malloc((size_t)p * sizeof(double));
    double *evecs = (double *)malloc((size_t)(p * p) * sizeof(double));
    eigen_symmetric(R, p, evals, evecs, 200);

    /* Noise subspace: eigenvectors n_sigs..p-1 (smallest eigenvalues) */
    int n_noise = p - n_sigs;

    for (int fi = 0; fi < half; fi++) {
        double f = (double)fi / (double)nfft;
        double w = 2.0 * M_PI * f;

        /* Compute eᴴ(f) · E_n · E_nᴴ · e(f) = Σ_k |eᴴ·v_k|² */
        double denom = 0.0;
        for (int k = 0; k < n_noise; k++) {
            int kk = n_sigs + k;  /* noise eigenvector index */
            double re_dot = 0.0, im_dot = 0.0;
            for (int m = 0; m < p; m++) {
                double e_re = cos(w * (double)m);
                double e_im = -sin(w * (double)m);
                double v = evecs[m * p + kk];
                re_dot += e_re * v;
                im_dot += e_im * v;
            }
            denom += re_dot * re_dot + im_dot * im_dot;
        }

        spec[fi] = 10.0 * log10(1.0 / (denom + 1e-30));
    }

    free(R); free(evals); free(evecs);
}

void music_frequencies(const double *x, int n, int p, int n_sigs,
                       double *freqs)
{
    /* Compute MUSIC spectrum at high resolution */
    int nfft = 4096;
    int half = nfft / 2;
    double *spec = (double *)malloc((size_t)half * sizeof(double));
    music_spectrum(x, n, p, n_sigs, spec, nfft);

    /* Find top n_sigs peaks */
    for (int s = 0; s < n_sigs; s++) freqs[s] = 0.0;

    for (int s = 0; s < n_sigs; s++) {
        double max_val = -1e30;
        int max_idx = 0;
        for (int i = 1; i < half - 1; i++) {
            /* Is local peak? */
            if (spec[i] > spec[i - 1] && spec[i] > spec[i + 1] &&
                spec[i] > max_val) {
                /* Check not too close to already found peaks */
                int too_close = 0;
                for (int ss = 0; ss < s; ss++) {
                    int prev_bin = (int)(freqs[ss] * (double)nfft);
                    if (abs(i - prev_bin) < nfft / 100) {
                        too_close = 1;
                        break;
                    }
                }
                if (!too_close) {
                    max_val = spec[i];
                    max_idx = i;
                }
            }
        }
        freqs[s] = (double)max_idx / (double)nfft;
    }

    /* Sort frequencies ascending */
    for (int i = 0; i < n_sigs - 1; i++)
        for (int j = i + 1; j < n_sigs; j++)
            if (freqs[j] < freqs[i]) {
                double tmp = freqs[i];
                freqs[i] = freqs[j];
                freqs[j] = tmp;
            }

    free(spec);
}

/* ================================================================== */
/*  Capon (Minimum Variance) Spectral Estimator                        */
/* ================================================================== */

/*
 * Invert p×p symmetric positive-definite matrix using Cholesky decomp.
 * A_inv = (LLᵀ)⁻¹. Input A is destroyed. Result in A_inv.
 */
static int invert_spd(const double *R, int p, double *R_inv)
{
    /* Cholesky: R = LLᵀ */
    double *L = (double *)calloc((size_t)(p * p), sizeof(double));
    if (!L) return -1;

    for (int i = 0; i < p; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = R[i * p + j];
            for (int k = 0; k < j; k++)
                sum -= L[i * p + k] * L[j * p + k];
            if (i == j) {
                if (sum <= 0.0) { free(L); return -1; }
                L[i * p + j] = sqrt(sum);
            } else {
                L[i * p + j] = sum / L[j * p + j];
            }
        }
    }

    /* Invert L (lower triangular) → L_inv */
    double *Li = (double *)calloc((size_t)(p * p), sizeof(double));
    if (!Li) { free(L); return -1; }

    for (int i = 0; i < p; i++) {
        Li[i * p + i] = 1.0 / L[i * p + i];
        for (int j = i + 1; j < p; j++) {
            double sum = 0.0;
            for (int k = i; k < j; k++)
                sum -= L[j * p + k] * Li[k * p + i];
            Li[j * p + i] = sum / L[j * p + j];
        }
    }

    /* R_inv = L_inv^T · L_inv */
    for (int i = 0; i < p; i++)
        for (int j = 0; j < p; j++) {
            double sum = 0.0;
            for (int k = 0; k < p; k++)
                sum += Li[k * p + i] * Li[k * p + j];
            R_inv[i * p + j] = sum;
        }

    free(L); free(Li);
    return 0;
}

void capon_spectrum(const double *x, int n, int p,
                    double *spec, int nfft)
{
    int half = nfft / 2;

    double *R = (double *)malloc((size_t)(p * p) * sizeof(double));
    double *R_inv = (double *)malloc((size_t)(p * p) * sizeof(double));
    correlation_matrix(x, n, p, R);

    if (invert_spd(R, p, R_inv) != 0) {
        /* Fallback: return flat spectrum */
        for (int i = 0; i < half; i++) spec[i] = 0.0;
        free(R); free(R_inv);
        return;
    }

    for (int fi = 0; fi < half; fi++) {
        double f = (double)fi / (double)nfft;
        double w = 2.0 * M_PI * f;

        /* eᴴ · R⁻¹ · e */
        /* First compute R⁻¹ · e */
        double *Rie = (double *)malloc((size_t)p * sizeof(double));
        for (int i = 0; i < p; i++) {
            double re_sum = 0.0;
            for (int j = 0; j < p; j++) {
                double e_re = cos(w * (double)j);
                double e_im = sin(w * (double)j);  /* conjugate of e^{-jωj} */
                /* R_inv is real, so R_inv·e has both real and imag parts */
                /* For real R_inv: (R⁻¹·e)_i = Σ_j R⁻¹[i][j]·e^{-jωj} */
                re_sum += R_inv[i * p + j] * e_re;
                (void)e_im;
            }
            Rie[i] = re_sum;
        }

        /* For real symmetric R_inv, eᴴ·R⁻¹·e = Σ_{i,j} R_inv[i][j]·cos(ω(i-j)) */
        double denom = 0.0;
        for (int i = 0; i < p; i++)
            for (int j = 0; j < p; j++)
                denom += R_inv[i * p + j] * cos(w * (double)(i - j));

        spec[fi] = 10.0 * log10(1.0 / (fabs(denom) + 1e-30));

        free(Rie);
    }

    free(R); free(R_inv);
}

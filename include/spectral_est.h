/**
 * @file spectral_est.h
 * @brief Parametric spectral estimation — MUSIC, Pisarenko, min-variance.
 *
 * ── MUSIC Algorithm ──────────────────────────────────────────────
 *
 *   1. Form data autocorrelation matrix R_xx (p×p).
 *   2. Eigendecompose R_xx → signal subspace + noise subspace.
 *   3. P_MUSIC(f) = 1 / (eᴴ(f) · E_n · E_nᴴ · e(f))
 *
 *   where e(f) = [1, e^{-j2πf}, ..., e^{-j2πf(p-1)}]ᵀ
 *         E_n = noise eigenvectors (those with smallest eigenvalues)
 *
 *   ┌─────────────────────────────────────────────────┐
 *   │   R_xx eigenvalues:                              │
 *   │   λ₁ ≥ λ₂ ≥ ... ≥ λ_M > λ_{M+1} ≈ ... ≈ λ_p  │
 *   │   ╰──── signal ────╯    ╰──── noise ────────╯   │
 *   │                                                  │
 *   │   M signal eigenvalues >> σ² (noise floor)       │
 *   └─────────────────────────────────────────────────┘
 *
 * ── Pisarenko Harmonic Decomposition ─────────────────────────────
 *
 *   Special case of MUSIC with noise subspace dimension = 1.
 *   For M sinusoids, use correlation matrix of size (M+1)×(M+1).
 *   Minimum eigenvalue eigenvector gives frequency estimates.
 *
 * ── Minimum Variance (Capon) ─────────────────────────────────────
 *
 *   P_MV(f) = 1 / (eᴴ(f) · R_xx⁻¹ · e(f))
 *
 *   Better resolution than periodogram, but requires matrix inverse.
 */

#ifndef SPECTRAL_EST_H
#define SPECTRAL_EST_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Estimate autocorrelation matrix from signal data.
 *
 * R[i][j] = (1/N) Σ x[n+i]·x[n+j], stored row-major.
 *
 * @param x      Input signal (length n).
 * @param n      Signal length.
 * @param p      Matrix size (p × p).
 * @param R      Output matrix (length p*p, row-major).
 */
void correlation_matrix(const double *x, int n, int p, double *R);

/**
 * @brief Eigendecompose a real symmetric matrix using Jacobi iteration.
 *
 * @param A      Input matrix (p×p, row-major) — DESTROYED on output.
 * @param p      Matrix size.
 * @param evals  Output eigenvalues (length p), sorted descending.
 * @param evecs  Output eigenvectors (p×p, column-major: evec[i] = evecs[i*p..]).
 * @param max_iter Maximum Jacobi sweeps (typically 100).
 * @return       0 on success, -1 on convergence failure.
 */
int eigen_symmetric(double *A, int p,
                    double *evals, double *evecs, int max_iter);

/**
 * @brief MUSIC pseudospectrum estimation.
 *
 * @param x       Input signal (length n).
 * @param n       Signal length.
 * @param p       Autocorrelation matrix dimension.
 * @param n_sigs  Number of signal sources (sinusoids).
 * @param spec    Output MUSIC pseudospectrum in dB (length nfft/2).
 * @param nfft    Number of frequency evaluation points.
 */
void music_spectrum(const double *x, int n, int p, int n_sigs,
                    double *spec, int nfft);

/**
 * @brief MUSIC frequency estimation — return peak frequencies.
 *
 * @param x       Input signal (length n).
 * @param n       Signal length.
 * @param p       Matrix dimension (must be > n_sigs).
 * @param n_sigs  Number of sinusoidal components.
 * @param freqs   Output normalised frequencies (length n_sigs, range 0–0.5).
 */
void music_frequencies(const double *x, int n, int p, int n_sigs,
                       double *freqs);

/**
 * @brief Capon (Minimum Variance) spectral estimator.
 *
 * @param x       Input signal (length n).
 * @param n       Signal length.
 * @param p       Matrix dimension.
 * @param spec    Output spectrum in dB (length nfft/2).
 * @param nfft    Number of frequency evaluation points.
 */
void capon_spectrum(const double *x, int n, int p,
                    double *spec, int nfft);

#ifdef __cplusplus
}
#endif

#endif /* SPECTRAL_EST_H */

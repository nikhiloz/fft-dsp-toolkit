/**
 * @file lpc.h
 * @brief Linear Prediction Coding — AR modelling, Levinson-Durbin recursion.
 *
 * ── Linear Prediction Model ──────────────────────────────────────
 *
 *   x̂[n] = - Σ_{k=1}^{p} a[k] · x[n-k]
 *
 *   Prediction error: e[n] = x[n] - x̂[n]
 *   Minimise E{e²} → Wiener-Hopf (normal) equations.
 *
 * ── Autocorrelation Method ───────────────────────────────────────
 *
 *      ┌ r[0]   r[1]   ... r[p-1] ┐   ┌ a[1] ┐     ┌ r[1] ┐
 *      │ r[1]   r[0]   ... r[p-2] │   │ a[2] │     │ r[2] │
 *      │  :      :    ·     :     │ · │  :   │ = - │  :   │
 *      └ r[p-1] r[p-2] ... r[0]  ┘   └ a[p] ┘     └ r[p] ┘
 *
 *   Toeplitz structure → Levinson-Durbin solves in O(p²).
 *
 * ── Levinson-Durbin Recursion ────────────────────────────────────
 *
 *   For each order m = 1..p:
 *     1. k[m] = -(r[m] + Σ a_{m-1}[i]·r[m-i]) / E[m-1]
 *     2. a_m[i] = a_{m-1}[i] + k[m]·a_{m-1}[m-i]
 *     3. a_m[m] = k[m]
 *     4. E[m] = (1 - k[m]²)·E[m-1]
 *
 *   k[m] = reflection coefficients (PARCOR), |k[m]| < 1 ↔ stable.
 */

#ifndef LPC_H
#define LPC_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Compute autocorrelation r[0..p] of signal x.
 *
 * @param x      Input signal (length n).
 * @param n      Signal length.
 * @param r      Output autocorrelation (length p+1).
 * @param p      Maximum lag (prediction order).
 */
void lpc_autocorrelation(const double *x, int n, double *r, int p);

/**
 * @brief Levinson-Durbin recursion — solve LPC normal equations.
 *
 * Given autocorrelation r[0..p], computes LP coefficients a[1..p],
 * reflection coefficients k[1..p], and prediction error energy E.
 *
 * @param r      Autocorrelation values (length p+1).
 * @param p      Prediction order.
 * @param a      Output LP coefficients a[1..p] (array of length p).
 * @param k_out  Output reflection coefficients (length p, may be NULL).
 * @param E_out  Output prediction error energy (may be NULL).
 * @return       0 on success, -1 if r[0] == 0 or unstable.
 */
int levinson_durbin(const double *r, int p,
                    double *a, double *k_out, double *E_out);

/**
 * @brief Compute LPC coefficients from a signal (autocorrelation method).
 *
 * Convenience wrapper: computes autocorrelation then Levinson-Durbin.
 *
 * @param x      Input signal (length n).
 * @param n      Signal length.
 * @param p      Prediction order.
 * @param a      Output LP coefficients a[1..p] (length p).
 * @param E_out  Output prediction error power (may be NULL).
 * @return       0 on success, -1 on failure.
 */
int lpc_coefficients(const double *x, int n, int p,
                     double *a, double *E_out);

/**
 * @brief Compute LP prediction error (residual) signal.
 *
 * e[n] = x[n] + Σ_{k=1}^{p} a[k]·x[n-k]
 *
 * @param x      Input signal (length n).
 * @param n      Signal length.
 * @param a      LP coefficients a[1..p] (length p).
 * @param p      Prediction order.
 * @param e      Output residual (length n).
 */
void lpc_residual(const double *x, int n, const double *a, int p,
                  double *e);

/**
 * @brief Synthesise signal from LP coefficients and excitation.
 *
 * x[n] = e[n] - Σ_{k=1}^{p} a[k]·x[n-k]
 *
 * @param e      Excitation (residual) signal (length n).
 * @param n      Signal length.
 * @param a      LP coefficients (length p).
 * @param p      Prediction order.
 * @param x      Output synthesised signal (length n).
 */
void lpc_synthesise(const double *e, int n, const double *a, int p,
                    double *x);

/**
 * @brief Compute LP spectral envelope (AR power spectrum).
 *
 * S(f) = E / |A(e^{j2πf})|²
 *
 * where A(z) = 1 + a[1]z⁻¹ + ... + a[p]z⁻ᵖ.
 *
 * @param a       LP coefficients a[1..p] (length p).
 * @param p       Prediction order.
 * @param E       Prediction error energy.
 * @param spec    Output spectral envelope in dB (length nfft/2).
 * @param nfft    FFT size (power of 2).
 */
void lpc_spectrum(const double *a, int p, double E,
                  double *spec, int nfft);

#ifdef __cplusplus
}
#endif

#endif /* LPC_H */

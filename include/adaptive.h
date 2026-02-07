/**
 * @file adaptive.h
 * @brief Adaptive filters — LMS, NLMS, RLS algorithms.
 *
 * ── LMS (Least Mean Squares) ─────────────────────────────────────
 *
 *   x[n] ──►┌──────────┐──► y[n] ──►(+)──► e[n]
 *            │  w[k]    │          -│
 *            │  FIR     │     d[n] ─┘
 *            └──────────┘
 *                 ▲
 *    w[n+1] = w[n] + μ·e[n]·x[n]
 *
 *   - μ (step size): controls convergence speed vs stability
 *   - Stable if 0 < μ < 2/(L·P_x), where L=taps, P_x=input power
 *   - Converges to Wiener solution: w_opt = R_xx⁻¹ · r_xd
 *
 * ── NLMS (Normalised LMS) ───────────────────────────────────────
 *
 *    w[n+1] = w[n] + (μ / (‖x‖² + ε)) · e[n] · x[n]
 *
 *   - Step size normalised by input power → robust to signal level
 *   - ε prevents division by zero (regularisation)
 *
 * ── RLS (Recursive Least Squares) ───────────────────────────────
 *
 *    k[n] = P[n-1]·x[n] / (λ + x[n]ᵀ·P[n-1]·x[n])
 *    e[n] = d[n] - w[n-1]ᵀ·x[n]
 *    w[n] = w[n-1] + k[n]·e[n]
 *    P[n] = (1/λ)·(P[n-1] - k[n]·x[n]ᵀ·P[n-1])
 *
 *   - λ (forgetting factor): 0.95–1.0; lower = faster tracking
 *   - Converges in ~L iterations (much faster than LMS)
 *   - Cost: O(L²) per sample vs O(L) for LMS
 */

#ifndef ADAPTIVE_H
#define ADAPTIVE_H

#ifdef __cplusplus
extern "C" {
#endif

/* ── LMS State ──────────────────────────────────────────────────── */

typedef struct {
    double *w;          /**< Filter coefficients (length taps) */
    double *x_buf;      /**< Input delay line (length taps) */
    int taps;           /**< Number of filter taps */
    double mu;          /**< Step size */
    int pos;            /**< Current position in delay line */
} LmsState;

/**
 * @brief Initialise LMS adaptive filter.
 * @param s      State structure to initialise.
 * @param taps   Number of filter taps.
 * @param mu     Step size (0 < mu < 2/(taps*P_x)).
 * @return       0 on success, -1 on allocation failure.
 */
int lms_init(LmsState *s, int taps, double mu);

/**
 * @brief Process one sample through LMS adaptive filter.
 * @param s      Initialised LMS state.
 * @param x      Input sample.
 * @param d      Desired (reference) sample.
 * @param y_out  Output: filter output y[n] = wᵀx.
 * @param e_out  Output: error e[n] = d[n] - y[n].
 */
void lms_update(LmsState *s, double x, double d,
                double *y_out, double *e_out);

/**
 * @brief Free LMS state resources.
 */
void lms_free(LmsState *s);

/* ── NLMS State ─────────────────────────────────────────────────── */

typedef struct {
    double *w;
    double *x_buf;
    int taps;
    double mu;          /**< Normalised step size (0 < mu < 2) */
    double eps;         /**< Regularisation constant */
    int pos;
} NlmsState;

/**
 * @brief Initialise NLMS adaptive filter.
 * @param s      State structure.
 * @param taps   Number of taps.
 * @param mu     Normalised step size (typically 0.1–1.0).
 * @param eps    Regularisation (small positive, e.g. 1e-8).
 * @return       0 on success, -1 on failure.
 */
int nlms_init(NlmsState *s, int taps, double mu, double eps);

/**
 * @brief Process one sample through NLMS.
 */
void nlms_update(NlmsState *s, double x, double d,
                 double *y_out, double *e_out);

/**
 * @brief Free NLMS state resources.
 */
void nlms_free(NlmsState *s);

/* ── RLS State ──────────────────────────────────────────────────── */

typedef struct {
    double *w;          /**< Filter coefficients (length taps) */
    double *x_buf;      /**< Input delay line (length taps) */
    double *P;          /**< Inverse correlation matrix (taps × taps) */
    double *k;          /**< Gain vector (length taps) */
    int taps;
    double lambda;      /**< Forgetting factor (0.95–1.0) */
    int pos;
} RlsState;

/**
 * @brief Initialise RLS adaptive filter.
 * @param s       State structure.
 * @param taps    Number of taps.
 * @param lambda  Forgetting factor (typically 0.99–1.0).
 * @param delta   Initial P = (1/delta)*I (typically 100–1000).
 * @return        0 on success, -1 on failure.
 */
int rls_init(RlsState *s, int taps, double lambda, double delta);

/**
 * @brief Process one sample through RLS.
 */
void rls_update(RlsState *s, double x, double d,
                double *y_out, double *e_out);

/**
 * @brief Free RLS state resources.
 */
void rls_free(RlsState *s);

/* ── Convenience Functions ──────────────────────────────────────── */

/**
 * @brief Run LMS on an entire signal (batch mode).
 *
 * @param x      Input signal (length n).
 * @param d      Desired signal (length n).
 * @param n      Signal length.
 * @param taps   Number of filter taps.
 * @param mu     Step size.
 * @param y      Output: filter outputs (length n).
 * @param e      Output: error signal (length n).
 * @param w_final Output: final filter coefficients (length taps, may be NULL).
 */
void lms_filter(const double *x, const double *d, int n,
                int taps, double mu,
                double *y, double *e, double *w_final);

/**
 * @brief Run NLMS on an entire signal (batch mode).
 */
void nlms_filter(const double *x, const double *d, int n,
                 int taps, double mu, double eps,
                 double *y, double *e, double *w_final);

/**
 * @brief Run RLS on an entire signal (batch mode).
 */
void rls_filter(const double *x, const double *d, int n,
                int taps, double lambda, double delta,
                double *y, double *e, double *w_final);

#ifdef __cplusplus
}
#endif

#endif /* ADAPTIVE_H */

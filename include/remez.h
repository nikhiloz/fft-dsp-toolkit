/**
 * @file remez.h
 * @brief Parks-McClellan (Remez exchange) optimal FIR filter design.
 *
 * ── Remez Exchange Algorithm ─────────────────────────────────────
 *
 *   Given: band edges, desired response, weights
 *
 *   1. Initialise extremal frequencies (Chebyshev spacing)
 *   2. LOOP:
 *      a. Solve for optimal polynomial via Lagrange interpolation
 *      b. Compute error E(ω) = W(ω)·[D(ω) - H(ω)]
 *      c. Find new extremal set (peaks of |E|)
 *      d. If converged (δ stable), done; else goto 2
 *   3. Extract filter coefficients from optimal polynomial
 *
 * ── Equiripple Property ──────────────────────────────────────────
 *
 *   |E(ω)|        Alternation theorem:
 *    δ  ╌╌╌╌┬╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌   num extrema = L + 2
 *            │  ╱╲    ╱╲    ╱╲    (L = filter half-length)
 *        ────┼─╱──╲──╱──╲──╱──╲─
 *            │╱    ╲╱    ╲╱    ╲
 *   -δ  ╌╌╌╌┴╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
 *
 *   Error oscillates with equal amplitude — optimal minimax design.
 */

#ifndef REMEZ_H
#define REMEZ_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Band specification for Remez filter design.
 */
typedef struct {
    double low;       /**< Band lower edge (normalised, 0–0.5) */
    double high;      /**< Band upper edge (normalised, 0–0.5) */
    double desired;   /**< Desired gain in band (0.0 or 1.0 typically) */
    double weight;    /**< Relative weight (higher = stricter) */
} RemezBand;

/**
 * @brief Design an optimal equiripple FIR filter using the Remez exchange algorithm.
 *
 * @param h         Output filter coefficients (length taps).
 * @param taps      Number of taps (odd for Type I, even for Type II).
 * @param bands     Array of band specifications.
 * @param n_bands   Number of bands (typically 2 for LP/HP, 3 for BP/BS).
 * @param max_iter  Maximum Remez iterations (typically 25-50).
 * @return          0 on success, -1 on convergence failure.
 */
int remez_fir(double *h, int taps, const RemezBand *bands, int n_bands,
              int max_iter);

/**
 * @brief Design a lowpass FIR filter using Remez (convenience wrapper).
 *
 * @param h         Output coefficients (length taps).
 * @param taps      Number of taps.
 * @param fpass     Passband edge (normalised, 0–0.5).
 * @param fstop     Stopband edge (normalised, 0–0.5).
 * @param wpass     Passband weight (1.0 = standard).
 * @param wstop     Stopband weight (1.0 = standard).
 * @return          0 on success.
 */
int remez_lowpass(double *h, int taps, double fpass, double fstop,
                  double wpass, double wstop);

/**
 * @brief Design a bandpass FIR filter using Remez.
 *
 * @param h         Output coefficients (length taps).
 * @param taps      Number of taps.
 * @param fstop1    Lower stopband edge.
 * @param fpass1    Lower passband edge.
 * @param fpass2    Upper passband edge.
 * @param fstop2    Upper stopband edge.
 * @return          0 on success.
 */
int remez_bandpass(double *h, int taps,
                   double fstop1, double fpass1,
                   double fpass2, double fstop2);

#ifdef __cplusplus
}
#endif

#endif /* REMEZ_H */

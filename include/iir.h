/**
 * @file iir.h
 * @brief IIR (Infinite Impulse Response) digital filters — design, structures, processing.
 *
 * Provides biquad (second-order section) processing in multiple forms,
 * cascaded SOS chains, and classic analog-prototype IIR filter design
 * using the bilinear transform.
 *
 *   ┌─────────────────────────────────────────────────────────┐
 *   │              IIR Filter Signal Flow                     │
 *   │                                                         │
 *   │  x[n] ──►[ b0 ]──┬──► Σ ──────────────────► y[n]      │
 *   │           [ b1 ]◄─┘      ▲       ▲                     │
 *   │           [ b2 ]  z⁻¹   │       │                      │
 *   │                   ──►  Feedforward  Feedback            │
 *   │           [-a1 ]──────────┘       │                     │
 *   │           [-a2 ]  z⁻¹ ───────────┘                     │
 *   │                                                         │
 *   │  H(z) = (b0 + b1·z⁻¹ + b2·z⁻²)                       │
 *   │       / ( 1 + a1·z⁻¹ + a2·z⁻²)                        │
 *   └─────────────────────────────────────────────────────────┘
 *
 * TUTORIAL CROSS-REFERENCES:
 *   Frequency response    → chapters/06-frequency-response.md
 *   IIR filter design     → chapters/11-iir-filter-design.md
 *   Filter structures     → chapters/12-filter-structures.md
 *   FIR filters           → include/filter.h
 *
 * DESIGN NOTE:
 *   All design functions use normalized frequency where
 *   cutoff ∈ (0.0, 0.5):  0.0 = DC,  0.5 = Nyquist (fs/2).
 *   This matches the convention in filter.h (fir_lowpass).
 */

#ifndef IIR_H
#define IIR_H

#include "dsp_utils.h"

/* ── Biquad coefficient structure ────────────────────────────────── */

/**
 * @brief Second-Order Section (biquad) coefficients.
 *
 * Transfer function (a0 = 1 by normalisation):
 *
 *   H(z) = (b0 + b1·z⁻¹ + b2·z⁻²) / (1 + a1·z⁻¹ + a2·z⁻²)
 */
typedef struct {
    double b0, b1, b2;  /**< Numerator (feed-forward) coefficients  */
    double a1, a2;       /**< Denominator (feedback) coefficients    */
} Biquad;

/* ── Biquad state structures ─────────────────────────────────────── */

/**
 * @brief State for Direct Form I biquad processing.
 *
 * Stores two samples of input history and two samples of output history.
 *
 *   y[n] = b0·x[n] + b1·x[n-1] + b2·x[n-2]
 *                   - a1·y[n-1] - a2·y[n-2]
 */
typedef struct {
    double x1, x2;  /**< Previous input samples:  x[n-1], x[n-2] */
    double y1, y2;  /**< Previous output samples: y[n-1], y[n-2] */
} BiquadDF1State;

/**
 * @brief State for Direct Form II Transposed biquad processing.
 *
 * Only two state variables (vs four in DF1) — more memory efficient.
 *
 *   y[n] = b0·x[n] + s1
 *   s1   = b1·x[n] - a1·y[n] + s2
 *   s2   = b2·x[n] - a2·y[n]
 */
typedef struct {
    double s1, s2;  /**< Internal state variables */
} BiquadDF2TState;

/* ── Cascaded SOS (Second-Order Sections) ────────────────────────── */

/** Maximum number of biquad sections in a cascade (supports up to order 16) */
#define MAX_SOS_SECTIONS 8

/**
 * @brief Cascaded biquad filter — the recommended IIR structure.
 *
 * Factoring a high-order IIR into cascaded 2nd-order sections avoids
 * the coefficient-sensitivity and overflow problems of a single
 * high-order Direct Form implementation.
 *
 *   x[n] → [SOS 0] → [SOS 1] → ... → [SOS K-1] → × gain → y[n]
 */
typedef struct {
    Biquad sections[MAX_SOS_SECTIONS];      /**< Biquad coefficients        */
    BiquadDF1State states[MAX_SOS_SECTIONS]; /**< Per-section filter state   */
    int n_sections;                          /**< Number of active sections  */
    double gain;                             /**< Overall cascade gain       */
} SOSCascade;

/* ══════════════════════════════════════════════════════════════════
 *  General IIR filtering
 * ══════════════════════════════════════════════════════════════════ */

/**
 * @brief Apply a general IIR filter to a signal (direct form, zero initial state).
 *
 * Implements the difference equation:
 *   y[n] = (1/a[0]) * ( b[0]x[n] + b[1]x[n-1] + ...
 *                       − a[1]y[n-1] − a[2]y[n-2] − ... )
 *
 * @param b      Numerator coefficients,   length b_len
 * @param b_len  Number of numerator coefficients
 * @param a      Denominator coefficients, length a_len (a[0] may be != 1)
 * @param a_len  Number of denominator coefficients
 * @param in     Input signal,  length n
 * @param out    Output signal, length n (caller allocates)
 * @param n      Number of samples
 */
void iir_filter(const double *b, int b_len,
                const double *a, int a_len,
                const double *in, double *out, int n);

/**
 * @brief Compute the impulse response of an IIR system.
 *
 * Feeds δ[n] through the IIR filter and records the first @p n_out samples.
 *
 * @param b      Numerator coefficients
 * @param b_len  Number of numerator coefficients
 * @param a      Denominator coefficients
 * @param a_len  Number of denominator coefficients
 * @param h      Output impulse response buffer (caller allocates, length n_out)
 * @param n_out  Number of impulse response samples to compute
 */
void iir_impulse_response(const double *b, int b_len,
                          const double *a, int a_len,
                          double *h, int n_out);

/* ══════════════════════════════════════════════════════════════════
 *  Biquad per-sample processing
 * ══════════════════════════════════════════════════════════════════ */

/**
 * @brief Initialise Direct Form I biquad state to zero.
 * @param s  State to clear
 */
void biquad_df1_init(BiquadDF1State *s);

/**
 * @brief Initialise Direct Form II Transposed biquad state to zero.
 * @param s  State to clear
 */
void biquad_df2t_init(BiquadDF2TState *s);

/**
 * @brief Process one sample through a biquad using Direct Form I.
 *
 *   y[n] = b0·x[n] + b1·x[n-1] + b2·x[n-2]
 *                   - a1·y[n-1] - a2·y[n-2]
 *
 * @param bq  Biquad coefficients
 * @param s   Persistent filter state (updated on each call)
 * @param x   Input sample
 * @return    Output sample y[n]
 */
double biquad_process_df1(const Biquad *bq, BiquadDF1State *s, double x);

/**
 * @brief Process one sample through a biquad using Direct Form II Transposed.
 *
 * More numerically stable than DF1 for fixed-point or high-Q filters.
 *
 *   y[n] = b0·x[n] + s1
 *   s1   = b1·x[n] - a1·y[n] + s2
 *   s2   = b2·x[n] - a2·y[n]
 *
 * @param bq  Biquad coefficients
 * @param s   Persistent state (updated on each call)
 * @param x   Input sample
 * @return    Output sample y[n]
 */
double biquad_process_df2t(const Biquad *bq, BiquadDF2TState *s, double x);

/**
 * @brief Process a block of samples through a biquad (Direct Form I).
 *
 * @param bq   Biquad coefficients
 * @param s    Persistent state (updated across calls)
 * @param in   Input buffer,  length n
 * @param out  Output buffer, length n (caller allocates)
 * @param n    Number of samples
 */
void biquad_process_block(const Biquad *bq, BiquadDF1State *s,
                          const double *in, double *out, int n);

/* ══════════════════════════════════════════════════════════════════
 *  Cascaded SOS processing
 * ══════════════════════════════════════════════════════════════════ */

/**
 * @brief Initialise all states in a cascade to zero.
 * @param sos  Cascade to reset
 */
void sos_init(SOSCascade *sos);

/**
 * @brief Process one sample through the full SOS cascade.
 *
 * Signal flows through each biquad section in sequence, then is
 * multiplied by the overall gain.
 *
 * @param sos  Cascade (states updated in-place)
 * @param x    Input sample
 * @return     Output sample
 */
double sos_process_sample(SOSCascade *sos, double x);

/**
 * @brief Process a block of samples through the SOS cascade.
 *
 * @param sos  Cascade (states updated in-place)
 * @param in   Input buffer, length n
 * @param out  Output buffer, length n (caller allocates)
 * @param n    Number of samples
 */
void sos_process_block(SOSCascade *sos,
                       const double *in, double *out, int n);

/* ══════════════════════════════════════════════════════════════════
 *  IIR Filter Design — Butterworth
 *
 *  Butterworth filters are "maximally flat" in the passband:
 *    |H(jΩ)|² = 1 / (1 + (Ω/Ωc)^{2N})
 *
 *  Design method: analog prototype → bilinear transform → digital.
 * ══════════════════════════════════════════════════════════════════ */

/**
 * @brief Design a Butterworth lowpass filter.
 *
 * @param order   Filter order (1..2*MAX_SOS_SECTIONS)
 * @param cutoff  Normalised cutoff frequency, (0.0, 0.5) where 0.5 = Nyquist
 * @param sos     Output cascade (caller provides)
 * @return        0 on success, -1 on invalid parameters
 */
int butterworth_lowpass(int order, double cutoff, SOSCascade *sos);

/**
 * @brief Design a Butterworth highpass filter.
 *
 * Uses the lowpass-to-highpass frequency transformation:
 * substitute z⁻¹ → −z⁻¹ in the lowpass prototype.
 *
 * @param order   Filter order
 * @param cutoff  Normalised cutoff frequency (0.0, 0.5)
 * @param sos     Output cascade
 * @return        0 on success, -1 on invalid parameters
 */
int butterworth_highpass(int order, double cutoff, SOSCascade *sos);

/* ══════════════════════════════════════════════════════════════════
 *  IIR Filter Design — Chebyshev Type I
 *
 *  Chebyshev I filters allow passband ripple for a steeper rolloff:
 *    |H(jΩ)|² = 1 / (1 + ε² · T_N²(Ω/Ωc))
 *  where T_N is the Nth-order Chebyshev polynomial.
 * ══════════════════════════════════════════════════════════════════ */

/**
 * @brief Design a Chebyshev Type I lowpass filter.
 *
 * @param order     Filter order
 * @param ripple_db Maximum passband ripple in dB (e.g. 0.5, 1.0, 3.0)
 * @param cutoff    Normalised cutoff frequency (0.0, 0.5)
 * @param sos       Output cascade
 * @return          0 on success, -1 on invalid parameters
 */
int chebyshev1_lowpass(int order, double ripple_db, double cutoff,
                       SOSCascade *sos);

/* ══════════════════════════════════════════════════════════════════
 *  Frequency Response Evaluation
 *
 *  These compute H(e^{jω}) on n_points evenly spaced frequencies
 *  from ω = 0 to ω = π (DC to Nyquist).
 * ══════════════════════════════════════════════════════════════════ */

/**
 * @brief Compute frequency response of a general IIR transfer function.
 *
 * Evaluates H(z) = B(z)/A(z) at z = e^{jω} for n_points frequencies.
 *
 * @param b         Numerator coefficients
 * @param b_len     Length of b
 * @param a         Denominator coefficients
 * @param a_len     Length of a
 * @param mag       Output magnitude array, length n_points (caller allocates)
 * @param phase     Output phase array (radians), length n_points (caller allocates)
 * @param n_points  Number of frequency evaluation points (DC to Nyquist)
 */
void freq_response(const double *b, int b_len,
                   const double *a, int a_len,
                   double *mag, double *phase, int n_points);

/**
 * @brief Compute frequency response of a single biquad section.
 *
 * @param bq        Biquad coefficients
 * @param mag       Output magnitude, length n_points
 * @param phase     Output phase (radians), length n_points
 * @param n_points  Number of evaluation points
 */
void biquad_freq_response(const Biquad *bq,
                          double *mag, double *phase, int n_points);

/**
 * @brief Compute frequency response of a full SOS cascade.
 *
 * Multiplies the individual section responses together.
 *
 * @param sos       Cascade
 * @param mag       Output magnitude, length n_points
 * @param phase     Output phase (radians), length n_points
 * @param n_points  Number of evaluation points
 */
void sos_freq_response(const SOSCascade *sos,
                       double *mag, double *phase, int n_points);

/**
 * @brief Compute group delay at a single frequency.
 *
 * Group delay τ(ω) = −dφ/dω, approximated by central difference.
 *
 * @param b      Numerator coefficients
 * @param b_len  Length of b
 * @param a      Denominator coefficients
 * @param a_len  Length of a
 * @param omega  Angular frequency in radians (0 to π)
 * @return       Group delay in samples
 */
double group_delay_at(const double *b, int b_len,
                      const double *a, int a_len,
                      double omega);

#endif /* IIR_H */

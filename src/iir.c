/**
 * @file iir.c
 * @brief IIR digital filter implementation — processing, design, analysis.
 *
 * IMPLEMENTATION OVERVIEW:
 *
 *   1. General IIR filtering  — direct-form difference equation
 *   2. Biquad processing      — DF1 and DF2-Transposed per-sample
 *   3. Cascaded SOS            — chain of biquads for high-order filters
 *   4. Butterworth design      — maximally-flat magnitude response
 *   5. Chebyshev Type I design — equiripple passband, steeper rolloff
 *   6. Frequency response      — evaluate H(e^{jω}) on unit circle
 *
 * TUTORIAL CROSS-REFERENCES:
 *   Frequency response theory → chapters/06-frequency-response.md
 *   IIR design walkthrough    → chapters/11-iir-filter-design.md
 *   Filter structure details  → chapters/12-filter-structures.md
 *   FIR filters for contrast  → src/filter.c
 *
 * DESIGN METHOD (Butterworth / Chebyshev):
 *
 *   Step 1: Pre-warp the digital cutoff frequency to analog domain
 *           Ωa = 2·tan(ωd / 2)   where ωd = 2π·cutoff
 *
 *   Step 2: Place analog prototype poles (circle for Butterworth,
 *           ellipse for Chebyshev) scaled to Ωa
 *
 *   Step 3: Apply bilinear transform to each pole (or pole pair):
 *           z = (2 + s) / (2 − s)
 *
 *   Step 4: Convert z-domain pole pairs → biquad coefficients
 *
 *        ┌──────────────────────────────────────────────────┐
 *        │   Analog Domain        Digital Domain            │
 *        │                                                  │
 *        │   s-plane              z-plane                   │
 *        │   jΩ ↑                 Im ↑                      │
 *        │      |  ×p             |  ×                      │
 *        │   ───┼────→ σ      ───┼──(○)──→ Re              │
 *        │      |  ×p*            |  ×                      │
 *        │                        unit circle               │
 *        │                                                  │
 *        │   bilinear: s = 2(z−1)/(z+1)                    │
 *        │   maps  jΩ axis → unit circle                    │
 *        │   maps  left half → inside circle                │
 *        └──────────────────────────────────────────────────┘
 */

#define _GNU_SOURCE
#include "iir.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ══════════════════════════════════════════════════════════════════
 *  Section 1: General IIR Filtering
 *
 *  Direct-form difference equation:
 *
 *    a[0]·y[n] = b[0]·x[n] + b[1]·x[n-1] + ... + b[M]·x[n-M]
 *              − a[1]·y[n-1] − a[2]·y[n-2] − ... − a[N]·y[n-N]
 *
 *  When a_len == 1 (only a[0]), this reduces to an FIR filter.
 * ══════════════════════════════════════════════════════════════════ */

void iir_filter(const double *b, int b_len,
                const double *a, int a_len,
                const double *in, double *out, int n)
{
    double a0 = (a_len > 0) ? a[0] : 1.0;
    if (fabs(a0) < 1e-30) a0 = 1.0; /* safety */

    for (int i = 0; i < n; i++) {
        /* Feed-forward (numerator) sum: Σ b[k]·x[n-k] */
        double sum = 0.0;
        for (int k = 0; k < b_len; k++) {
            int idx = i - k;
            if (idx >= 0)
                sum += b[k] * in[idx];
        }

        /* Feedback (denominator) sum: − Σ a[k]·y[n-k], k=1..a_len-1 */
        for (int k = 1; k < a_len; k++) {
            int idx = i - k;
            if (idx >= 0)
                sum -= a[k] * out[idx];
        }

        out[i] = sum / a0;
    }
}

/**
 * Feed a unit impulse into the IIR system and record the response.
 * Useful for seeing how a filter "rings" and verifying stability.
 */
void iir_impulse_response(const double *b, int b_len,
                          const double *a, int a_len,
                          double *h, int n_out)
{
    double *impulse = (double *)calloc((size_t)n_out, sizeof(double));
    if (!impulse) return;
    impulse[0] = 1.0;

    iir_filter(b, b_len, a, a_len, impulse, h, n_out);
    free(impulse);
}

/* ══════════════════════════════════════════════════════════════════
 *  Section 2: Biquad Per-Sample Processing
 *
 *  A biquad is a 2nd-order IIR section. Processing one sample at
 *  a time is the basis for real-time audio and streaming DSP.
 * ══════════════════════════════════════════════════════════════════ */

void biquad_df1_init(BiquadDF1State *s)
{
    s->x1 = s->x2 = s->y1 = s->y2 = 0.0;
}

void biquad_df2t_init(BiquadDF2TState *s)
{
    s->s1 = s->s2 = 0.0;
}

/**
 * Direct Form I: straightforward implementation of the difference equation.
 *
 *   y[n] = b0·x[n] + b1·x[n-1] + b2·x[n-2]
 *                   - a1·y[n-1] - a2·y[n-2]
 *
 * Requires 4 state variables (x1, x2, y1, y2).
 * Advantage:  easy to understand, immune to limit cycles under some conditions.
 * Disadvantage: needs more memory than DF2T.
 */
double biquad_process_df1(const Biquad *bq, BiquadDF1State *s, double x)
{
    double y = bq->b0 * x + bq->b1 * s->x1 + bq->b2 * s->x2
                           - bq->a1 * s->y1 - bq->a2 * s->y2;

    /* Shift delay lines */
    s->x2 = s->x1;
    s->x1 = x;
    s->y2 = s->y1;
    s->y1 = y;

    return y;
}

/**
 * Direct Form II Transposed: minimum-state biquad.
 *
 *   y[n] = b0·x[n] + s1
 *   s1   = b1·x[n] - a1·y[n] + s2
 *   s2   = b2·x[n] - a2·y[n]
 *
 * Only 2 state variables (s1, s2).
 * Recommended for floating-point: better numerical properties for
 * narrow-band (high-Q) filters because the state variables store
 * differences rather than large accumulated sums.
 *
 * See chapters/12-filter-structures.md § "Direct Form II Transposed".
 */
double biquad_process_df2t(const Biquad *bq, BiquadDF2TState *s, double x)
{
    double y = bq->b0 * x + s->s1;

    s->s1 = bq->b1 * x - bq->a1 * y + s->s2;
    s->s2 = bq->b2 * x - bq->a2 * y;

    return y;
}

/**
 * Block processing wrapper around Direct Form I.
 * Processes n samples in one call — the state persists across calls
 * for streaming applications.
 */
void biquad_process_block(const Biquad *bq, BiquadDF1State *s,
                          const double *in, double *out, int n)
{
    for (int i = 0; i < n; i++) {
        out[i] = biquad_process_df1(bq, s, in[i]);
    }
}

/* ══════════════════════════════════════════════════════════════════
 *  Section 3: Cascaded SOS Processing
 *
 *  A high-order IIR filter is factored into K biquad sections:
 *
 *    x[n] → [SOS 0] → [SOS 1] → ··· → [SOS K−1] → ×gain → y[n]
 *
 *  This avoids coefficient-sensitivity and overflow problems that
 *  plague single-section high-order Direct-Form implementations.
 *  See O&S §6.8 and chapters/12-filter-structures.md.
 * ══════════════════════════════════════════════════════════════════ */

void sos_init(SOSCascade *sos)
{
    memset(sos, 0, sizeof(SOSCascade));
    sos->gain = 1.0;
}

double sos_process_sample(SOSCascade *sos, double x)
{
    double y = x;
    for (int i = 0; i < sos->n_sections; i++) {
        y = biquad_process_df1(&sos->sections[i], &sos->states[i], y);
    }
    return y * sos->gain;
}

void sos_process_block(SOSCascade *sos,
                       const double *in, double *out, int n)
{
    for (int i = 0; i < n; i++) {
        out[i] = sos_process_sample(sos, in[i]);
    }
}

/* ══════════════════════════════════════════════════════════════════
 *  Section 4: Butterworth Filter Design
 *
 *  Butterworth = "maximally flat" magnitude in the passband.
 *  All 2N−1 derivatives of |H(jΩ)|² are zero at Ω = 0.
 *
 *  Analog magnitude squared response:
 *    |H_a(jΩ)|² = 1 / (1 + (Ω/Ωc)^{2N})
 *
 *  Analog prototype poles lie on a circle of radius Ωc:
 *    s_k = Ωc · e^{j·π·(2k + N + 1) / (2N)},  k = 0..N−1
 *
 *  Only the left-half-plane poles (Re(s) < 0) are used.
 *
 *      Im(s)                      Butterworth poles (N=4)
 *       ↑                         on circle of radius Ωc
 *       |     ×                   ×  = pole location
 *       |  ×     ×                Only left-half used
 *   ────┼──────────→ Re(s)        for a stable filter
 *       |  ×     ×
 *       |     ×
 *       |
 * ══════════════════════════════════════════════════════════════════ */

int butterworth_lowpass(int order, double cutoff, SOSCascade *sos)
{
    if (order < 1 || order > 2 * MAX_SOS_SECTIONS ||
        cutoff <= 0.0 || cutoff >= 0.5)
        return -1;

    sos_init(sos);

    /* Step 1: Pre-warp the digital cutoff to analog frequency.
     *   ωd = 2π · cutoff   (digital angular frequency)
     *   Ωa = 2 · tan(ωd/2) (analog prototype frequency)
     *
     * The factor 2 comes from the bilinear transform with T = 1:
     *   s = 2(z−1)/(z+1)  →  Ωa = 2·tan(ω/2)
     */
    double wd = 2.0 * M_PI * cutoff;
    double Wa = 2.0 * tan(wd / 2.0);

    int n_pairs = order / 2;
    int has_real = (order % 2);

    /* Step 2 + 3: For each conjugate pole pair, compute analog pole
     * and apply bilinear transform to get biquad coefficients. */
    for (int k = 0; k < n_pairs; k++) {
        /* Pole angle on the unit circle (left-half plane guaranteed
         * because theta ∈ (π/2, π) for all valid k, N combinations) */
        double theta = M_PI * (2.0 * k + order + 1.0) / (2.0 * order);

        /* Analog pole: p = Wa · e^{jθ} */
        double sigma = Wa * cos(theta);   /* Real part (< 0 → stable) */
        double omega = Wa * sin(theta);   /* Imaginary part            */

        /* Bilinear transform: z = (2 + s) / (2 − s)
         *
         * For conjugate pair p, p*:
         *   z-domain denominator = (z − zp)(z − zp*)
         *                        = z² − 2·Re(zp)·z + |zp|²
         *
         * In z⁻¹ form (divide by z²):
         *   1 − 2·Re(zp)·z⁻¹ + |zp|²·z⁻²
         */
        double dre = 2.0 - sigma;           /* denominator real     */
        double dim = -omega;                 /* denominator imag     */
        double dsq = dre * dre + dim * dim;  /* |2 − p|²            */

        /* z-domain pole (complex): zp = (2 + p) / (2 − p) */
        double zp_re = ((2.0 + sigma) * dre + omega * dim) / dsq;
        double zp_im = (omega * dre - (2.0 + sigma) * dim) / dsq;

        /* Biquad denominator: 1 + a1·z⁻¹ + a2·z⁻² */
        double a1 = -2.0 * zp_re;
        double a2 = zp_re * zp_re + zp_im * zp_im;

        /* For lowpass, numerator zeros are at z = −1 (Nyquist null):
         *   (1 + z⁻¹)² = 1 + 2z⁻¹ + z⁻²
         *
         * Gain K chosen so H(z=1) = 1 for each section (unity DC gain):
         *   H(1) = K·(1+1)² / (1 + a1 + a2) = 4K / (1 + a1 + a2) = 1
         *   K = (1 + a1 + a2) / 4
         */
        double K = (1.0 + a1 + a2) / 4.0;

        Biquad *bq = &sos->sections[sos->n_sections];
        bq->b0 = K;
        bq->b1 = 2.0 * K;
        bq->b2 = K;
        bq->a1 = a1;
        bq->a2 = a2;
        biquad_df1_init(&sos->states[sos->n_sections]);
        sos->n_sections++;
    }

    /* Step 4: Handle real pole (odd order only).
     *
     * For odd N, there is one pole on the negative real axis: s = −Wa.
     * Bilinear: z = (2 + (−Wa)) / (2 − (−Wa)) = (2 − Wa) / (2 + Wa)
     *
     * First-order section stored as biquad with b2 = a2 = 0.
     */
    if (has_real) {
        double zp = (2.0 - Wa) / (2.0 + Wa);

        /* Numerator zero at z = −1 → (1 + z⁻¹)
         * DC gain: K·2 / (1 − zp) = 1  →  K = (1 − zp) / 2 */
        double K = (1.0 - zp) / 2.0;  /* actually 1 + (-zp) = 1 - zp */

        Biquad *bq = &sos->sections[sos->n_sections];
        bq->b0 = K;
        bq->b1 = K;
        bq->b2 = 0.0;
        bq->a1 = -zp;
        bq->a2 = 0.0;
        biquad_df1_init(&sos->states[sos->n_sections]);
        sos->n_sections++;
    }

    return 0;
}

/**
 * Butterworth highpass via the lowpass-to-highpass transform.
 *
 * In the analog domain this is Ω → Ωc²/Ω. After bilinear transform,
 * the effect is: replace z⁻¹ with −z⁻¹ in the lowpass transfer function.
 *
 * Practically: negate the odd-power numerator coefficients (b1) and
 * denominator coefficients (a1) of each biquad section.
 *
 * Numerator zeros move from z = −1 (lowpass) to z = +1 (highpass),
 * creating a DC null instead of a Nyquist null.
 */
int butterworth_highpass(int order, double cutoff, SOSCascade *sos)
{
    /* Design the complementary lowpass first */
    double lp_cutoff = 0.5 - cutoff;
    if (lp_cutoff <= 0.0 || lp_cutoff >= 0.5)
        return -1;

    /* Actually, the correct approach is to pre-warp the highpass cutoff
     * and use the analog highpass prototype directly. */
    if (order < 1 || order > 2 * MAX_SOS_SECTIONS ||
        cutoff <= 0.0 || cutoff >= 0.5)
        return -1;

    sos_init(sos);

    /* Pre-warp the highpass cutoff */
    double wd = 2.0 * M_PI * cutoff;
    double Wa = 2.0 * tan(wd / 2.0);

    int n_pairs = order / 2;
    int has_real = (order % 2);

    for (int k = 0; k < n_pairs; k++) {
        double theta = M_PI * (2.0 * k + order + 1.0) / (2.0 * order);

        /* Analog prototype pole (same as lowpass): p = Wa · e^{jθ} */
        double sigma = Wa * cos(theta);
        double omega = Wa * sin(theta);

        /* For highpass, invert the analog pole: p_hp = Wa² / p_lp
         * This maps lowpass poles to highpass poles. */
        double p_sq = sigma * sigma + omega * omega;
        double hp_sigma = Wa * Wa * sigma / p_sq;
        double hp_omega = Wa * Wa * omega / p_sq;  /* Note: negated for conjugation */

        /* Bilinear transform */
        double dre = 2.0 - hp_sigma;
        double dim = -hp_omega;
        double dsq = dre * dre + dim * dim;

        double zp_re = ((2.0 + hp_sigma) * dre + hp_omega * dim) / dsq;
        double zp_im = (hp_omega * dre - (2.0 + hp_sigma) * dim) / dsq;

        double a1 = -2.0 * zp_re;
        double a2 = zp_re * zp_re + zp_im * zp_im;

        /* Highpass: numerator zeros at z = +1 (DC null):
         *   (1 − z⁻¹)² = 1 − 2z⁻¹ + z⁻²
         *
         * Gain K: at Nyquist (z = −1), H should = 1:
         *   H(−1) = K·(1−(−1))² / (1 + a1·(−1) + a2) = 4K / (1 − a1 + a2)
         *   K = (1 − a1 + a2) / 4
         */
        double K = (1.0 - a1 + a2) / 4.0;

        Biquad *bq = &sos->sections[sos->n_sections];
        bq->b0 = K;
        bq->b1 = -2.0 * K;
        bq->b2 = K;
        bq->a1 = a1;
        bq->a2 = a2;
        biquad_df1_init(&sos->states[sos->n_sections]);
        sos->n_sections++;
    }

    if (has_real) {
        /* Real pole: p = −Wa, highpass inversion: p_hp = −Wa²/Wa = −Wa
         * Actually for highpass real pole: p_hp = Wa (positive real → right half?)
         * No — for highpass from lowpass prototype, use Ω → Ωc/Ω:
         *   p_hp = Wa²/(−Wa) = −Wa  (same pole for the prototype)
         * But we need the *highpass* first-order section. */
        double zp = (2.0 - Wa) / (2.0 + Wa);

        /* Highpass: zero at z = +1, DC null.
         * K·(1 − (−1)) / (1 − (−zp)) = K·2 / (1 + zp) → should be 1 at Nyquist
         * K = (1 + zp) / 2 */
        double K = (1.0 + zp) / 2.0;

        Biquad *bq = &sos->sections[sos->n_sections];
        bq->b0 = K;
        bq->b1 = -K;
        bq->b2 = 0.0;
        bq->a1 = -zp;
        bq->a2 = 0.0;
        biquad_df1_init(&sos->states[sos->n_sections]);
        sos->n_sections++;
    }

    return 0;
}

/* ══════════════════════════════════════════════════════════════════
 *  Section 5: Chebyshev Type I Filter Design
 *
 *  Chebyshev I has equiripple in the passband and monotonic stopband.
 *  The poles lie on an ellipse (not a circle like Butterworth).
 *
 *  Steps:
 *    1. Compute ripple parameter ε from ripple_db
 *    2. Compute the "ellipse factor" a = (1/N)·asinh(1/ε)
 *    3. Analog poles: s_k = Ωc·(−sinh(a)sin(θ_k) + j·cosh(a)cos(θ_k))
 *    4. Apply bilinear transform (same as Butterworth)
 *
 *       Im(s)
 *        ↑          Chebyshev poles (N=4)
 *        |     ×    on ELLIPSE
 *        |  ×    ×  × = pole
 *    ────┼──────────→ Re(s)
 *        |  ×    ×  semi-minor = Ωc·sinh(a)
 *        |     ×    semi-major = Ωc·cosh(a)
 * ══════════════════════════════════════════════════════════════════ */

int chebyshev1_lowpass(int order, double ripple_db, double cutoff,
                       SOSCascade *sos)
{
    if (order < 1 || order > 2 * MAX_SOS_SECTIONS ||
        cutoff <= 0.0 || cutoff >= 0.5 || ripple_db <= 0.0)
        return -1;

    sos_init(sos);

    /* Ripple parameter: ε = √(10^{R/10} − 1) */
    double eps = sqrt(pow(10.0, ripple_db / 10.0) - 1.0);

    /* Ellipse factor: a = (1/N) · arcsinh(1/ε) */
    double a = asinh(1.0 / eps) / (double)order;

    /* Pre-warp cutoff */
    double wd = 2.0 * M_PI * cutoff;
    double Wa = 2.0 * tan(wd / 2.0);

    int n_pairs = order / 2;
    int has_real = (order % 2);

    for (int k = 0; k < n_pairs; k++) {
        /* Pole angle (same pattern as Butterworth) */
        double theta = M_PI * (2.0 * k + 1.0) / (2.0 * order);

        /* Analog Chebyshev pole on ellipse:
         *   Re(s) = −Wa · sinh(a) · sin(θ)     (semi-minor axis)
         *   Im(s) =  Wa · cosh(a) · cos(θ)     (semi-major axis)
         */
        double sigma = -Wa * sinh(a) * sin(theta);
        double omega =  Wa * cosh(a) * cos(theta);

        /* Bilinear transform: z = (2 + s)/(2 − s) */
        double dre = 2.0 - sigma;
        double dim = -omega;
        double dsq = dre * dre + dim * dim;

        double zp_re = ((2.0 + sigma) * dre + omega * dim) / dsq;
        double zp_im = (omega * dre - (2.0 + sigma) * dim) / dsq;

        double a1 = -2.0 * zp_re;
        double a2 = zp_re * zp_re + zp_im * zp_im;

        /* Lowpass numerator zeros at z = −1 */
        double K = (1.0 + a1 + a2) / 4.0;

        Biquad *bq = &sos->sections[sos->n_sections];
        bq->b0 = K;
        bq->b1 = 2.0 * K;
        bq->b2 = K;
        bq->a1 = a1;
        bq->a2 = a2;
        biquad_df1_init(&sos->states[sos->n_sections]);
        sos->n_sections++;
    }

    if (has_real) {
        /* Real pole: s = −Wa·sinh(a) */
        double sp = -Wa * sinh(a);
        double zp = (2.0 + sp) / (2.0 - sp);

        double K = (1.0 + (-zp)) / 2.0;

        Biquad *bq = &sos->sections[sos->n_sections];
        bq->b0 = K;
        bq->b1 = K;
        bq->b2 = 0.0;
        bq->a1 = -zp;
        bq->a2 = 0.0;
        biquad_df1_init(&sos->states[sos->n_sections]);
        sos->n_sections++;
    }

    /* Chebyshev gain correction:
     * The passband ripple means the DC gain is NOT automatically 1.
     * For even order, the gain at DC equals 1/√(1+ε²).
     * We correct so the cascade has unity peak (at DC for odd, or at
     * the ripple maximum for even order). */
    if (!has_real) {
        /* Even order: scale so DC gain = 1/√(1+ε²) (passband edge) */
        /* Actually, we'll just normalise DC to match expected level */
        double dc_gain = 1.0;
        for (int i = 0; i < sos->n_sections; i++) {
            Biquad *bq = &sos->sections[i];
            double section_dc = (bq->b0 + bq->b1 + bq->b2) /
                                (1.0 + bq->a1 + bq->a2);
            dc_gain *= section_dc;
        }
        if (fabs(dc_gain) > 1e-30) {
            /* For even-order Chebyshev, DC should be 1/√(1+ε²) */
            double target_dc = 1.0 / sqrt(1.0 + eps * eps);
            sos->gain = target_dc / dc_gain;
        }
    }

    return 0;
}

/* ══════════════════════════════════════════════════════════════════
 *  Section 6: Frequency Response Evaluation
 *
 *  The frequency response of H(z) is obtained by evaluating on
 *  the unit circle: z = e^{jω} for ω ∈ [0, π].
 *
 *  For H(z) = B(z)/A(z):
 *    B(e^{jω}) = Σ b[k] · e^{−jkω}
 *    A(e^{jω}) = Σ a[k] · e^{−jkω}
 *    H(e^{jω}) = B(e^{jω}) / A(e^{jω})
 * ══════════════════════════════════════════════════════════════════ */

/**
 * Evaluate polynomial P(z) = p[0] + p[1]·z⁻¹ + ... + p[n-1]·z^{-(n-1)}
 * at z = e^{jω}.
 */
static Complex eval_poly_on_circle(const double *p, int len, double omega)
{
    Complex result = {0.0, 0.0};
    for (int k = 0; k < len; k++) {
        /* e^{−jkω} = cos(kω) − j·sin(kω) */
        result.re += p[k] * cos(k * omega);
        result.im -= p[k] * sin(k * omega);
    }
    return result;
}

void freq_response(const double *b, int b_len,
                   const double *a, int a_len,
                   double *mag, double *phase, int n_points)
{
    for (int i = 0; i < n_points; i++) {
        double omega = M_PI * (double)i / (double)(n_points - 1);

        Complex B = eval_poly_on_circle(b, b_len, omega);
        Complex A = {1.0, 0.0};
        if (a_len > 0)
            A = eval_poly_on_circle(a, a_len, omega);

        /* H = B / A */
        double a_mag_sq = A.re * A.re + A.im * A.im;
        Complex H;
        if (a_mag_sq < 1e-30) {
            H.re = 1e10; H.im = 0.0;
        } else {
            H.re = (B.re * A.re + B.im * A.im) / a_mag_sq;
            H.im = (B.im * A.re - B.re * A.im) / a_mag_sq;
        }

        if (mag)   mag[i]   = complex_mag(H);
        if (phase) phase[i] = complex_phase(H);
    }
}

void biquad_freq_response(const Biquad *bq,
                          double *mag, double *phase, int n_points)
{
    double b[3] = {bq->b0, bq->b1, bq->b2};
    double a[3] = {1.0, bq->a1, bq->a2};
    freq_response(b, 3, a, 3, mag, phase, n_points);
}

void sos_freq_response(const SOSCascade *sos,
                       double *mag, double *phase, int n_points)
{
    /* Initialise to gain only */
    for (int i = 0; i < n_points; i++) {
        mag[i]   = fabs(sos->gain);
        phase[i] = (sos->gain < 0) ? M_PI : 0.0;
    }

    /* Multiply each section's response */
    double *sec_mag   = (double *)calloc((size_t)n_points, sizeof(double));
    double *sec_phase = (double *)calloc((size_t)n_points, sizeof(double));
    if (!sec_mag || !sec_phase) {
        free(sec_mag);
        free(sec_phase);
        return;
    }

    for (int s = 0; s < sos->n_sections; s++) {
        biquad_freq_response(&sos->sections[s], sec_mag, sec_phase, n_points);
        for (int i = 0; i < n_points; i++) {
            mag[i]   *= sec_mag[i];
            phase[i] += sec_phase[i];
        }
    }

    free(sec_mag);
    free(sec_phase);
}

/**
 * Approximate group delay at a single frequency using
 * central finite difference:
 *
 *   τ(ω) = −dφ/dω ≈ −(φ(ω+Δ) − φ(ω−Δ)) / (2Δ)
 */
double group_delay_at(const double *b, int b_len,
                      const double *a, int a_len,
                      double omega)
{
    double delta = 1e-5;
    double phase_plus, phase_minus;

    /* Evaluate phase at ω + δ */
    freq_response(b, b_len, a, a_len, NULL, &phase_plus, 1);
    /* Actually need to evaluate at specific omega, not at index 0 */

    /* Direct evaluation at specific omega values */
    Complex B_plus  = eval_poly_on_circle(b, b_len, omega + delta);
    Complex A_plus  = (a_len > 0) ? eval_poly_on_circle(a, a_len, omega + delta)
                                  : (Complex){1.0, 0.0};
    double a_sq = A_plus.re * A_plus.re + A_plus.im * A_plus.im;
    Complex H_plus;
    H_plus.re = (B_plus.re * A_plus.re + B_plus.im * A_plus.im) / (a_sq + 1e-30);
    H_plus.im = (B_plus.im * A_plus.re - B_plus.re * A_plus.im) / (a_sq + 1e-30);
    phase_plus = complex_phase(H_plus);

    Complex B_minus = eval_poly_on_circle(b, b_len, omega - delta);
    Complex A_minus = (a_len > 0) ? eval_poly_on_circle(a, a_len, omega - delta)
                                  : (Complex){1.0, 0.0};
    a_sq = A_minus.re * A_minus.re + A_minus.im * A_minus.im;
    Complex H_minus;
    H_minus.re = (B_minus.re * A_minus.re + B_minus.im * A_minus.im) / (a_sq + 1e-30);
    H_minus.im = (B_minus.im * A_minus.re - B_minus.re * A_minus.im) / (a_sq + 1e-30);
    phase_minus = complex_phase(H_minus);

    /* Unwrap phase difference */
    double dp = phase_plus - phase_minus;
    while (dp > M_PI) dp -= 2.0 * M_PI;
    while (dp < -M_PI) dp += 2.0 * M_PI;

    return -dp / (2.0 * delta);
}

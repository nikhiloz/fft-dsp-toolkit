/**
 * @file 06-frequency-response.c
 * @brief Chapter 6 Demo: Frequency Response, Poles & Zeros
 *
 * Interactive demonstrations:
 *   1. Frequency response of a simple 2-point average (FIR lowpass)
 *   2. Magnitude and phase response — FIR vs IIR
 *   3. Pole-zero placement controls the frequency response shape
 *   4. Effect of pole radius on resonance (Q factor)
 *   5. Group delay — FIR (constant) vs IIR (frequency-dependent)
 *   6. Stability from the pole-zero plot
 *   7. All-pass filter — unity magnitude, phase-only modification
 *
 * Build:  make release && ./build/bin/ch06
 *
 * References:
 *   - Oppenheim & Schafer, Discrete-Time Signal Processing, Ch 2-3, 5
 *   - Proakis & Manolakis, Digital Signal Processing, Ch 3-4
 *   - Lyons, Understanding DSP, Ch 5-7
 *
 * Prerequisites: Ch 03 (Complex Numbers), Ch 05 (Z-Transform)
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dsp_utils.h"
#include "fft.h"
#include "iir.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ── Helpers ─────────────────────────────────────────────────────── */

static void print_separator(const char *title)
{
    printf("\n╔══════════════════════════════════════════════════════════╗\n");
    printf("║  %-56s║\n", title);
    printf("╚══════════════════════════════════════════════════════════╝\n");
}

/**
 * Evaluate transfer function H(z) = B(z)/A(z) at a point on the unit circle.
 * z = e^{jω} = cos(ω) + j·sin(ω)
 */
static Complex eval_hz(const double *b, int b_len,
                       const double *a, int a_len,
                       double omega)
{
    /* Numerator: B(e^{jω}) = Σ b[k]·e^{-jkω} */
    Complex num = {0.0, 0.0};
    for (int k = 0; k < b_len; k++) {
        num.re += b[k] * cos(k * omega);
        num.im -= b[k] * sin(k * omega);
    }

    /* Denominator: A(e^{jω}) = 1 + Σ a[k]·e^{-jkω}, k=1.. */
    Complex den = {1.0, 0.0};
    for (int k = 1; k < a_len; k++) {
        den.re += a[k] * cos(k * omega);
        den.im -= a[k] * sin(k * omega);
    }

    /* H = num / den */
    double den_sq = den.re * den.re + den.im * den.im;
    if (den_sq < 1e-30) {
        Complex inf = {1e6, 0.0};
        return inf;
    }
    Complex H;
    H.re = (num.re * den.re + num.im * den.im) / den_sq;
    H.im = (num.im * den.re - num.re * den.im) / den_sq;
    return H;
}

/**
 * Print an ASCII magnitude response plot.
 * bar_width controls the scale of the horizontal bars.
 */
static void plot_mag_response(const char *label,
                              const double *b, int b_len,
                              const double *a, int a_len,
                              int n_points)
{
    printf("\n  %s — Magnitude Response:\n", label);
    printf("  %6s  %8s  %8s  %s\n", "f/fs", "ω/π", "|H| dB", "");

    for (int i = 0; i <= n_points; i++) {
        double omega = M_PI * (double)i / (double)n_points;
        Complex H = eval_hz(b, b_len, a, a_len, omega);
        double mag = complex_mag(H);
        double mag_db = (mag > 1e-10) ? 20.0 * log10(mag) : -80.0;

        /* Clamp for display */
        if (mag_db < -60.0) mag_db = -60.0;
        if (mag_db > 20.0) mag_db = 20.0;

        /* ASCII bar: map [-60, +20] → [0, 40] chars */
        int bar = (int)((mag_db + 60.0) / 80.0 * 40.0);
        if (bar < 0) bar = 0;
        if (bar > 40) bar = 40;

        printf("  %6.3f  %8.3f  %+8.1f  ",
               (double)i / (2.0 * n_points),
               (double)i / (double)n_points,
               mag_db);
        for (int j = 0; j < bar; j++) printf("█");
        printf("\n");
    }
}

/* ── Demo 1: Simple FIR lowpass frequency response ───────────────── */

static void demo_fir_response(void)
{
    print_separator("Demo 1: Frequency Response of 2-Point Average");

    /*
     * H(z) = (1 + z⁻¹) / 2 = 0.5 + 0.5·z⁻¹
     *
     * This is the simplest lowpass FIR.
     * Evaluating on the unit circle z = e^{jω}:
     *
     *   H(e^{jω}) = 0.5·(1 + e^{-jω})
     *             = 0.5 · e^{-jω/2} · 2·cos(ω/2)
     *             = cos(ω/2) · e^{-jω/2}
     *
     *   |H(e^{jω})| = |cos(ω/2)|
     *   ∠H(e^{jω})  = −ω/2         (linear phase!)
     *
     *        Magnitude                Phase
     *   1.0 |████████           0 |────────
     *       |  ████████           |  ──────
     *   0.5 |    ████████    −π/4 |    ────
     *       |      ████████       |      ──
     *   0.0 |────────────── −π/2 |────────
     *       DC     f_s/4   Nyquist
     */

    printf("\n  H(z) = (1 + z⁻¹) / 2    (2-point average → lowpass)\n\n");
    printf("  Closed-form: |H(e^{jω})| = |cos(ω/2)|\n");
    printf("               ∠H(e^{jω}) = −ω/2  (linear phase!)\n\n");

    double b[] = {0.5, 0.5};
    int b_len = 2;
    double a[] = {1.0};        /* FIR: denominator = 1 */
    int a_len = 1;

    printf("  %6s  %8s  %10s  %10s  %10s\n",
           "f/fs", "ω/π", "|H|", "|cos(ω/2)|", "|H| dB");

    for (int i = 0; i <= 16; i++) {
        double omega = M_PI * (double)i / 16.0;
        Complex H = eval_hz(b, b_len, a, a_len, omega);
        double mag = complex_mag(H);
        double expected = fabs(cos(omega / 2.0));
        double mag_db = (mag > 1e-10) ? 20.0 * log10(mag) : -100.0;

        printf("  %6.3f  %8.3f  %10.4f  %10.4f  %+10.1f",
               (double)i / 32.0, (double)i / 16.0,
               mag, expected, mag_db);
        if (fabs(mag - expected) < 0.001) printf("  ✓");
        printf("\n");
    }
    printf("\n  ✓ Computed magnitude matches closed-form |cos(ω/2)| perfectly.\n");
    printf("  ✓ Gain = 0 dB at DC, −∞ dB at Nyquist → lowpass!\n");
}

/* ── Demo 2: FIR vs IIR magnitude/phase comparison ───────────────── */

static void demo_fir_vs_iir(void)
{
    print_separator("Demo 2: Magnitude & Phase — FIR vs IIR");

    /* FIR: 5-tap moving average */
    double b_fir[] = {0.2, 0.2, 0.2, 0.2, 0.2};
    double a_fir[] = {1.0};

    /* IIR: 1st-order lowpass y[n] = 0.1·x[n] + 0.9·y[n-1] */
    double b_iir[] = {0.1};
    double a_iir[] = {1.0, -0.9};

    printf("\n  FIR: h = {0.2, 0.2, 0.2, 0.2, 0.2}  (5-tap average)\n");
    printf("  IIR: H(z) = 0.1 / (1 − 0.9z⁻¹)    (1st-order, pole at 0.9)\n\n");

    printf("  %6s  %10s %10s  %10s %10s\n",
           "ω/π", "|FIR| dB", "∠FIR °", "|IIR| dB", "∠IIR °");

    for (int i = 0; i <= 16; i++) {
        double omega = M_PI * (double)i / 16.0;

        Complex H_fir = eval_hz(b_fir, 5, a_fir, 1, omega);
        Complex H_iir = eval_hz(b_iir, 1, a_iir, 2, omega);

        double mag_fir_db = 20.0 * log10(complex_mag(H_fir) + 1e-30);
        double mag_iir_db = 20.0 * log10(complex_mag(H_iir) + 1e-30);
        double phase_fir = complex_phase(H_fir) * 180.0 / M_PI;
        double phase_iir = complex_phase(H_iir) * 180.0 / M_PI;

        printf("  %6.3f  %+10.1f %+10.1f  %+10.1f %+10.1f\n",
               (double)i / 16.0,
               mag_fir_db, phase_fir, mag_iir_db, phase_iir);
    }

    printf("\n  Key observations:\n");
    printf("  ✓ FIR has LINEAR phase — constant group delay.\n");
    printf("  ✓ IIR has NON-LINEAR phase — group delay varies with frequency.\n");
    printf("  ✓ IIR achieves steep rolloff with fewer coefficients.\n");
}

/* ── Demo 3: Pole-zero placement → frequency response ────────────── */

static void demo_pole_zero_effect(void)
{
    print_separator("Demo 3: Pole-Zero Placement Shapes the Response");

    printf("\n  Zeros PULL the response DOWN (create nulls).\n");
    printf("  Poles PUSH the response UP  (create peaks).\n\n");

    /*
     *  z-plane pole-zero geometry:
     *
     *        Im(z)
     *         ↑
     *     ○   |           ○ = zero
     *         |           × = pole
     *   ──────(○)────→ Re(z)
     *     ×   |           ( ) = unit circle
     *         |
     *
     *  Rule of thumb:
     *    - A zero ON the unit circle → complete null at that frequency
     *    - A pole NEAR the unit circle → sharp peak at that frequency
     *    - Distance from unit circle → width of the null/peak
     */

    struct {
        const char *name;
        double b[3]; int bl;
        double a[3]; int al;
        const char *desc;
    } systems[] = {
        {"Zero at z=1 (DC null)",
         {1.0, -1.0, 0.0}, 2, {1.0, 0.0, 0.0}, 1,
         "Highpass: zero blocks DC"},
        {"Zero at z=-1 (Nyquist null)",
         {1.0, 1.0, 0.0}, 2, {1.0, 0.0, 0.0}, 1,
         "Lowpass: zero blocks Nyquist"},
        {"Pole at z=0.9 (near DC)",
         {0.1, 0.0, 0.0}, 1, {1.0, -0.9, 0.0}, 2,
         "Boosts low frequencies (lowpass)"},
        {"Pole at z=-0.9 (near Nyquist)",
         {0.1, 0.0, 0.0}, 1, {1.0, 0.9, 0.0}, 2,
         "Boosts high frequencies (highpass)"},
    };
    int n_sys = 4;

    for (int s = 0; s < n_sys; s++) {
        printf("  [%d] %s\n      %s\n", s + 1, systems[s].name, systems[s].desc);
        plot_mag_response(systems[s].name,
                          systems[s].b, systems[s].bl,
                          systems[s].a, systems[s].al, 16);
    }
}

/* ── Demo 4: Pole radius and Q factor ────────────────────────────── */

static void demo_pole_radius(void)
{
    print_separator("Demo 4: Pole Radius Controls Resonance (Q Factor)");

    /*
     * Complex conjugate pole pair at angle θ = π/4 (= fs/8).
     * H(z) = 1 / (1 − 2r·cos(θ)·z⁻¹ + r²·z⁻²)
     *
     * As r → 1, the peak at θ gets sharper (higher Q).
     *
     *     r = 0.5     r = 0.9     r = 0.99
     *      ___          /\          │
     *     /   \        /  \         │
     *    /     \      /    \        │
     *   ────────    ─/──────\─   ──┤├──
     */

    double theta = M_PI / 4.0;
    double radii[] = {0.5, 0.8, 0.95, 0.99};
    int n_radii = 4;

    printf("\n  Resonant frequency: ω = π/4 (f = fs/8)\n");
    printf("  Varying pole radius r:\n");

    for (int ri = 0; ri < n_radii; ri++) {
        double r = radii[ri];

        double a[] = {1.0, -2.0 * r * cos(theta), r * r};
        double b[] = {1.0};

        printf("\n  r = %.2f  (distance from unit circle: %.4f):\n",
               r, 1.0 - r);

        double peak_db = -200.0;
        for (int i = 0; i <= 16; i++) {
            double w = M_PI * (double)i / 16.0;
            Complex H = eval_hz(b, 1, a, 3, w);
            double mag_db = 20.0 * log10(complex_mag(H) + 1e-30);
            if (mag_db > peak_db) peak_db = mag_db;

            if (mag_db < -40.0) mag_db = -40.0;
            int bar = (int)((mag_db + 40.0) / 60.0 * 30.0);
            if (bar < 0) bar = 0;
            if (bar > 30) bar = 30;

            char marker = (i == 4) ? '*' : ' ';
            printf("    %5.3f  %+7.1f dB  ", (double)i / 16.0, mag_db);
            for (int j = 0; j < bar; j++) printf("█");
            printf("%c\n", marker);
        }
        printf("    Peak gain: %+.1f dB\n", peak_db);
    }
    printf("\n  * = resonant frequency (π/4)\n");
    printf("  ✓ As r → 1, peak gain ↑ and bandwidth ↓ (higher Q).\n");
}

/* ── Demo 5: Group delay — FIR vs IIR ────────────────────────────── */

static void demo_group_delay(void)
{
    print_separator("Demo 5: Group Delay — FIR (Constant) vs IIR");

    /*
     * Group delay: τ(ω) = −dφ(ω)/dω
     *
     * FIR with symmetric coefficients → linear phase → constant τ
     * IIR → non-linear phase → τ varies with frequency
     *
     *   FIR τ(ω):          IIR τ(ω):
     *                          /\
     *   ─────────────      ───/──\───
     *   τ = (N-1)/2        varies!
     */

    /* FIR: symmetric 5-tap */
    double b_fir[] = {0.1, 0.25, 0.3, 0.25, 0.1};
    double a_unit[] = {1.0};

    /* IIR: 2nd-order lowpass (pole at 0.8, angle π/4) */
    double b_iir[] = {0.05};
    double a_iir[] = {1.0, -2.0 * 0.8 * cos(M_PI / 4.0), 0.64};

    printf("\n  FIR: h = {0.1, 0.25, 0.3, 0.25, 0.1}  (symmetric, linear phase)\n");
    printf("       Expected group delay = (5−1)/2 = 2.0 samples (constant)\n\n");
    printf("  IIR: H(z) = 0.05 / (1 − 1.131z⁻¹ + 0.64z⁻²)\n");
    printf("       Group delay varies with frequency!\n\n");

    printf("  %6s  %12s  %12s\n", "ω/π", "τ_FIR", "τ_IIR");

    for (int i = 1; i < 16; i++) {
        double omega = M_PI * (double)i / 16.0;

        double gd_fir = group_delay_at(b_fir, 5, a_unit, 1, omega);
        double gd_iir = group_delay_at(b_iir, 1, a_iir, 3, omega);

        printf("  %6.3f  %12.4f  %12.4f", (double)i / 16.0, gd_fir, gd_iir);
        if (fabs(gd_fir - 2.0) < 0.05) printf("  ✓ FIR constant");
        printf("\n");
    }

    printf("\n  ✓ FIR group delay ≈ 2.0 at all frequencies (linear phase).\n");
    printf("  ✓ IIR group delay varies — peaks near the resonant frequency.\n");
    printf("  ✓ Linear phase preserves waveform shape (no dispersion).\n");
}

/* ── Demo 6: Stability from pole locations ───────────────────────── */

static void demo_stability(void)
{
    print_separator("Demo 6: Stability Analysis from Pole Locations");

    /*
     * BIBO Stability Rule (discrete-time):
     *
     *   All poles of H(z) must lie STRICTLY INSIDE the
     *   unit circle (|z| < 1) for the system to be stable.
     *
     *     Im(z)
     *      ↑    Stable region: inside
     *      |   ┌─────────────────┐
     *      |   │    ×  × stable  │
     *      |───│──(unit circle)──│──→ Re(z)
     *      |   │    ×  × stable  │
     *      |   └─────────────────┘
     *      |         × unstable (outside)
     */

    struct {
        const char *name;
        double a[3];
        int a_len;
        const char *verdict;
    } systems[] = {
        {"Single pole at z = 0.5",
         {1.0, -0.5, 0.0}, 2, "✓ STABLE (|0.5| < 1)"},
        {"Single pole at z = 0.99",
         {1.0, -0.99, 0.0}, 2, "✓ STABLE but near margin (|0.99| < 1)"},
        {"Single pole at z = 1.0",
         {1.0, -1.0, 0.0}, 2, "⚠ MARGINALLY STABLE (|1.0| = 1, integrator)"},
        {"Single pole at z = 1.01",
         {1.0, -1.01, 0.0}, 2, "✗ UNSTABLE (|1.01| > 1)"},
        {"Conjugate pair |z| = 0.9",
         {1.0, -1.272, 0.81}, 3, "✓ STABLE (|0.9| < 1)"},
    };
    int n_sys = 5;

    printf("\n  Rule: ALL poles must have |z| < 1 for BIBO stability.\n\n");
    printf("  %-35s  %-8s  %s\n", "System", "Max|pole|", "Verdict");
    printf("  %-35s  %-8s  %s\n",
           "───────────────────────────────────",
           "────────", "───────────────────");

    for (int s = 0; s < n_sys; s++) {
        /* Compute impulse response to show stability visually */
        double b[] = {1.0};
        double h[32];
        iir_impulse_response(b, 1, systems[s].a, systems[s].a_len, h, 32);

        /* Find max |pole| from the impulse response trend */
        double max_h = 0.0;
        for (int i = 0; i < 32; i++) {
            if (fabs(h[i]) > max_h) max_h = fabs(h[i]);
        }

        /* Estimate pole radius from coefficients */
        double pole_r;
        if (systems[s].a_len == 2) {
            pole_r = fabs(systems[s].a[1]);
        } else {
            pole_r = sqrt(fabs(systems[s].a[2]));
        }

        printf("  %-35s  %-8.4f  %s\n",
               systems[s].name, pole_r, systems[s].verdict);

        /* Show first few impulse response samples */
        printf("    h[n]: ");
        for (int i = 0; i < 8; i++) {
            printf("%+.3f ", h[i]);
        }
        printf("...\n");
        double h_end = fabs(h[31]);
        if (h_end < 0.01) printf("    → decays to ~0 ✓\n");
        else if (h_end > 100) printf("    → diverges! ✗\n");
        else printf("    → |h[31]| = %.3f\n", h_end);
    }
}

/* ── Demo 7: All-pass filter ─────────────────────────────────────── */

static void demo_allpass(void)
{
    print_separator("Demo 7: All-Pass Filter — Unity Magnitude, Phase Only");

    /*
     * First-order all-pass:
     *   H(z) = (a + z⁻¹) / (1 + a·z⁻¹)
     *
     * The zero at z = −a is the RECIPROCAL (and conjugate) of the
     * pole at z = −a. This forces |H(e^{jω})| = 1 for all ω.
     *
     * All-pass filters are used for:
     *   - Phase equalization
     *   - Building lattice structures
     *   - Creating minimum-phase systems
     */

    double a_val = 0.5;  /* pole at z = −0.5, zero at z = −2.0 */

    printf("\n  H(z) = (%.1f + z⁻¹) / (1 + %.1f·z⁻¹)\n", a_val, a_val);
    printf("  Pole at z = −%.1f, zero at z = −%.1f\n", a_val, 1.0 / a_val);
    printf("  (zero is the reciprocal of the pole → all-pass)\n\n");

    double b[] = {a_val, 1.0};
    double a[] = {1.0, a_val};

    printf("  %6s  %10s  %10s  %s\n", "ω/π", "|H|", "∠H (deg)", "");

    for (int i = 0; i <= 16; i++) {
        double omega = M_PI * (double)i / 16.0;
        Complex H = eval_hz(b, 2, a, 2, omega);
        double mag = complex_mag(H);
        double phase_deg = complex_phase(H) * 180.0 / M_PI;

        printf("  %6.3f  %10.6f  %+10.1f  ", (double)i / 16.0, mag, phase_deg);
        if (fabs(mag - 1.0) < 0.001) printf("✓ |H|=1");
        printf("\n");
    }

    printf("\n  ✓ |H(e^{jω})| = 1.000 at ALL frequencies — magnitude is flat.\n");
    printf("  ✓ Phase varies from 0° to −180° — this is the useful part.\n");
    printf("  ✓ All-pass filters modify only phase, not magnitude.\n");
}

/* ── Main ────────────────────────────────────────────────────────── */

int main(void)
{
    printf("╔══════════════════════════════════════════════════════════╗\n");
    printf("║  Chapter 6: Frequency Response, Poles & Zeros          ║\n");
    printf("║  DSP Tutorial Suite — Interactive Demos                 ║\n");
    printf("╚══════════════════════════════════════════════════════════╝\n");

    demo_fir_response();
    demo_fir_vs_iir();
    demo_pole_zero_effect();
    demo_pole_radius();
    demo_group_delay();
    demo_stability();
    demo_allpass();

    printf("\n════════════════════════════════════════════════════════════\n");
    printf("  All demos complete. See chapters/06-frequency-response.md\n");
    printf("  for theory, derivations, and exercises.\n");
    printf("════════════════════════════════════════════════════════════\n\n");

    return 0;
}

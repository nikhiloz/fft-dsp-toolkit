/**
 * Chapter 5 Demo: The Z-Transform
 *
 * Interactive demonstrations:
 *   1. Z-Transform of an FIR filter (polynomial evaluation)
 *   2. Pole-zero map and stability analysis
 *   3. Frequency response from H(z) — magnitude and phase
 *   4. Effect of pole radius on resonance sharpness
 *   5. Comparison: FIR (all-zero) vs. IIR (poles + zeros)
 *   6. Difference equation → impulse response via Z-domain
 *   7. Time delay property: z^{-k} in action
 *
 * Build:  make chapters && ./build/bin/ch05s
 *
 * References:
 *   - Oppenheim & Willsky, Signals and Systems, Ch 10
 *   - Oppenheim & Schafer, Discrete-Time Signal Processing, Ch 3
 *   - Proakis & Manolakis, Digital Signal Processing, Ch 3
 *   - Lyons, Understanding DSP, Ch 6-7
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dsp_utils.h"
#include "signal_gen.h"
#include "convolution.h"

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
 * Evaluate polynomial B(z) = b[0] + b[1]*z^{-1} + ... + b[M]*z^{-M}
 * at a complex point z.
 */
static Complex eval_poly(const double *b, int len, Complex z)
{
    /* Compute z^{-1} */
    double mag_sq = z.re * z.re + z.im * z.im;
    Complex z_inv;
    if (mag_sq < 1e-30) {
        z_inv.re = 0.0; z_inv.im = 0.0; /* avoid division by zero */
    } else {
        z_inv.re =  z.re / mag_sq;
        z_inv.im = -z.im / mag_sq;
    }

    /* Horner's method in z^{-1}: B = b[M] * z^{-1} + b[M-1], iterate */
    Complex result;
    result.re = b[len - 1];
    result.im = 0.0;

    for (int i = len - 2; i >= 0; i--) {
        result = complex_mul(result, z_inv);
        result.re += b[i];
    }

    return result;
}

/**
 * Evaluate transfer function H(z) = B(z) / A(z)
 * where A(z) = 1 + a[1]*z^{-1} + ... (a[0] is implicitly 1)
 * and B(z) = b[0] + b[1]*z^{-1} + ...
 */
static Complex eval_h(const double *b, int b_len,
                       const double *a, int a_len,
                       Complex z)
{
    Complex num = eval_poly(b, b_len, z);

    /* For FIR (a_len == 0), denominator = 1 */
    if (a_len <= 0) return num;

    /* Build full denominator: 1 + a[0]*z^{-1} + a[1]*z^{-2} + ... */
    double a_full[32];
    a_full[0] = 1.0;
    for (int i = 0; i < a_len && i < 31; i++) {
        a_full[i + 1] = a[i];
    }
    Complex den = eval_poly(a_full, a_len + 1, z);

    /* H(z) = num / den */
    double den_mag_sq = den.re * den.re + den.im * den.im;
    if (den_mag_sq < 1e-30) {
        Complex inf_val = {1e10, 0.0};
        return inf_val;
    }
    Complex result;
    result.re = (num.re * den.re + num.im * den.im) / den_mag_sq;
    result.im = (num.im * den.re - num.re * den.im) / den_mag_sq;
    return result;
}

/* ── Demo 1: Z-Transform of an FIR filter ────────────────────────── */

static void demo_fir_ztransform(void)
{
    print_separator("Demo 1: Z-Transform of FIR {1, -2, 1}");

    double b[] = {1.0, -2.0, 1.0}; /* H(z) = 1 - 2z^{-1} + z^{-2} */
    int b_len = 3;

    printf("\n  h[n] = {1, -2, 1}\n");
    printf("  H(z) = 1 - 2z^{-1} + z^{-2}\n");
    printf("       = (1 - z^{-1})^2  → double zero at z = 1\n");
    printf("\n  Evaluating H(z) at several points:\n\n");
    printf("  %8s  %12s %12s  %10s  %10s\n",
           "z", "Re{H(z)}", "Im{H(z)}", "|H(z)|", "∠H(z)°");

    double test_z[][2] = {
        {1.0, 0.0},    /* z = 1 (DC) */
        {-1.0, 0.0},   /* z = -1 (Nyquist) */
        {0.0, 1.0},    /* z = j */
        {0.0, -1.0},   /* z = -j */
        {2.0, 0.0},    /* z = 2 (outside unit circle) */
        {0.5, 0.0},    /* z = 0.5 (inside unit circle) */
    };
    int n_tests = 6;

    for (int i = 0; i < n_tests; i++) {
        Complex z = {test_z[i][0], test_z[i][1]};
        Complex Hz = eval_poly(b, b_len, z);
        printf("  %4.1f%+4.1fj  %+12.4f %+12.4f  %10.4f  %10.1f\n",
               z.re, z.im, Hz.re, Hz.im,
               complex_mag(Hz), complex_phase(Hz) * 180.0 / M_PI);
    }

    printf("\n  ✓ H(1) = 0 → this filter blocks DC (constant signals).\n");
    printf("  ✓ H(-1) = 4 → maximum gain at Nyquist (high-pass behaviour).\n");
}

/* ── Demo 2: Pole-zero map and stability ─────────────────────────── */

static void demo_pole_zero(void)
{
    print_separator("Demo 2: Pole-Zero Map & Stability");

    printf("\n  System: H(z) = 1 / (1 - 0.9 z^{-1})\n");
    printf("  This is a 1st-order IIR with one pole at z = 0.9\n\n");

    double pole_radii[] = {0.5, 0.9, 0.99, 1.0, 1.1};
    int n_poles = 5;

    printf("  %-10s  %-8s  %-10s  %-40s\n",
           "Pole at z=", "|pole|", "Stable?", "Impulse response behaviour");
    printf("  %-10s  %-8s  %-10s  %-40s\n",
           "----------", "------", "-------", "----------------------------");

    for (int i = 0; i < n_poles; i++) {
        double p = pole_radii[i];
        const char *stable = (fabs(p) < 1.0) ? "✓ YES" :
                             (fabs(p) == 1.0) ? "⚠ MARGINAL" : "✗ NO";
        const char *behaviour;
        if (fabs(p) < 1.0) behaviour = "Decays: h[n] = p^n → 0";
        else if (fabs(p) == 1.0) behaviour = "Constant: h[n] = 1 (never decays)";
        else behaviour = "Grows: h[n] = p^n → ∞ (unstable!)";

        printf("  %-10.2f  %-8.2f  %-10s  %-40s\n",
               p, fabs(p), stable, behaviour);
    }

    /* Show impulse response for p = 0.9 */
    printf("\n  Impulse response for p = 0.9 (first 16 samples):\n");
    for (int n = 0; n < 16; n++) {
        double h_n = pow(0.9, n);
        printf("    h[%2d] = %.6f  ", n, h_n);
        int bar = (int)(h_n * 30.0);
        for (int b = 0; b < bar; b++) printf("█");
        printf("\n");
    }
    printf("    → converges to 0 (stable!)\n");
}

/* ── Demo 3: Frequency response from H(z) ────────────────────────── */

static void demo_freq_response(void)
{
    print_separator("Demo 3: Frequency Response — Magnitude & Phase");

    /* Simple lowpass: H(z) = (1 + z^{-1}) / 2 */
    double b[] = {0.5, 0.5};
    int b_len = 2;

    printf("\n  H(z) = (1 + z^{-1}) / 2  (2-point average → lowpass)\n");
    printf("  Frequency response = H(e^{jω}) on the unit circle\n\n");

    printf("  %8s  %8s  %10s  %8s  %30s\n",
           "f/f_s", "ω/π", "|H|", "|H| dB", "Magnitude");
    printf("  %8s  %8s  %10s  %8s  %30s\n",
           "-----", "----", "---", "------", "---------");

    int n_points = 17;
    for (int i = 0; i <= n_points; i++) {
        double omega = M_PI * (double)i / (double)n_points;
        Complex z = complex_from_polar(1.0, omega);
        Complex Hz = eval_poly(b, b_len, z);
        double mag = complex_mag(Hz);
        double mag_db = (mag > 1e-10) ? 20.0 * log10(mag) : -100.0;

        /* ASCII bar */
        int bar = (int)(mag * 30.0);
        char bar_str[32];
        int j;
        for (j = 0; j < bar && j < 30; j++) bar_str[j] = '#';
        bar_str[j] = '\0';

        printf("  %8.3f  %8.3f  %10.4f  %8.1f  |%-30s|\n",
               (double)i / (2.0 * n_points),
               (double)i / (double)n_points,
               mag, mag_db, bar_str);
    }
    printf("\n  ✓ Magnitude drops to 0 at Nyquist (ω = π) → lowpass!\n");
}

/* ── Demo 4: Pole radius and resonance sharpness ─────────────────── */

static void demo_resonance(void)
{
    print_separator("Demo 4: Pole Radius Controls Resonance Sharpness");

    /* Complex conjugate pole pair at angle θ = π/4 (f = f_s/8) */
    double theta = M_PI / 4.0;
    double radii[] = {0.5, 0.8, 0.95, 0.99};
    int n_radii = 4;

    printf("\n  Complex pole pair at angle θ = π/4 (f = f_s/8)\n");
    printf("  Varying pole radius to show selectivity:\n");

    for (int ri = 0; ri < n_radii; ri++) {
        double r = radii[ri];

        /* H(z) = 1 / (1 - 2r·cos(θ)·z^{-1} + r²·z^{-2}) */
        double a[2];
        a[0] = -2.0 * r * cos(theta);
        a[1] = r * r;

        double b[] = {1.0};
        int b_len = 1;

        printf("\n  r = %.2f:\n", r);
        printf("  %6s  %8s  %s\n", "ω/π", "|H| dB", "");

        for (int i = 0; i <= 16; i++) {
            double omega = M_PI * (double)i / 16.0;
            Complex z = complex_from_polar(1.0, omega);
            Complex Hz = eval_h(b, b_len, a, 2, z);
            double mag = complex_mag(Hz);
            double mag_db = 20.0 * log10(mag + 1e-30);

            /* Clamp for display */
            if (mag_db < -40.0) mag_db = -40.0;
            if (mag_db > 40.0) mag_db = 40.0;

            int bar = (int)((mag_db + 40.0) / 80.0 * 40.0);
            if (bar < 0) bar = 0;
            if (bar > 40) bar = 40;

            char marker = (i == 4) ? '*' : ' '; /* θ = π/4 → i=4 */
            printf("  %6.3f  %+8.1f  ", (double)i / 16.0, mag_db);
            for (int b2 = 0; b2 < bar; b2++) printf("█");
            printf("%c\n", marker);
        }
    }
    printf("\n  * marks ω = π/4 (the resonant frequency)\n");
    printf("  ✓ As r → 1, the peak gets sharper (more selective).\n");
}

/* ── Demo 5: FIR vs IIR comparison ───────────────────────────────── */

static void demo_fir_vs_iir(void)
{
    print_separator("Demo 5: FIR (All-Zero) vs IIR (Poles + Zeros)");

    /* FIR lowpass: 5-tap average */
    double b_fir[] = {0.2, 0.2, 0.2, 0.2, 0.2};
    int b_fir_len = 5;

    /* IIR lowpass: 1st-order, cutoff ~ same region */
    double b_iir[] = {0.1};
    double a_iir[] = {-0.9}; /* pole at 0.9 */
    int b_iir_len = 1, a_iir_len = 1;

    printf("\n  FIR: H(z) = 0.2(1 + z^{-1} + z^{-2} + z^{-3} + z^{-4})\n");
    printf("  IIR: H(z) = 0.1 / (1 - 0.9 z^{-1})\n\n");

    printf("  %6s  %10s  %10s\n", "ω/π", "FIR |H| dB", "IIR |H| dB");

    for (int i = 0; i <= 16; i++) {
        double omega = M_PI * (double)i / 16.0;
        Complex z = complex_from_polar(1.0, omega);

        Complex Hz_fir = eval_poly(b_fir, b_fir_len, z);
        Complex Hz_iir = eval_h(b_iir, b_iir_len, a_iir, a_iir_len, z);

        double mag_fir = 20.0 * log10(complex_mag(Hz_fir) + 1e-30);
        double mag_iir = 20.0 * log10(complex_mag(Hz_iir) + 1e-30);

        printf("  %6.3f  %+10.1f  %+10.1f\n",
               (double)i / 16.0, mag_fir, mag_iir);
    }

    printf("\n  ✓ IIR achieves steeper rolloff with fewer coefficients.\n");
    printf("  ✓ FIR has linear phase (symmetric coefficients).\n");
}

/* ── Demo 6: Difference equation → impulse response ──────────────── */

static void demo_diff_eq(void)
{
    print_separator("Demo 6: Difference Equation → Impulse Response");

    /* y[n] = x[n] + 0.8*y[n-1] - 0.5*y[n-2]
     * B = {1}, A = {1, -0.8, 0.5} → poles from 1 - 0.8z^{-1} + 0.5z^{-2} = 0 */
    double a1 = 0.8, a2 = -0.5;

    printf("\n  Difference equation: y[n] = x[n] + 0.8·y[n-1] - 0.5·y[n-2]\n");
    printf("  Transfer function: H(z) = 1 / (1 - 0.8 z^{-1} + 0.5 z^{-2})\n");

    /* Poles: z = (0.8 ± sqrt(0.64 - 2.0)) / 2 — complex conjugate pair */
    double disc = a1 * a1 + 4.0 * a2; /* 0.64 - 2.0 = -1.36 */
    if (disc < 0) {
        double real_part = a1 / 2.0;
        double imag_part = sqrt(-disc) / 2.0;
        double pole_mag = sqrt(real_part * real_part + imag_part * imag_part);
        printf("  Poles: %.4f ± %.4fj  (|pole| = %.4f)\n",
               real_part, imag_part, pole_mag);
        printf("  Stable? %s (|pole| %s 1)\n",
               pole_mag < 1.0 ? "✓ YES" : "✗ NO",
               pole_mag < 1.0 ? "<" : "≥");
    }

    /* Compute impulse response by running the difference equation */
    int n = 32;
    double h[32];
    double y_prev1 = 0.0, y_prev2 = 0.0;

    printf("\n  Impulse response (computed by running the equation):\n");
    for (int i = 0; i < n; i++) {
        double x_i = (i == 0) ? 1.0 : 0.0;
        h[i] = x_i + a1 * y_prev1 + a2 * y_prev2;
        y_prev2 = y_prev1;
        y_prev1 = h[i];

        if (i < 20) {
            printf("    h[%2d] = %+10.6f  ", i, h[i]);
            int bar = (int)(h[i] * 15.0);
            if (bar >= 0) {
                printf("  |");
                for (int b2 = 0; b2 < bar && b2 < 30; b2++) printf("█");
            } else {
                int spaces = 15 + bar;
                if (spaces < 0) spaces = 0;
                for (int b2 = 0; b2 < spaces; b2++) printf(" ");
                for (int b2 = 0; b2 < -bar && b2 < 30; b2++) printf("█");
                printf("|");
            }
            printf("\n");
        }
    }
    printf("    ...\n");
    printf("\n  ✓ Oscillating, decaying response → stable system.\n");
}

/* ── Demo 7: Time delay property ─────────────────────────────────── */

static void demo_delay(void)
{
    print_separator("Demo 7: Time Delay Property — z^{-k}");

    double x[16];
    gen_impulse(x, 16, 3); /* impulse at n=3 */

    /* Delay by 4: multiply X(z) by z^{-4} → convolve with δ[n-4] */
    double delay_filter[5];
    int i;
    for (i = 0; i < 4; i++) delay_filter[i] = 0.0;
    delay_filter[4] = 1.0; /* δ[n-4] */

    double y[16 + 5 - 1];
    int y_len = convolve(x, 16, delay_filter, 5, y);

    printf("\n  x[n] = δ[n-3]  (impulse at n=3)\n");
    printf("  Multiply by z^{-4} → delay by 4 → impulse at n=7\n\n");

    printf("  %s  %10s  %10s\n", "n", "x[n]", "y[n] (delayed)");
    for (i = 0; i < 12 && i < y_len; i++) {
        double xi = (i < 16) ? x[i] : 0.0;
        printf("  %d  %10.1f  %10.1f\n", i, xi, y[i]);
    }
    (void)y_len;
    printf("\n  ✓ Impulse moved from n=3 to n=7 (delayed by 4).\n");
}

/* ── Main ────────────────────────────────────────────────────────── */

int main(void)
{
    printf("╔══════════════════════════════════════════════════════════╗\n");
    printf("║  Chapter 5: The Z-Transform                            ║\n");
    printf("║  DSP Tutorial Suite — Interactive Demos                 ║\n");
    printf("╚══════════════════════════════════════════════════════════╝\n");

    demo_fir_ztransform();
    demo_pole_zero();
    demo_freq_response();
    demo_resonance();
    demo_fir_vs_iir();
    demo_diff_eq();
    demo_delay();

    printf("\n════════════════════════════════════════════════════════════\n");
    printf("  All demos complete. See chapters/05-z-transform.md\n");
    printf("  for theory, derivations, and exercises.\n");
    printf("════════════════════════════════════════════════════════════\n\n");

    return 0;
}

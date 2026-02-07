/**
 * Chapter 4 Demo: LTI Systems & Discrete Convolution
 *
 * Interactive demonstrations:
 *   1. Convolution with a unit impulse (identity property)
 *   2. Convolution with a shifted impulse (delay)
 *   3. Moving-average as convolution
 *   4. Echo effect via convolution
 *   5. Commutativity: x*h == h*x
 *   6. Associativity: cascading two filters
 *   7. Cross-correlation for delay estimation
 *   8. Auto-correlation of a periodic signal
 *   9. BIBO stability check
 *  10. Signal energy and power
 *
 * Build:  make chapters && ./build/bin/ch04s
 *
 * References:
 *   - Oppenheim & Willsky, Signals and Systems, Ch 1-2
 *   - Oppenheim & Schafer, Discrete-Time Signal Processing, Ch 2
 *   - Proakis & Manolakis, Digital Signal Processing, Ch 2
 *   - Smith, Scientist & Engineer's Guide to DSP, Ch 5-7
 *   - Lyons, Understanding DSP, Ch 5
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "convolution.h"
#include "signal_gen.h"
#include "dsp_utils.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ── Helpers ─────────────────────────────────────────────────────── */

static void print_array(const char *name, const double *x, int n, int count)
{
    if (count > n) count = n;
    printf("\n── %s (%d of %d samples) ──\n", name, count, n);
    for (int i = 0; i < count; i++) {
        printf("  [%2d] = %+8.4f  ", i, x[i]);
        int bar = (int)(x[i] * 15.0);
        if (bar >= 0) {
            printf("  |");
            for (int b = 0; b < bar && b < 30; b++) printf("█");
        } else {
            int spaces = 15 + bar;
            if (spaces < 0) spaces = 0;
            for (int b = 0; b < spaces; b++) printf(" ");
            for (int b = 0; b < -bar && b < 30; b++) printf("█");
            printf("|");
        }
        printf("\n");
    }
}

static void print_separator(const char *title)
{
    printf("\n╔══════════════════════════════════════════════════════════╗\n");
    printf("║  %-56s║\n", title);
    printf("╚══════════════════════════════════════════════════════════╝\n");
}

/* ── Demo 1: Convolution with unit impulse ───────────────────────── */

static void demo_impulse_identity(void)
{
    print_separator("Demo 1: Convolution with δ[n] (Identity Property)");

    double x[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    double delta[] = {1.0};
    double y[5 + 1 - 1]; /* x_len + h_len - 1 */

    int y_len = convolve(x, 5, delta, 1, y);

    printf("\n  x[n] = {1, 2, 3, 4, 5}\n");
    printf("  h[n] = δ[n] = {1}\n");
    printf("  Output length: %d (should be 5)\n", y_len);

    print_array("x * δ  (should equal x)", y, y_len, y_len);

    printf("\n  ✓ Convolving with δ[n] gives back the original signal.\n");
}

/* ── Demo 2: Convolution with shifted impulse ────────────────────── */

static void demo_shifted_impulse(void)
{
    print_separator("Demo 2: Convolution with δ[n-2] (Delay by 2)");

    double x[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    double delta_shifted[] = {0.0, 0.0, 1.0}; /* δ[n-2] */
    double y[5 + 3 - 1];

    int y_len = convolve(x, 5, delta_shifted, 3, y);

    print_array("x[n]", x, 5, 5);
    print_array("x * δ[n-2]  (delayed by 2)", y, y_len, y_len);

    printf("\n  ✓ Signal appears shifted right by 2 samples.\n");
}

/* ── Demo 3: Moving average as convolution ───────────────────────── */

static void demo_moving_average(void)
{
    print_separator("Demo 3: 3-Point Moving Average as Convolution");

    double x[] = {0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0};
    int x_len = 8;
    double h[] = {1.0/3.0, 1.0/3.0, 1.0/3.0}; /* 3-pt average */
    int h_len = 3;
    double y[8 + 3 - 1];

    int y_len = convolve(x, x_len, h, h_len, y);

    printf("\n  h[n] = {1/3, 1/3, 1/3}  (uniform 3-point average)\n");
    print_array("Input  (rectangular pulse)", x, x_len, x_len);
    print_array("Output (smoothed)", y, y_len, y_len);

    printf("\n  ✓ Sharp edges are smoothed. Output is longer than input.\n");
}

/* ── Demo 4: Echo effect ─────────────────────────────────────────── */

static void demo_echo(void)
{
    print_separator("Demo 4: Echo Effect via Convolution");

    /* Simple click signal */
    double x[32];
    gen_impulse(x, 32, 0);

    /* Echo impulse response: direct + 0.6*delayed by 8 + 0.3*delayed by 16 */
    double h[17];
    memset(h, 0, sizeof(h));
    h[0]  = 1.0;   /* direct */
    h[8]  = 0.6;   /* first echo at 8 samples */
    h[16] = 0.3;   /* second echo at 16 samples */

    double y[32 + 17 - 1];
    int y_len = convolve(x, 32, h, 17, y);

    printf("\n  h[n] = δ[n] + 0.6·δ[n-8] + 0.3·δ[n-16]\n");
    printf("  This creates two echoes at delays of 8 and 16 samples.\n");
    print_array("Input (single impulse)", x, 32, 20);
    print_array("Output (impulse + echoes)", y, y_len, 32);
}

/* ── Demo 5: Commutativity ───────────────────────────────────────── */

static void demo_commutativity(void)
{
    print_separator("Demo 5: Commutativity — x*h == h*x");

    double x[] = {1.0, -1.0, 2.0, 0.5};
    double h[] = {0.5, 1.0, 0.5};

    double y1[4 + 3 - 1], y2[4 + 3 - 1];
    int len1 = convolve(x, 4, h, 3, y1);
    int len2 = convolve(h, 3, x, 4, y2);

    print_array("x * h", y1, len1, len1);
    print_array("h * x", y2, len2, len2);

    double max_diff = 0.0;
    for (int i = 0; i < len1; i++) {
        double d = fabs(y1[i] - y2[i]);
        if (d > max_diff) max_diff = d;
    }
    printf("\n  Max |x*h - h*x| = %.2e  → %s\n",
           max_diff, max_diff < 1e-12 ? "✓ Identical!" : "✗ Different");
}

/* ── Demo 6: Associativity (cascading) ───────────────────────────── */

static void demo_associativity(void)
{
    print_separator("Demo 6: Associativity — Cascading Filters");

    double x[] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0};
    int x_len = 8;
    double h1[] = {0.5, 0.5};             /* 2-point average */
    double h2[] = {1.0/3, 1.0/3, 1.0/3}; /* 3-point average */
    int h1_len = 2, h2_len = 3;

    /* Method A: cascade (x * h1) * h2 */
    double temp[8 + 2 - 1];
    int temp_len = convolve(x, x_len, h1, h1_len, temp);
    double yA[9 + 3 - 1];
    int yA_len = convolve(temp, temp_len, h2, h2_len, yA);

    /* Method B: combine h_total = h1 * h2, then x * h_total */
    double h_total[2 + 3 - 1];
    int ht_len = convolve(h1, h1_len, h2, h2_len, h_total);
    double yB[8 + 4 - 1];
    int yB_len = convolve(x, x_len, h_total, ht_len, yB);

    print_array("(x * h1) * h2", yA, yA_len, yA_len);
    print_array("x * (h1 * h2)", yB, yB_len, yB_len);

    double max_diff = 0.0;
    int min_len = yA_len < yB_len ? yA_len : yB_len;
    for (int i = 0; i < min_len; i++) {
        double d = fabs(yA[i] - yB[i]);
        if (d > max_diff) max_diff = d;
    }
    printf("\n  Max error = %.2e → %s\n",
           max_diff, max_diff < 1e-12 ? "✓ Identical!" : "Check");
    printf("  ✓ Two cascaded filters = one combined filter (h1 * h2).\n");
}

/* ── Demo 7: Cross-correlation for delay estimation ──────────────── */

static void demo_cross_correlation(void)
{
    print_separator("Demo 7: Cross-Correlation — Delay Estimation");

    int n = 32;
    double x[32], y[32];

    /* x = brief pulse */
    memset(x, 0, sizeof(x));
    x[4] = 1.0; x[5] = 0.8; x[6] = 0.5; x[7] = 0.2;

    /* y = same pulse delayed by 10 samples + noise */
    memset(y, 0, sizeof(y));
    y[14] = 1.0; y[15] = 0.8; y[16] = 0.5; y[17] = 0.2;

    double r[32 + 32 - 1];
    int r_len = cross_correlate(x, n, y, n, r);

    /* Find peak */
    int peak_idx = 0;
    double peak_val = r[0];
    for (int i = 1; i < r_len; i++) {
        if (r[i] > peak_val) {
            peak_val = r[i];
            peak_idx = i;
        }
    }
    int estimated_delay = peak_idx - (n - 1); /* centre is at n-1 */

    printf("\n  x: pulse at samples 4-7\n");
    printf("  y: same pulse at samples 14-17 (delay = 10)\n");
    printf("  Cross-correlation peak at index %d → lag = %d\n",
           peak_idx, estimated_delay);
    printf("  ✓ Estimated delay: %d samples %s\n", estimated_delay,
           estimated_delay == 10 ? "(correct!)" : "(check)");

    print_array("Cross-correlation (around peak)", r + peak_idx - 5, 11, 11);
}

/* ── Demo 8: Auto-correlation of a periodic signal ───────────────── */

static void demo_auto_correlation(void)
{
    print_separator("Demo 8: Auto-Correlation of a Periodic Signal");

    int n = 64;
    double x[64];
    gen_sine(x, n, 1.0, 8.0, 64.0, 0.0); /* 8 Hz at 64 Hz → period = 8 */

    double r[2 * 64 - 1];
    int r_len = auto_correlate(x, n, r);
    (void)r_len;

    /* Centre of auto-correlation is at index n-1 */
    int centre = n - 1;

    printf("\n  Signal: 8 Hz sine at 64 Hz sample rate → period = 8 samples\n");
    printf("  Auto-correlation shows peaks at multiples of 8:\n\n");

    /* Show lags 0 to 32 */
    for (int lag = 0; lag <= 32; lag++) {
        double val = r[centre + lag];
        char marker = ' ';
        if (lag % 8 == 0) marker = '*';
        printf("  lag %2d: %+8.4f %c\n", lag, val, marker);
    }
    printf("\n  * = multiples of 8 (expected peaks for period-8 signal)\n");
}

/* ── Demo 9: BIBO stability check ────────────────────────────────── */

static void demo_bibo_stability(void)
{
    print_separator("Demo 9: BIBO Stability Check");

    double h_fir[] = {0.25, 0.5, 0.25}; /* stable FIR */

    /* Simulated "IIR" — decaying exponential impulse response */
    double h_iir[64];
    gen_exponential(h_iir, 64, 1.0, 0.9); /* |base| < 1 → decays */

    /* Unstable: growing exponential */
    double h_unstable[64];
    gen_exponential(h_unstable, 64, 1.0, 1.1); /* |base| > 1 → grows */

    printf("\n  FIR {0.25, 0.5, 0.25}:  BIBO stable? %s\n",
           is_bibo_stable(h_fir, 3) ? "✓ YES" : "✗ NO");
    printf("  Decaying exp (base=0.9):  BIBO stable? %s\n",
           is_bibo_stable(h_iir, 64) ? "✓ YES" : "✗ NO");
    printf("  Growing exp (base=1.1):   BIBO stable? %s\n",
           is_bibo_stable(h_unstable, 64) ? "✓ YES (finite truncation)" : "✗ NO");
    printf("\n  Note: The growing exponential is \"stable\" only because we truncated"
           "\n  it to 64 samples. A true IIR with |pole| > 1 would be unstable.\n");
}

/* ── Demo 10: Signal energy and power ────────────────────────────── */

static void demo_energy_power(void)
{
    print_separator("Demo 10: Signal Energy & Power");

    int n = 64;
    double sine[64], noise[64], impulse[64];

    gen_sine(sine, n, 1.0, 4.0, 64.0, 0.0);
    gen_white_noise(noise, n, 1.0, 42);
    gen_impulse(impulse, n, 0);

    printf("\n  %-20s  Energy     Power      RMS\n", "Signal");
    printf("  %-20s  --------   --------   --------\n", "------");

    printf("  %-20s  %8.4f   %8.4f   %8.4f\n", "Sine (A=1)",
           signal_energy(sine, n), signal_power(sine, n), rms(sine, n));
    printf("  %-20s  %8.4f   %8.4f   %8.4f\n", "White noise",
           signal_energy(noise, n), signal_power(noise, n), rms(noise, n));
    printf("  %-20s  %8.4f   %8.4f   %8.4f\n", "Impulse",
           signal_energy(impulse, n), signal_power(impulse, n), rms(impulse, n));

    printf("\n  Sine energy = N/2 = %d/2 = %.1f  (theoretical for A=1)\n",
           n, (double)n / 2.0);
    printf("  Sine power  = 1/2 = 0.5000  (A²/2 for A=1)\n");
}

/* ── Main ────────────────────────────────────────────────────────── */

int main(void)
{
    printf("╔══════════════════════════════════════════════════════════╗\n");
    printf("║  Chapter 4: LTI Systems & Discrete Convolution         ║\n");
    printf("║  DSP Tutorial Suite — Interactive Demos                 ║\n");
    printf("╚══════════════════════════════════════════════════════════╝\n");

    demo_impulse_identity();
    demo_shifted_impulse();
    demo_moving_average();
    demo_echo();
    demo_commutativity();
    demo_associativity();
    demo_cross_correlation();
    demo_auto_correlation();
    demo_bibo_stability();
    demo_energy_power();

    printf("\n════════════════════════════════════════════════════════════\n");
    printf("  All demos complete. See chapters/04-lti-systems.md\n");
    printf("  for theory, derivations, and exercises.\n");
    printf("════════════════════════════════════════════════════════════\n\n");

    return 0;
}

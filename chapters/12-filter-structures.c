/**
 * @file 12-filter-structures.c
 * @brief Chapter 12 Demo: Filter Structures — DF1, DF2T, Cascaded Biquads
 *
 * Interactive demonstrations:
 *   1. Direct Form I vs Direct Form II Transposed — same biquad, two forms
 *   2. Numerical precision comparison — DF1 vs DF2T under extreme conditions
 *   3. Cascaded biquads (SOS) vs single high-order section
 *   4. Per-sample vs block processing
 *   5. State inspection: watching the delay lines in real time
 *   6. Coefficient sensitivity: small perturbation → big output change
 *
 * Build:  make release && ./build/bin/ch12
 *
 * References:
 *   - Oppenheim & Schafer, Discrete-Time Signal Processing, Ch 6
 *   - Proakis & Manolakis, Digital Signal Processing, Ch 9
 *   - Lyons, Understanding DSP, Ch 6-7
 *
 * Prerequisites: Ch 06 (Frequency Response), Ch 11 (IIR Filter Design)
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dsp_utils.h"
#include "iir.h"
#include "signal_gen.h"

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

/* ── Demo 1: DF1 vs DF2T — same biquad, two implementations ─────── */

static void demo_df1_vs_df2t(void)
{
    print_separator("Demo 1: Direct Form I vs Direct Form II Transposed");

    /*
     *  Direct Form I (4 state variables):
     *
     *   x[n] ──►[b0]──┬──► (+) ──────────────────► y[n]
     *               z⁻¹│       ▲           ▲
     *   x[n-1]──►[b1]──┘       │           │
     *               z⁻¹│   [-a1]◄──y[n-1]  │
     *   x[n-2]──►[b2]──┘       │           │
     *                       [-a2]◄──y[n-2]──┘
     *
     *
     *  Direct Form II Transposed (2 state variables):
     *
     *   x[n]──►[b0]──►(+)──────────────────────► y[n]
     *      │           ▲                    │
     *      │          s1                 [-a1]
     *      │           ▲                    │
     *      ├──►[b1]──►(+)◄─────────────────┘
     *      │           ▲                    │
     *      │          s2                 [-a2]
     *      │           ▲                    │
     *      └──►[b2]──►(+)◄─────────────────┘
     */

    /* Example: 2nd-order resonant filter (pole at angle π/4, radius 0.9) */
    Biquad bq;
    double theta = M_PI / 4.0;
    double r = 0.9;
    bq.b0 = 1.0;
    bq.b1 = 0.0;
    bq.b2 = 0.0;
    bq.a1 = -2.0 * r * cos(theta);
    bq.a2 = r * r;

    printf("\n  Biquad: resonance at ω = π/4, pole radius = 0.9\n");
    printf("  b = {%.3f, %.3f, %.3f}, a = {1, %.6f, %.6f}\n\n",
           bq.b0, bq.b1, bq.b2, bq.a1, bq.a2);

    /* Process an impulse through both forms */
    BiquadDF1State  df1_state;
    BiquadDF2TState df2t_state;
    biquad_df1_init(&df1_state);
    biquad_df2t_init(&df2t_state);

    printf("  %4s  %12s  %12s  %12s\n",
           "n", "DF1 y[n]", "DF2T y[n]", "Difference");

    double max_diff = 0.0;
    for (int i = 0; i < 20; i++) {
        double x = (i == 0) ? 1.0 : 0.0;

        double y_df1  = biquad_process_df1(&bq, &df1_state, x);
        double y_df2t = biquad_process_df2t(&bq, &df2t_state, x);

        double diff = fabs(y_df1 - y_df2t);
        if (diff > max_diff) max_diff = diff;

        printf("  %4d  %+12.8f  %+12.8f  %12.2e\n",
               i, y_df1, y_df2t, diff);
    }

    printf("\n  Max difference: %.2e\n", max_diff);
    printf("  ✓ Both forms produce identical results (within floating-point).\n");
    printf("  ✓ DF1 uses 4 state variables; DF2T uses only 2.\n");
}

/* ── Demo 2: Numerical precision under extreme conditions ────────── */

static void demo_precision(void)
{
    print_separator("Demo 2: Numerical Precision — High-Q Narrow-Band");

    /*
     * A narrow-band filter has poles very close to the unit circle.
     * This stresses numerical precision because small coefficient
     * errors produce large response changes.
     *
     * High Q = poles close to unit circle = narrow bandwidth
     *
     * Pole radius 0.999 → 3-dB bandwidth ≈ 0.001·π ≈ 0.05% of Nyquist
     */

    double r = 0.999;
    double theta = M_PI / 4.0;

    Biquad bq;
    bq.b0 = (1.0 - r) * (1.0 - r);  /* Unity peak gain (approx) */
    bq.b1 = 0.0;
    bq.b2 = 0.0;
    bq.a1 = -2.0 * r * cos(theta);
    bq.a2 = r * r;

    printf("\n  Narrow-band resonator: r = %.3f, θ = π/4\n", r);
    printf("  Pole distance from unit circle: %.4f\n", 1.0 - r);
    printf("  Approximate 3-dB bandwidth: %.5f π rad/sample\n\n",
           2.0 * (1.0 - r));

    /* Process a burst of sine at the resonant frequency */
    BiquadDF1State  df1;
    BiquadDF2TState df2t;
    biquad_df1_init(&df1);
    biquad_df2t_init(&df2t);

    double max_diff = 0.0;
    int n = 200;

    printf("  Feeding resonant sine for %d samples, then silence...\n\n", n / 2);
    printf("  %4s  %14s  %14s  %12s\n", "n", "DF1", "DF2T", "Difference");

    for (int i = 0; i < n; i++) {
        double x = (i < n / 2) ? sin(theta * i) : 0.0;
        double y1 = biquad_process_df1(&bq, &df1, x);
        double y2 = biquad_process_df2t(&bq, &df2t, x);
        double diff = fabs(y1 - y2);
        if (diff > max_diff) max_diff = diff;

        if (i < 10 || (i >= n / 2 - 2 && i < n / 2 + 5) || i >= n - 5) {
            printf("  %4d  %+14.10f  %+14.10f  %12.2e\n", i, y1, y2, diff);
        } else if (i == 10) {
            printf("  ...\n");
        }
    }

    printf("\n  Max DF1 vs DF2T difference: %.2e\n", max_diff);
    printf("  ✓ For double-precision float, both forms agree well.\n");
    printf("  ✓ In fixed-point (16-bit), DF2T is generally more robust.\n");
}

/* ── Demo 3: Cascaded SOS vs single high-order section ───────────── */

static void demo_sos_vs_direct(void)
{
    print_separator("Demo 3: Cascaded SOS vs Single High-Order Section");

    /*
     *  6th-order IIR: can be implemented as:
     *
     *  (A) Single section with 6 feedforward + 6 feedback coefficients
     *      → Large intermediate values, coefficient sensitivity
     *
     *  (B) Cascade of 3 biquad sections
     *      → Each section handles only 2 poles/zeros
     *      → Much more numerically stable
     *
     *       ┌─────┐   ┌─────┐   ┌─────┐
     *  x ──►│SOS 0│──►│SOS 1│──►│SOS 2│──► y
     *       └─────┘   └─────┘   └─────┘
     */

    /* Design a 6th-order Butterworth */
    SOSCascade sos;
    butterworth_lowpass(6, 0.15, &sos);

    printf("\n  6th-order Butterworth lowpass at fc = 0.15\n\n");

    /* Show the SOS cascade coefficients */
    printf("  Cascaded SOS representation (%d sections):\n\n", sos.n_sections);
    printf("  %-4s  %10s %10s %10s | %10s %10s\n",
           "Sec", "b0", "b1", "b2", "a1", "a2");
    for (int i = 0; i < sos.n_sections; i++) {
        const Biquad *bq = &sos.sections[i];
        printf("  [%d]   %10.6f %10.6f %10.6f | %10.6f %10.6f\n",
               i, bq->b0, bq->b1, bq->b2, bq->a1, bq->a2);
    }

    /* Expand to a single high-order transfer function for comparison */
    /* B(z) = product of all numerators, A(z) = product of all denominators */
    printf("\n  Equivalent single-section coefficients (expanded):\n");

    /* Polynomial multiplication: start with [1.0], multiply each section */
    double b_expanded[16] = {0.0};
    double a_expanded[16] = {0.0};
    b_expanded[0] = 1.0;
    a_expanded[0] = 1.0;
    int b_len = 1, a_len = 1;

    for (int s = 0; s < sos.n_sections; s++) {
        const Biquad *bq = &sos.sections[s];
        double bk[3] = {bq->b0, bq->b1, bq->b2};
        double ak[3] = {1.0, bq->a1, bq->a2};
        int bk_len = (bq->b2 != 0.0) ? 3 : 2;
        int ak_len = (bq->a2 != 0.0) ? 3 : 2;

        /* Multiply b_expanded by bk */
        double new_b[16] = {0.0};
        for (int i = 0; i < b_len; i++)
            for (int j = 0; j < bk_len; j++)
                new_b[i + j] += b_expanded[i] * bk[j];
        b_len = b_len + bk_len - 1;
        for (int i = 0; i < 16; i++) b_expanded[i] = new_b[i];

        /* Multiply a_expanded by ak */
        double new_a[16] = {0.0};
        for (int i = 0; i < a_len; i++)
            for (int j = 0; j < ak_len; j++)
                new_a[i + j] += a_expanded[i] * ak[j];
        a_len = a_len + ak_len - 1;
        for (int i = 0; i < 16; i++) a_expanded[i] = new_a[i];
    }

    printf("  B(z) = ");
    for (int i = 0; i < b_len; i++)
        printf("%+.6f z^{-%d} ", b_expanded[i], i);
    printf("\n  A(z) = ");
    for (int i = 0; i < a_len; i++)
        printf("%+.6f z^{-%d} ", a_expanded[i], i);
    printf("\n");

    /* Process an impulse through both and compare */
    printf("\n  Impulse response comparison:\n");
    printf("  %4s  %14s  %14s  %12s\n",
           "n", "Cascaded SOS", "Single section", "Difference");

    /* Reset SOS states */
    for (int i = 0; i < sos.n_sections; i++)
        biquad_df1_init(&sos.states[i]);

    double h_single[32];
    iir_filter(b_expanded, b_len, a_expanded, a_len,
               (double[]){1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
               h_single, 32);

    double max_diff = 0.0;
    for (int i = 0; i < 32; i++) {
        double x = (i == 0) ? 1.0 : 0.0;
        double y_sos = sos_process_sample(&sos, x);
        double diff = fabs(y_sos - h_single[i]);
        if (diff > max_diff) max_diff = diff;

        if (i < 12 || i >= 28) {
            printf("  %4d  %+14.10f  %+14.10f  %12.2e\n",
                   i, y_sos, h_single[i], diff);
        } else if (i == 12) {
            printf("  ...\n");
        }
    }

    printf("\n  Max difference: %.2e (due to polynomial expansion rounding)\n",
           max_diff);
    printf("  ✓ SOS cascade avoids expanding to high-order polynomials.\n");
    printf("  ✓ In practice, use SOS — never expand to single section.\n");
}

/* ── Demo 4: Per-sample vs block processing ──────────────────────── */

static void demo_per_sample_vs_block(void)
{
    print_separator("Demo 4: Per-Sample vs Block Processing");

    int n = 64;
    double input[64], out_sample[64], out_block[64];

    /* Generate test signal */
    gen_sine(input, n, 0.1, 1.0, 0.0, 1.0);

    /* Design filter */
    SOSCascade sos1, sos2;
    butterworth_lowpass(4, 0.2, &sos1);
    butterworth_lowpass(4, 0.2, &sos2);

    /* Method 1: Per-sample in a loop */
    for (int i = 0; i < n; i++)
        out_sample[i] = sos_process_sample(&sos1, input[i]);

    /* Method 2: Block processing */
    sos_process_block(&sos2, input, out_block, n);

    /* Compare */
    printf("\n  Both methods should produce identical results.\n\n");
    printf("  %4s  %12s  %12s  %12s\n",
           "n", "Per-sample", "Block", "Difference");

    double max_d = 0.0;
    for (int i = 0; i < n; i++) {
        double d = fabs(out_sample[i] - out_block[i]);
        if (d > max_d) max_d = d;
        if (i < 8 || i >= 60) {
            printf("  %4d  %+12.8f  %+12.8f  %12.2e\n",
                   i, out_sample[i], out_block[i], d);
        } else if (i == 8) {
            printf("  ...\n");
        }
    }

    printf("\n  Max diff: %.2e\n", max_d);
    printf("  ✓ Results are bit-identical — same computation underneath.\n");
    printf("  ✓ Per-sample: suitable for real-time (one sample at a time).\n");
    printf("  ✓ Block: better cache performance for offline processing.\n");
}

/* ── Demo 5: Watching the delay lines ────────────────────────────── */

static void demo_state_inspection(void)
{
    print_separator("Demo 5: State Inspection — Delay Lines in Action");

    /*
     * DF1 state: x[n-1], x[n-2], y[n-1], y[n-2]
     *
     * Watch how the state evolves as samples flow through:
     *
     *  n=0: x₁=0    x₂=0    y₁=0    y₂=0
     *  n=1: x₁=1.0  x₂=0    y₁=0.5  y₂=0
     *  n=2: x₁=0    x₂=1.0  y₁=0.2  y₂=0.5
     *  ...
     */

    Biquad bq;
    bq.b0 = 0.5;
    bq.b1 = 0.0;
    bq.b2 = 0.0;
    bq.a1 = -0.5;
    bq.a2 = 0.0;

    BiquadDF1State state;
    biquad_df1_init(&state);

    printf("\n  Biquad: H(z) = 0.5 / (1 − 0.5z⁻¹)  (1st-order IIR)\n\n");
    printf("  Input: impulse (1, 0, 0, 0, ...)\n\n");
    printf("  %4s  %6s  %8s  %10s %10s %10s %10s\n",
           "n", "x[n]", "y[n]", "state.x1", "state.x2", "state.y1", "state.y2");

    for (int i = 0; i < 12; i++) {
        double x = (i == 0) ? 1.0 : 0.0;
        double y = biquad_process_df1(&bq, &state, x);

        printf("  %4d  %6.1f  %8.6f  %10.6f %10.6f %10.6f %10.6f\n",
               i, x, y, state.x1, state.x2, state.y1, state.y2);
    }

    printf("\n  ✓ y[n] = 0.5^{n+1} — exponential decay (stable, pole at 0.5).\n");
    printf("  ✓ State stores the memory needed for the feedback loop.\n");
    printf("  ✓ Without state, the filter has no memory → output = b0·x[n].\n");
}

/* ── Demo 6: Coefficient sensitivity ─────────────────────────────── */

static void demo_sensitivity(void)
{
    print_separator("Demo 6: Coefficient Sensitivity");

    /*
     * Small changes in IIR coefficients can cause large changes in
     * the frequency response — especially for high-order or narrow-band
     * filters. This is why SOS cascade is preferred over single-section.
     *
     * We perturb one coefficient by 0.1% and measure the change.
     */

    /* 4th-order Butterworth at cutoff 0.15 */
    SOSCascade sos_orig;
    butterworth_lowpass(4, 0.15, &sos_orig);

    /* Perturbed version: modify a2 of first section by 0.1% */
    SOSCascade sos_pert = sos_orig;
    sos_pert.sections[0].a2 *= 1.001;

    printf("\n  Original a2[0] = %.10f\n", sos_orig.sections[0].a2);
    printf("  Perturbed a2[0] = %.10f  (0.1%% change)\n\n",
           sos_pert.sections[0].a2);

    /* Compare frequency responses */
    int np = 33;
    double mag_orig[33], mag_pert[33], ph_orig[33], ph_pert[33];

    sos_freq_response(&sos_orig, mag_orig, ph_orig, np);
    sos_freq_response(&sos_pert, mag_pert, ph_pert, np);

    printf("  %6s  %10s  %10s  %10s\n",
           "f/fs", "Orig dB", "Pert dB", "Change dB");

    double max_change = 0.0;
    for (int i = 0; i < np; i++) {
        double db_orig = (mag_orig[i] > 1e-10) ? 20.0 * log10(mag_orig[i]) : -80.0;
        double db_pert = (mag_pert[i] > 1e-10) ? 20.0 * log10(mag_pert[i]) : -80.0;
        double change = db_pert - db_orig;
        if (fabs(change) > max_change) max_change = fabs(change);

        printf("  %6.3f  %+10.2f  %+10.2f  %+10.4f\n",
               (double)i / (2.0 * (np - 1)), db_orig, db_pert, change);
    }

    printf("\n  Maximum response change: %.4f dB from 0.1%% coefficient change\n",
           max_change);
    printf("  ✓ SOS limits sensitivity to 2nd-order sections.\n");
    printf("  ✓ A single 4th-order section would be even more sensitive.\n");
    printf("  ✓ For fixed-point: use SOS + DF2T for best robustness.\n");
}

/* ── Main ────────────────────────────────────────────────────────── */

int main(void)
{
    printf("╔══════════════════════════════════════════════════════════╗\n");
    printf("║  Chapter 12: Filter Structures                         ║\n");
    printf("║  DSP Tutorial Suite — Interactive Demos                 ║\n");
    printf("╚══════════════════════════════════════════════════════════╝\n");

    demo_df1_vs_df2t();
    demo_precision();
    demo_sos_vs_direct();
    demo_per_sample_vs_block();
    demo_state_inspection();
    demo_sensitivity();

    printf("\n════════════════════════════════════════════════════════════\n");
    printf("  All demos complete. See chapters/12-filter-structures.md\n");
    printf("  for theory, diagrams, and exercises.\n");
    printf("════════════════════════════════════════════════════════════\n\n");

    return 0;
}

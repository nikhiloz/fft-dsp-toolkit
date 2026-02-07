/**
 * @file 11-iir-filter-design.c
 * @brief Chapter 11 Demo: IIR Filter Design
 *
 * Interactive demonstrations:
 *   1. Butterworth lowpass — maximally flat magnitude
 *   2. Butterworth highpass — complement of lowpass
 *   3. Chebyshev Type I lowpass — equiripple passband, steeper rolloff
 *   4. Order comparison: how order affects transition width
 *   5. Bilinear transform frequency warping visualisation
 *   6. IIR vs FIR design comparison (same spec, different approaches)
 *   7. Filter the signal: apply designed IIR to a noisy waveform
 *
 * Build:  make release && ./build/bin/ch11
 *
 * References:
 *   - Oppenheim & Schafer, Discrete-Time Signal Processing, Ch 8
 *   - Proakis & Manolakis, Digital Signal Processing, Ch 8
 *   - Lyons, Understanding DSP, Ch 7
 *   - Smith, Scientist & Engineer's Guide to DSP, Ch 19-21
 *
 * Prerequisites: Ch 05 (Z-Transform), Ch 06 (Frequency Response),
 *                Ch 10 (FIR Filters for contrast)
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dsp_utils.h"
#include "fft.h"
#include "iir.h"
#include "filter.h"
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

/**
 * Print SOS cascade coefficients in a readable table.
 */
static void print_sos(const SOSCascade *sos)
{
    printf("  Sections: %d,  Gain: %.6f\n\n", sos->n_sections, sos->gain);
    printf("  %-4s  %10s %10s %10s  %10s %10s\n",
           "Sec", "b0", "b1", "b2", "a1", "a2");
    for (int i = 0; i < sos->n_sections; i++) {
        const Biquad *bq = &sos->sections[i];
        printf("  [%d]   %10.6f %10.6f %10.6f  %10.6f %10.6f\n",
               i, bq->b0, bq->b1, bq->b2, bq->a1, bq->a2);
    }
}

/**
 * Plot SOS cascade magnitude response as ASCII.
 */
static void plot_sos_response(const SOSCascade *sos, int n_points)
{
    double mag[256], phase[256];
    if (n_points > 256) n_points = 256;

    sos_freq_response(sos, mag, phase, n_points);

    printf("\n  %6s  %8s  %s\n", "f/fs", "|H| dB", "");
    for (int i = 0; i < n_points; i++) {
        double mag_db = (mag[i] > 1e-10) ? 20.0 * log10(mag[i]) : -80.0;
        if (mag_db < -60.0) mag_db = -60.0;
        if (mag_db > 5.0) mag_db = 5.0;

        int bar = (int)((mag_db + 60.0) / 65.0 * 35.0);
        if (bar < 0) bar = 0;
        if (bar > 35) bar = 35;

        printf("  %6.3f  %+8.1f  ",
               (double)i / (2.0 * (n_points - 1)),
               mag_db);
        for (int j = 0; j < bar; j++) printf("█");
        printf("\n");
    }
}

/* ── Demo 1: Butterworth Lowpass ─────────────────────────────────── */

static void demo_butterworth_lp(void)
{
    print_separator("Demo 1: Butterworth Lowpass — Maximally Flat");

    /*
     *  Butterworth magnitude response:
     *
     *    |H(jΩ)|² = 1 / (1 + (Ω/Ωc)^{2N})
     *
     *              N=2      N=4      N=8
     *    0 dB  ───────  ───────  ───────
     *              \         \        │      passband (flat)
     *    -3 dB  ....\.........\.......│....  ← cutoff
     *                \         \      │
     *   -40 dB        \         \     │      stopband
     *                   \        \    │
     *
     *    Higher order → steeper transition, but more computation.
     */

    int orders[] = {2, 4, 6};
    double cutoff = 0.2; /* fc = 0.2·fs → passband 0 to 0.1·fs */

    printf("\n  Butterworth lowpass at cutoff = %.2f (normalized)\n", cutoff);
    printf("  Comparing orders 2, 4, 6:\n");

    for (int oi = 0; oi < 3; oi++) {
        int N = orders[oi];
        SOSCascade sos;
        butterworth_lowpass(N, cutoff, &sos);

        printf("\n  ─── Order %d (%d biquad sections) ───\n", N, sos.n_sections);
        print_sos(&sos);

        /* Compute -3dB point verification */
        double mag_at_cutoff[1], phase_at_cutoff[1];
        sos_freq_response(&sos, mag_at_cutoff, phase_at_cutoff, 1);

        /* Actually compute at the cutoff frequency specifically */
        double test_mag[33];
        double test_phase[33];
        sos_freq_response(&sos, test_mag, test_phase, 33);

        /* Index for cutoff: cutoff = 0.2, ranges 0..0.5 in 33 pts */
        int cutoff_idx = (int)(cutoff * 2.0 * 32.0 + 0.5);
        if (cutoff_idx >= 33) cutoff_idx = 32;
        double mag_db_at_cut = 20.0 * log10(test_mag[cutoff_idx] + 1e-30);

        printf("\n  |H| at cutoff (f=%.2f): %.4f (%.1f dB) — target: −3 dB\n",
               cutoff, test_mag[cutoff_idx], mag_db_at_cut);

        plot_sos_response(&sos, 33);
    }

    printf("\n  ✓ Higher order → steeper rolloff (sharper transition).\n");
    printf("  ✓ Butterworth: maximally flat passband (no ripple).\n");
    printf("  ✓ −3 dB point matches the design cutoff frequency.\n");
}

/* ── Demo 2: Butterworth Highpass ────────────────────────────────── */

static void demo_butterworth_hp(void)
{
    print_separator("Demo 2: Butterworth Highpass");

    double cutoff = 0.3;
    SOSCascade sos;
    butterworth_highpass(4, cutoff, &sos);

    printf("\n  4th-order Butterworth highpass, cutoff = %.2f\n", cutoff);
    print_sos(&sos);
    plot_sos_response(&sos, 33);

    printf("\n  ✓ Blocks frequencies below cutoff, passes above.\n");
    printf("  ✓ Numerator zeros at z = +1 (DC null) instead of z = −1.\n");
}

/* ── Demo 3: Chebyshev Type I ────────────────────────────────────── */

static void demo_chebyshev1(void)
{
    print_separator("Demo 3: Chebyshev Type I — Steeper Rolloff with Ripple");

    /*
     *  Chebyshev vs Butterworth (same order):
     *
     *    0 dB  ─┬─ripple─┬─    ─────────────
     *           │        │         \
     *    -R dB  └────────┘          \          Chebyshev (left)
     *                    │\          \         Butterworth (right)
     *   -40 dB           │  \         \
     *                    │    \        \
     *
     *  Steeper rolloff, but at the cost of passband ripple.
     */

    double cutoff = 0.2;
    double ripples[] = {0.5, 1.0, 3.0};
    int n_ripples = 3;

    printf("\n  4th-order Chebyshev Type I lowpass, cutoff = %.2f\n", cutoff);
    printf("  Comparing different ripple levels:\n");

    for (int ri = 0; ri < n_ripples; ri++) {
        SOSCascade sos;
        chebyshev1_lowpass(4, ripples[ri], cutoff, &sos);

        printf("\n  ─── Ripple: %.1f dB ───\n", ripples[ri]);
        print_sos(&sos);
        plot_sos_response(&sos, 33);
    }

    printf("\n  ✓ More ripple → steeper transition band.\n");
    printf("  ✓ Chebyshev poles lie on an ELLIPSE (not circle).\n");
    printf("  ✓ Trade-off: passband flatness vs transition sharpness.\n");
}

/* ── Demo 4: Order vs transition width ───────────────────────────── */

static void demo_order_comparison(void)
{
    print_separator("Demo 4: Filter Order vs Transition Width");

    double cutoff = 0.2;
    int orders[] = {1, 2, 4, 8};
    int n_orders = 4;

    printf("\n  Butterworth lowpass at cutoff = %.2f\n", cutoff);
    printf("  How quickly does the stopband attenuation increase?\n\n");

    printf("  %5s", "ω/π");
    for (int oi = 0; oi < n_orders; oi++)
        printf("  N=%-4d", orders[oi]);
    printf("\n");

    for (int i = 0; i <= 16; i++) {
        double f_norm = (double)i / 32.0;  /* 0 to 0.5 */
        printf("  %5.3f", (double)i / 16.0);

        for (int oi = 0; oi < n_orders; oi++) {
            SOSCascade sos;
            butterworth_lowpass(orders[oi], cutoff, &sos);

            double mag[17], phase[17];
            sos_freq_response(&sos, mag, phase, 17);

            double mag_db = (mag[i] > 1e-10) ? 20.0 * log10(mag[i]) : -80.0;
            printf("  %+7.1f", mag_db);
        }
        printf("\n");
        (void)f_norm;
    }

    printf("\n  ✓ Each doubling of order roughly doubles the dB/octave rolloff.\n");
    printf("  ✓ N=1: 20 dB/decade, N=2: 40 dB/decade, N=4: 80 dB/decade.\n");
}

/* ── Demo 5: Bilinear transform frequency warping ────────────────── */

static void demo_warping(void)
{
    print_separator("Demo 5: Bilinear Transform — Frequency Warping");

    /*
     * The bilinear transform maps the entire analog frequency axis
     * (−∞ to +∞) onto the digital range (−π to +π). This compression
     * introduces non-linear "warping":
     *
     *   Ωa = 2·tan(ωd / 2)
     *
     *   ωd = 0       → Ωa = 0         (DC preserved)
     *   ωd = π/4     → Ωa ≈ 0.828     (slight compression)
     *   ωd = π/2     → Ωa = 2         (noticeable warping)
     *   ωd = π       → Ωa = ∞         (Nyquist → ∞!)
     *
     * We PRE-WARP the cutoff frequency to compensate:
     *   Design the analog filter at Ωa = 2·tan(ωd/2)
     *   so the −3 dB point lands at exactly ωd after transformation.
     */

    printf("\n  Bilinear transform: s = 2(z−1)/(z+1)\n");
    printf("  Frequency warping: Ωa = 2·tan(ωd/2)\n\n");

    printf("  %12s  %12s  %12s  %s\n",
           "ωd (digital)", "ωd/π", "Ωa (analog)", "Warping");

    double test_wd[] = {0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 0.95, 1.0};
    int n_test = 9;

    for (int i = 0; i < n_test; i++) {
        double wd_pi = test_wd[i];
        double wd = wd_pi * M_PI;
        double Wa = 2.0 * tan(wd / 2.0);
        double ratio = (wd > 1e-10) ? Wa / wd : 1.0;

        printf("  %12.4f  %12.3f  %12.4f  ×%.3f\n",
               wd, wd_pi, Wa, ratio);
    }

    printf("\n  ✓ Near DC, warping is negligible (ratio ≈ 1.0).\n");
    printf("  ✓ Near Nyquist, warping is severe (ratio → ∞).\n");
    printf("  ✓ Pre-warping the cutoff ensures the −3 dB point is correct.\n");
}

/* ── Demo 6: IIR vs FIR comparison ───────────────────────────────── */

static void demo_iir_vs_fir(void)
{
    print_separator("Demo 6: IIR vs FIR — Same Spec, Different Trade-offs");

    double cutoff = 0.15;  /* lowpass cutoff */

    /* IIR: 4th-order Butterworth */
    SOSCascade sos;
    butterworth_lowpass(4, cutoff, &sos);

    /* FIR: windowed sinc with enough taps for similar transition width */
    int fir_taps = 31;
    double h_fir[31];
    fir_lowpass(h_fir, fir_taps, cutoff);

    printf("\n  Target: lowpass at cutoff = %.2f\n", cutoff);
    printf("  IIR: 4th-order Butterworth (9 multiplies/sample)\n");
    printf("  FIR: 31-tap windowed sinc  (31 multiplies/sample)\n\n");

    /* Compare at several frequencies */
    printf("  %6s  %10s  %10s  %10s\n",
           "f/fs", "IIR dB", "FIR dB", "Winner");

    double iir_mag[33], iir_phase[33];
    sos_freq_response(&sos, iir_mag, iir_phase, 33);

    /* Compute FIR frequency response manually */
    for (int i = 0; i <= 16; i++) {
        double omega = M_PI * (double)i / 16.0;

        /* IIR magnitude (already computed) */
        int idx = i * 2;
        if (idx >= 33) idx = 32;
        double iir_db = (iir_mag[idx] > 1e-10) ?
                        20.0 * log10(iir_mag[idx]) : -80.0;

        /* FIR magnitude: evaluate H(e^{jω}) = Σ h[k] e^{-jkω} */
        double fir_re = 0.0, fir_im = 0.0;
        for (int k = 0; k < fir_taps; k++) {
            fir_re += h_fir[k] * cos(k * omega);
            fir_im -= h_fir[k] * sin(k * omega);
        }
        double fir_mag = sqrt(fir_re * fir_re + fir_im * fir_im);
        double fir_db = (fir_mag > 1e-10) ? 20.0 * log10(fir_mag) : -80.0;

        const char *winner;
        if (i == 0 || (double)i / 16.0 < cutoff * 2.0 - 0.05)
            winner = "≈ tied";
        else if (fabs(iir_db) > fabs(fir_db) + 5.0)
            winner = "IIR (steeper)";
        else if (fabs(fir_db) > fabs(iir_db) + 5.0)
            winner = "FIR (steeper)";
        else
            winner = "≈ tied";

        printf("  %6.3f  %+10.1f  %+10.1f  %s\n",
               (double)i / 32.0, iir_db, fir_db, winner);
    }

    printf("\n  IIR advantages:\n");
    printf("    ✓ Steeper rolloff with fewer multiplies per sample.\n");
    printf("    ✓ Lower latency (smaller group delay in passband).\n");
    printf("  FIR advantages:\n");
    printf("    ✓ Linear phase (constant group delay).\n");
    printf("    ✓ Always stable.\n");
    printf("    ✓ Easier to design for arbitrary specifications.\n");
}

/* ── Demo 7: Filter a noisy signal ───────────────────────────────── */

static void demo_filter_signal(void)
{
    print_separator("Demo 7: IIR Filter Applied to Noisy Signal");

    int n = 128;
    double signal[128], noise[128], noisy[128], filtered[128];

    /* Generate a 0.05 fs sine wave + noise */
    gen_sine(signal, n, 0.05, 1.0, 0.0, 1.0);
    gen_white_noise(noise, n, 0.3, 42);

    for (int i = 0; i < n; i++)
        noisy[i] = signal[i] + noise[i];

    /* Design 4th-order Butterworth lowpass at 0.1 fs */
    SOSCascade sos;
    butterworth_lowpass(4, 0.1, &sos);

    /* Apply filter */
    sos_process_block(&sos, noisy, filtered, n);

    printf("\n  Signal: sine at 0.05·fs + white noise (σ = 0.3)\n");
    printf("  Filter: 4th-order Butterworth lowpass at 0.1·fs\n\n");

    /* Compute RMS of noise before/after */
    double noise_rms_before = 0.0, noise_rms_after = 0.0;
    /* Skip first 16 samples (filter transient) */
    for (int i = 16; i < n; i++) {
        double err_before = noisy[i] - signal[i];
        double err_after  = filtered[i] - signal[i];
        noise_rms_before += err_before * err_before;
        noise_rms_after  += err_after * err_after;
    }
    noise_rms_before = sqrt(noise_rms_before / (n - 16));
    noise_rms_after  = sqrt(noise_rms_after / (n - 16));

    printf("  Noise RMS before filtering: %.4f\n", noise_rms_before);
    printf("  Noise RMS after filtering:  %.4f\n", noise_rms_after);
    printf("  Noise reduction: %.1f dB\n",
           20.0 * log10(noise_rms_after / noise_rms_before));

    /* Show a few samples */
    printf("\n  %4s  %10s  %10s  %10s\n", "n", "Original", "Noisy", "Filtered");
    for (int i = 20; i < 36; i++) {
        printf("  %4d  %+10.4f  %+10.4f  %+10.4f\n",
               i, signal[i], noisy[i], filtered[i]);
    }

    printf("\n  ✓ IIR filter effectively removes high-frequency noise.\n");
    printf("  ✓ Filtered output closely tracks the original sine wave.\n");
}

/* ── Main ────────────────────────────────────────────────────────── */

int main(void)
{
    printf("╔══════════════════════════════════════════════════════════╗\n");
    printf("║  Chapter 11: IIR Filter Design                         ║\n");
    printf("║  DSP Tutorial Suite — Interactive Demos                 ║\n");
    printf("╚══════════════════════════════════════════════════════════╝\n");

    demo_butterworth_lp();
    demo_butterworth_hp();
    demo_chebyshev1();
    demo_order_comparison();
    demo_warping();
    demo_iir_vs_fir();
    demo_filter_signal();

    printf("\n════════════════════════════════════════════════════════════\n");
    printf("  All demos complete. See chapters/11-iir-filter-design.md\n");
    printf("  for theory, derivations, and exercises.\n");
    printf("════════════════════════════════════════════════════════════\n\n");

    return 0;
}

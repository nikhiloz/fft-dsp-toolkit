/**
 * Chapter 2 Demo: Sampling, Aliasing & the Nyquist Theorem
 *
 * Interactive demonstrations:
 *   1. Sampling a "continuous" sine at different rates
 *   2. Aliasing: two frequencies producing identical samples
 *   3. Under-sampling vs. proper sampling
 *   4. Critical sampling (f = f_s/2)
 *   5. Normalised frequency mapping
 *   6. Quantization effects (bit-depth reduction)
 *   7. Sinc reconstruction from samples
 *
 * Build:  make chapters && ./build/bin/ch02s
 *
 * References:
 *   - Oppenheim & Willsky, Signals and Systems, Ch 7
 *   - Oppenheim & Schafer, Discrete-Time Signal Processing, Ch 4
 *   - Proakis & Manolakis, Digital Signal Processing, Ch 1
 *   - Lyons, Understanding DSP, Ch 1-2
 *   - Smith, Scientist & Engineer's Guide to DSP, Ch 3
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "signal_gen.h"
#include "dsp_utils.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define N 64  /* default number of samples */

/* ── Helpers ─────────────────────────────────────────────────────── */

static void print_signal(const char *name, const double *x, int n, int count)
{
    if (count > n) count = n;
    printf("\n── %s (%d samples shown) ──\n", name, count);
    for (int i = 0; i < count; i++) {
        printf("  x[%2d] = %+8.4f  ", i, x[i]);
        int bar = (int)(x[i] * 20.0);
        if (bar >= 0) {
            printf("  |");
            for (int b = 0; b < bar && b < 40; b++) printf("█");
        } else {
            int spaces = 20 + bar;
            if (spaces < 0) spaces = 0;
            for (int b = 0; b < spaces; b++) printf(" ");
            for (int b = 0; b < -bar && b < 40; b++) printf("█");
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

/* ── Demo 1: Proper sampling ─────────────────────────────────────── */

static void demo_proper_sampling(void)
{
    print_separator("Demo 1: Proper Sampling (f_s >> 2*f_signal)");

    double fs = 1000.0;    /* sample rate = 1000 Hz */
    double freq = 100.0;   /* signal = 100 Hz sine */
    double signal[N];

    gen_sine(signal, N, 1.0, freq, fs, 0.0);

    printf("\n  Signal: %.0f Hz sine wave\n", freq);
    printf("  Sample rate: %.0f Hz\n", fs);
    printf("  Nyquist freq: %.0f Hz\n", fs / 2.0);
    printf("  Ratio f_s/f = %.1f  (well above 2.0 — safe!)\n", fs / freq);

    print_signal("100 Hz sine @ 1000 Hz sample rate", signal, N, 20);
}

/* ── Demo 2: Aliasing — identical samples from different freqs ──── */

static void demo_aliasing(void)
{
    print_separator("Demo 2: Aliasing — Different Frequencies, Same Samples");

    double fs = 1000.0;
    double signal_a[N], signal_b[N], signal_c[N];

    /* These three frequencies all produce the same samples at fs=1000 */
    double f1 = 300.0;               /* original */
    double f2 = fs - f1;             /* 700 Hz: aliases to 300 Hz */
    double f3 = fs + f1;             /* 1300 Hz: also aliases to 300 Hz */

    gen_sine(signal_a, N, 1.0, f1, fs, 0.0);
    gen_sine(signal_b, N, 1.0, f2, fs, 0.0);
    gen_sine(signal_c, N, 1.0, f3, fs, 0.0);

    printf("\n  Sample rate: %.0f Hz (Nyquist = %.0f Hz)\n", fs, fs / 2.0);
    printf("  Three frequencies:\n");
    printf("    f1 = %.0f Hz  (original — below Nyquist)\n", f1);
    printf("    f2 = %.0f Hz  (= f_s - f1 — above Nyquist, aliases to f1)\n", f2);
    printf("    f3 = %.0f Hz  (= f_s + f1 — far above Nyquist, aliases to f1)\n", f3);

    print_signal("f1 = 300 Hz", signal_a, N, 10);
    print_signal("f2 = 700 Hz", signal_b, N, 10);
    print_signal("f3 = 1300 Hz", signal_c, N, 10);

    /* Show that samples are (approximately) identical */
    printf("\n  Comparison (should be ~identical for f1 and f3, negated for f2):\n");
    double max_diff_13 = 0.0, max_diff_12 = 0.0;
    for (int i = 0; i < N; i++) {
        double d13 = fabs(signal_a[i] - signal_c[i]);
        double d12 = fabs(signal_a[i] + signal_b[i]); /* f2 flips sign */
        if (d13 > max_diff_13) max_diff_13 = d13;
        if (d12 > max_diff_12) max_diff_12 = d12;
    }
    printf("    Max |f1[n] - f3[n]| = %.6e  (≈ 0 → identical!)\n", max_diff_13);
    printf("    Max |f1[n] + f2[n]| = %.6e  (≈ 0 → f2 = -f1)\n", max_diff_12);
}

/* ── Demo 3: Under-sampling vs. Adequate sampling ────────────────── */

static void demo_under_sampling(void)
{
    print_separator("Demo 3: Under-sampling vs. Adequate Sampling");

    double freq = 400.0;
    double signal_ok[N], signal_bad[N];

    /* Good: fs = 1000 Hz  (fs/f = 2.5 > 2) */
    gen_sine(signal_ok, N, 1.0, freq, 1000.0, 0.0);

    /* Bad: fs = 600 Hz  (fs/f = 1.5 < 2 — violates Nyquist!) */
    gen_sine(signal_bad, N, 1.0, freq, 600.0, 0.0);

    printf("\n  Signal: %.0f Hz sine wave\n", freq);
    printf("\n  ✓ Adequate sampling (f_s = 1000 Hz, Nyquist = 500 Hz > 400 Hz):\n");
    print_signal("400 Hz @ 1000 Hz", signal_ok, N, 16);

    printf("\n  ✗ Under-sampling (f_s = 600 Hz, Nyquist = 300 Hz < 400 Hz):\n");
    printf("    400 Hz aliases to |400 - 600| = 200 Hz\n");
    print_signal("400 Hz @ 600 Hz (looks like 200 Hz!)", signal_bad, N, 16);
}

/* ── Demo 4: Critical sampling (f = f_s/2) ──────────────────────── */

static void demo_critical_sampling(void)
{
    print_separator("Demo 4: Critical Sampling (f = f_s / 2)");

    double fs = 1000.0;
    double freq = fs / 2.0; /* exactly Nyquist */
    double signal[N];

    gen_sine(signal, N, 1.0, freq, fs, 0.0);

    printf("\n  f_s = %.0f Hz, f = f_s/2 = %.0f Hz\n", fs, freq);
    printf("  At exactly Nyquist, the signal alternates between values.\n");
    printf("  Phase determines whether we see oscillation or silence!\n");

    print_signal("sine at Nyquist (phase=0)", signal, N, 16);

    /* With phase = pi/2, the cosine at Nyquist alternates ±1 */
    double signal_cos[N];
    gen_cosine(signal_cos, N, 1.0, freq, fs, 0.0);
    print_signal("cosine at Nyquist (alternates ±1)", signal_cos, N, 16);
}

/* ── Demo 5: Normalised frequency ────────────────────────────────── */

static void demo_normalised_freq(void)
{
    print_separator("Demo 5: Normalised Frequency — Same Shape, Different Rates");

    double signal_a[N], signal_b[N];

    /* Two different physical frequencies at different sample rates
       but the SAME normalised frequency → identical samples */
    double f1 = 100.0, fs1 = 1000.0;  /* f/fs = 0.1 */
    double f2 = 500.0, fs2 = 5000.0;  /* f/fs = 0.1 — same! */

    gen_sine(signal_a, N, 1.0, f1, fs1, 0.0);
    gen_sine(signal_b, N, 1.0, f2, fs2, 0.0);

    printf("\n  Signal A: %.0f Hz @ %.0f Hz  → f/f_s = %.2f\n", f1, fs1, f1/fs1);
    printf("  Signal B: %.0f Hz @ %.0f Hz  → f/f_s = %.2f\n", f2, fs2, f2/fs2);
    printf("  Same normalised frequency → identical sample values!\n");

    print_signal("100 Hz @ 1000 Hz (f/f_s = 0.1)", signal_a, N, 10);
    print_signal("500 Hz @ 5000 Hz (f/f_s = 0.1)", signal_b, N, 10);

    double max_diff = 0.0;
    for (int i = 0; i < N; i++) {
        double d = fabs(signal_a[i] - signal_b[i]);
        if (d > max_diff) max_diff = d;
    }
    printf("\n  Max |A[n] - B[n]| = %.6e  (≈ 0 → identical!)\n", max_diff);
}

/* ── Demo 6: Quantization effects ────────────────────────────────── */

static double quantize(double x, int bits)
{
    double levels = (double)(1 << bits);
    double half = levels / 2.0;
    double q = floor(x * half + 0.5) / half;
    if (q > 1.0)  q = 1.0;
    if (q < -1.0) q = -1.0;
    return q;
}

static void demo_quantization(void)
{
    print_separator("Demo 6: Quantization Effects");

    double fs = 1000.0;
    double freq = 50.0;
    double original[N];
    gen_sine(original, N, 0.9, freq, fs, 0.0);

    int bit_depths[] = {3, 4, 8, 16};
    int num_depths = 4;

    for (int b = 0; b < num_depths; b++) {
        int bits = bit_depths[b];
        double quantized[N], noise[N];
        double noise_rms;

        for (int i = 0; i < N; i++) {
            quantized[i] = quantize(original[i], bits);
            noise[i] = original[i] - quantized[i];
        }
        noise_rms = rms(noise, N);

        double snr_db = 20.0 * log10(rms(original, N) / noise_rms);
        double snr_theory = 6.02 * bits + 1.76;

        printf("\n  %d-bit quantization:\n", bits);
        printf("    Levels: %d\n", 1 << bits);
        printf("    Quantization noise RMS: %.6f\n", noise_rms);
        printf("    Measured SNR:  %.1f dB\n", snr_db);
        printf("    Theoretical:   %.1f dB  (6.02B + 1.76)\n", snr_theory);

        if (bits <= 4) {
            print_signal("Original", original, N, 10);
            print_signal("Quantized", quantized, N, 10);
            print_signal("Quant. noise", noise, N, 10);
        }
    }
}

/* ── Demo 7: Sinc reconstruction ─────────────────────────────────── */

static double sinc(double x)
{
    if (fabs(x) < 1e-12) return 1.0;
    return sin(M_PI * x) / (M_PI * x);
}

static void demo_reconstruction(void)
{
    print_separator("Demo 7: Sinc Reconstruction from Samples");

    /* Sample a 50 Hz sine at 200 Hz (4 samples/cycle) */
    int ns = 16; /* number of discrete samples */
    double fs = 200.0;
    double freq = 50.0;
    double samples[16];

    gen_sine(samples, ns, 1.0, freq, fs, 0.0);

    printf("\n  Original: %.0f Hz sine, sampled at %.0f Hz (%d samples)\n",
           freq, fs, ns);
    print_signal("Discrete samples", samples, ns, ns);

    /* Reconstruct at 4x resolution using sinc interpolation */
    int interp_factor = 4;
    int nr = ns * interp_factor;
    double reconstructed[256]; /* ns*interp_factor <= 256 */

    printf("\n  Reconstructing at %dx resolution (%d points) via sinc interpolation...\n",
           interp_factor, nr);

    for (int k = 0; k < nr; k++) {
        double t = (double)k / interp_factor; /* time in sample units */
        double val = 0.0;
        for (int n = 0; n < ns; n++) {
            val += samples[n] * sinc(t - n);
        }
        reconstructed[k] = val;
    }

    /* Compare with the "true" continuous signal at those points */
    printf("\n  Comparing reconstructed vs. true continuous signal:\n");
    printf("  %6s  %10s  %10s  %10s\n", "t", "Recon", "True", "Error");
    double max_err = 0.0;
    for (int k = 0; k < nr; k++) {
        double t_sec = (double)k / (fs * interp_factor);
        double true_val = sin(2.0 * M_PI * freq * t_sec);
        double err = fabs(reconstructed[k] - true_val);
        if (err > max_err) max_err = err;
        if (k % interp_factor == 0 || k < 8) {
            printf("  %6.3f  %+10.6f  %+10.6f  %10.2e%s\n",
                   (double)k / interp_factor,
                   reconstructed[k], true_val, err,
                   (k % interp_factor == 0) ? "  ← sample point" : "");
        }
    }
    printf("  ...\n");
    printf("\n  Max reconstruction error: %.2e\n", max_err);
    printf("  (Edge effects cause larger error near boundaries —\n");
    printf("   the theorem assumes an infinite sequence.)\n");
}

/* ── Main ────────────────────────────────────────────────────────── */

int main(void)
{
    printf("╔══════════════════════════════════════════════════════════╗\n");
    printf("║  Chapter 2: Sampling, Aliasing & the Nyquist Theorem   ║\n");
    printf("║  DSP Tutorial Suite — Interactive Demos                 ║\n");
    printf("╚══════════════════════════════════════════════════════════╝\n");

    demo_proper_sampling();
    demo_aliasing();
    demo_under_sampling();
    demo_critical_sampling();
    demo_normalised_freq();
    demo_quantization();
    demo_reconstruction();

    printf("\n════════════════════════════════════════════════════════════\n");
    printf("  All demos complete. See chapters/02-sampling-and-aliasing.md\n");
    printf("  for theory, proofs, and exercises.\n");
    printf("════════════════════════════════════════════════════════════\n\n");

    return 0;
}

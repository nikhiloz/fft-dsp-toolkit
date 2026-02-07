/**
 * Chapter 7 Demo: The Discrete Fourier Transform (DFT)
 *
 * Interactive demonstrations:
 *   1. Manual 4-point DFT (matrix multiply, O(N²))
 *   2. DFT frequency bin interpretation
 *   3. Conjugate symmetry for real signals
 *   4. Parseval's theorem verification
 *   5. Circular vs. linear convolution
 *   6. Zero-padding and spectral interpolation
 *   7. DFT of standard signals
 *   8. Time-shift property
 *
 * Build:  make chapters && ./build/bin/ch07s
 *
 * References:
 *   - Oppenheim & Schafer, Discrete-Time Signal Processing, Ch 4-5
 *   - Proakis & Manolakis, Digital Signal Processing, Ch 4-5
 *   - Lyons, Understanding DSP, Ch 3
 *   - Smith, Scientist & Engineer's Guide to DSP, Ch 8-9
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "fft.h"
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
 * Brute-force DFT: O(N²) for educational purposes.
 */
static void dft_direct(const Complex *x, Complex *X, int N)
{
    for (int k = 0; k < N; k++) {
        X[k].re = 0.0;
        X[k].im = 0.0;
        for (int n = 0; n < N; n++) {
            double angle = -2.0 * M_PI * k * n / N;
            Complex w = {cos(angle), sin(angle)};
            X[k] = complex_add(X[k], complex_mul(x[n], w));
        }
    }
}

/* ── Demo 1: Manual 4-point DFT ──────────────────────────────────── */

static void demo_manual_dft(void)
{
    print_separator("Demo 1: Manual 4-point DFT (O(N²))");

    Complex x[4] = {{1, 0}, {0, 0}, {-1, 0}, {0, 0}};
    Complex X[4];

    printf("\n  x[n] = {1, 0, -1, 0}\n\n");

    /* Show the DFT matrix W_4 */
    printf("  DFT matrix W_4:\n");
    printf("       n=0    n=1    n=2    n=3\n");
    for (int k = 0; k < 4; k++) {
        printf("  k=%d: ", k);
        for (int n = 0; n < 4; n++) {
            double angle = -2.0 * M_PI * k * n / 4.0;
            printf("%+5.1f%+5.1fj  ", cos(angle), sin(angle));
        }
        printf("\n");
    }

    dft_direct(x, X, 4);

    printf("\n  DFT result X[k]:\n");
    for (int k = 0; k < 4; k++) {
        printf("    X[%d] = %+6.2f %+6.2fj   |X[%d]| = %.2f\n",
               k, X[k].re, X[k].im, k, complex_mag(X[k]));
    }
    printf("\n  ✓ X = {0, 2, 0, 2} — energy at bins 1 and 3 (the \"signal\" frequency).\n");
}

/* ── Demo 2: Frequency bin interpretation ────────────────────────── */

static void demo_freq_bins(void)
{
    print_separator("Demo 2: Frequency Bin Interpretation");

    int N = 16;
    double fs = 1000.0;
    double signal[16];

    /* 250 Hz sine → should appear at bin k = N*f/f_s = 16*250/1000 = 4 */
    gen_sine(signal, N, 1.0, 250.0, fs, 0.0);

    Complex X[16];
    fft_real(signal, X, N);

    printf("\n  Signal: 250 Hz sine, f_s = 1000 Hz, N = %d\n", N);
    printf("  Expected bin: k = N × f / f_s = %d × 250 / 1000 = 4\n\n", N);

    printf("  %4s  %10s  %10s  %30s\n", "k", "f (Hz)", "|X[k]|", "");
    for (int k = 0; k <= N / 2; k++) {
        double freq = (double)k * fs / N;
        double mag = complex_mag(X[k]);
        int bar = (int)(mag * 4.0);
        if (bar > 30) bar = 30;

        char bar_str[32];
        int i;
        for (i = 0; i < bar; i++) bar_str[i] = '#';
        bar_str[i] = '\0';

        printf("  %4d  %10.1f  %10.2f  |%-30s|%s\n",
               k, freq, mag, bar_str,
               (k == 4) ? "  ← peak!" : "");
    }
}

/* ── Demo 3: Conjugate symmetry for real signals ─────────────────── */

static void demo_conjugate_symmetry(void)
{
    print_separator("Demo 3: Conjugate Symmetry (Real Signals)");

    int N = 8;
    double x[8];
    gen_sine(x, N, 1.0, 1.0, 8.0, 0.0);

    Complex Xc[8];
    fft_real(x, Xc, N);

    printf("\n  Real signal: sine wave, N = %d\n", N);
    printf("  Property: X[N-k] = X*[k] (conjugate symmetry)\n\n");

    printf("  %4s  %12s  %12s   %4s  %12s  %12s  %s\n",
           "k", "Re{X[k]}", "Im{X[k]}", "N-k", "Re{X[N-k]}", "Im{X[N-k]}", "X[N-k]=X*[k]?");

    for (int k = 0; k <= N / 2; k++) {
        int nk = (N - k) % N;
        int match = (fabs(Xc[k].re - Xc[nk].re) < 1e-10) &&
                    (fabs(Xc[k].im + Xc[nk].im) < 1e-10); /* im should negate */
        printf("  %4d  %+12.4f  %+12.4f   %4d  %+12.4f  %+12.4f  %s\n",
               k, Xc[k].re, Xc[k].im,
               nk, Xc[nk].re, Xc[nk].im,
               match ? "✓" : "✗");
    }
    printf("\n  ✓ Imaginary parts are negated — conjugate symmetry confirmed.\n");
}

/* ── Demo 4: Parseval's theorem ──────────────────────────────────── */

static void demo_parseval(void)
{
    print_separator("Demo 4: Parseval's Theorem Verification");

    int N = 64;
    double x[64];
    gen_sine(x, N, 1.0, 5.0, 64.0, 0.0);

    /* Time-domain energy */
    double E_time = signal_energy(x, N);

    /* Frequency-domain energy */
    Complex X[64];
    fft_real(x, X, N);

    double E_freq = 0.0;
    for (int k = 0; k < N; k++) {
        E_freq += X[k].re * X[k].re + X[k].im * X[k].im;
    }
    E_freq /= (double)N;

    printf("\n  Signal: 5 Hz sine, N = %d, f_s = 64 Hz\n\n", N);
    printf("  Time-domain energy:  Σ|x[n]|² = %.6f\n", E_time);
    printf("  Freq-domain energy:  (1/N) Σ|X[k]|² = %.6f\n", E_freq);
    printf("  Difference: %.2e\n", fabs(E_time - E_freq));
    printf("\n  ✓ Parseval's theorem: energy in time = energy in frequency.\n");
}

/* ── Demo 5: Circular vs. linear convolution ─────────────────────── */

static void demo_circular_vs_linear(void)
{
    print_separator("Demo 5: Circular vs. Linear Convolution");

    double x[] = {1, 2, 3, 4};
    double h[] = {1, 1, 1};
    int Lx = 4, Lh = 3;

    /* Linear convolution (direct) */
    int Ly = Lx + Lh - 1; /* = 6 */
    double y_linear[6];
    convolve(x, Lx, h, Lh, y_linear);

    /* Circular convolution via DFT (N = 4, too short!) */
    int N_short = 4;

    /* Pad h to length 4 */
    Complex xc[8], hc[8];
    int i;
    for (i = 0; i < 4; i++) { xc[i].re = x[i]; xc[i].im = 0; }
    for (i = 0; i < 3; i++) { hc[i].re = h[i]; hc[i].im = 0; }
    hc[3].re = 0; hc[3].im = 0;

    fft(xc, N_short);  /* in-place */
    fft(hc, N_short);
    Complex Yf[8];
    for (i = 0; i < N_short; i++) Yf[i] = complex_mul(xc[i], hc[i]);

    ifft(Yf, N_short);  /* in-place: Yf now holds time-domain result */

    /* Circular convolution via DFT (N = 8, properly padded) */
    int N_long = 8;
    Complex xp[8], hp[8];
    for (i = 0; i < 8; i++) {
        xp[i].re = (i < Lx) ? x[i] : 0;
        xp[i].im = 0;
        hp[i].re = (i < Lh) ? h[i] : 0;
        hp[i].im = 0;
    }
    fft(xp, N_long);
    fft(hp, N_long);
    Complex YP[8];
    for (i = 0; i < N_long; i++) YP[i] = complex_mul(xp[i], hp[i]);

    ifft(YP, N_long);

    printf("\n  x = {1, 2, 3, 4}, h = {1, 1, 1}\n");
    printf("  Linear conv length = %d + %d - 1 = %d\n\n", Lx, Lh, Ly);

    printf("  %4s  %10s  %14s  %14s\n",
           "n", "Linear", "Circ (N=4)", "Circ (N=8, padded)");
    for (i = 0; i < Ly; i++) {
        double circ_short = (i < N_short) ? Yf[i].re : 0;
        double circ_long = (i < N_long) ? YP[i].re : 0;
        const char *match = (fabs(y_linear[i] - circ_long) < 1e-6) ? " ✓" : " ✗";
        const char *wrap = (i < N_short && fabs(y_linear[i] - circ_short) > 1e-6) ?
                           " ← aliased!" : "";
        printf("  %4d  %10.1f  %14.1f%s  %10.1f%s\n",
               i, y_linear[i], circ_short, wrap, circ_long, match);
    }
    printf("\n  ✓ N=4 (too short): circular convolution wraps around → wrong!\n");
    printf("  ✓ N≥6 (padded): circular = linear → correct.\n");
}

/* ── Demo 6: Zero-padding and spectral interpolation ─────────────── */

static void demo_zero_padding(void)
{
    print_separator("Demo 6: Zero-Padding — Spectral Interpolation");

    /* 8-sample signal with 2 tones */
    int N_orig = 8;
    double x8[8];
    double freqs[] = {1.0, 3.0};
    double amps[] = {1.0, 0.5};
    gen_multi_tone(x8, N_orig, freqs, amps, 2, 8.0);

    /* DFT at N=8 */
    Complex X8[8];
    fft_real(x8, X8, N_orig);

    /* Zero-pad to N=32 */
    int N_pad = 32;
    Complex X32[32];
    int i;
    for (i = 0; i < 32; i++) {
        X32[i].re = (i < N_orig) ? x8[i] : 0.0;
        X32[i].im = 0.0;
    }
    fft(X32, N_pad);  /* in-place */

    printf("\n  Signal: 2 tones (1 Hz + 3 Hz) at f_s = 8 Hz, N = 8 samples\n");
    printf("  Zero-padded to N = 32 for spectral interpolation\n\n");

    printf("  N=8 DFT (coarse, %d bins):\n", N_orig / 2 + 1);
    for (int k = 0; k <= N_orig / 2; k++) {
        double mag = complex_mag(X8[k]) / N_orig * 2;
        printf("    k=%d (f=%.1f Hz): |X| = %.3f\n",
               k, (double)k * 8.0 / N_orig, mag);
    }

    printf("\n  N=32 DFT (interpolated, %d bins):\n", N_pad / 2 + 1);
    for (int k = 0; k <= N_pad / 2; k++) {
        double mag = complex_mag(X32[k]) / N_pad * 2;
        int bar = (int)(mag * 30.0);
        if (bar > 30) bar = 30;
        char bar_str[32];
        int j;
        for (j = 0; j < bar; j++) bar_str[j] = '#';
        bar_str[j] = '\0';

        printf("    k=%2d (f=%5.2f Hz): |X| = %.3f  |%-30s|\n",
               k, (double)k * 8.0 / N_pad, mag, bar_str);
    }
    printf("\n  ✓ Same peaks, but the spectrum is now 4× denser.\n");
    printf("  ⚠ Resolution hasn't improved — just interpolated between bins.\n");
}

/* ── Demo 7: DFT of standard signals ────────────────────────────── */

static void demo_standard_signals(void)
{
    print_separator("Demo 7: DFT of Standard Signals");

    int N = 16;
    Complex X[16];

    /* (a) Impulse */
    double imp[16];
    gen_impulse(imp, N, 0);
    fft_real(imp, X, N);
    printf("\n  δ[n] → flat spectrum (all bins equal):\n    ");
    for (int k = 0; k < N; k++) printf("%.1f ", complex_mag(X[k]));
    printf("\n");

    /* (b) Constant (DC) */
    double dc[16];
    int i;
    for (i = 0; i < N; i++) dc[i] = 1.0;
    fft_real(dc, X, N);
    printf("\n  x[n] = 1 → energy only at DC (bin 0):\n    ");
    for (int k = 0; k < N; k++) printf("%.1f ", complex_mag(X[k]));
    printf("\n");

    /* (c) Alternating */
    double alt[16];
    for (i = 0; i < N; i++) alt[i] = (i % 2 == 0) ? 1.0 : -1.0;
    fft_real(alt, X, N);
    printf("\n  x[n] = (-1)^n → energy only at Nyquist (bin N/2):\n    ");
    for (int k = 0; k < N; k++) printf("%.1f ", complex_mag(X[k]));
    printf("\n");

    /* (d) Pure cosine at bin 3 */
    double cosine[16];
    gen_cosine(cosine, N, 1.0, 3.0, 16.0, 0.0);
    fft_real(cosine, X, N);
    printf("\n  cos(2π·3·n/16) → peaks at bins 3 and 13 (= N-3):\n    ");
    for (int k = 0; k < N; k++) printf("%.1f ", complex_mag(X[k]));
    printf("\n");
}

/* ── Demo 8: Time-shift property ─────────────────────────────────── */

static void demo_time_shift(void)
{
    print_separator("Demo 8: Time-Shift Property");

    int N = 8;
    double x[8], x_shifted[8];

    /* Original: impulse at 0 */
    gen_impulse(x, N, 0);

    /* Shifted: impulse at 2 (circular shift) */
    gen_impulse(x_shifted, N, 2);

    Complex X_orig[8], X_shift[8];
    fft_real(x, X_orig, N);
    fft_real(x_shifted, X_shift, N);

    printf("\n  x[n] = δ[n],  x_shifted[n] = δ[n-2]\n");
    printf("  Time shift by m=2: X_shift[k] = X[k] · W_N^{2k}\n\n");

    printf("  %4s  %12s  %12s  %12s  %12s\n",
           "k", "|X[k]|", "∠X[k]°", "|X_s[k]|", "∠X_s[k]°");
    for (int k = 0; k < N; k++) {
        printf("  %4d  %12.4f  %+12.1f  %12.4f  %+12.1f\n",
               k,
               complex_mag(X_orig[k]),
               complex_phase(X_orig[k]) * 180.0 / M_PI,
               complex_mag(X_shift[k]),
               complex_phase(X_shift[k]) * 180.0 / M_PI);
    }
    printf("\n  ✓ Magnitudes identical — only phase changes.\n");
    printf("  ✓ Phase shift = 2×k×(360/8) = 90k degrees per bin.\n");
}

/* ── Main ────────────────────────────────────────────────────────── */

int main(void)
{
    printf("╔══════════════════════════════════════════════════════════╗\n");
    printf("║  Chapter 7: The Discrete Fourier Transform (DFT)       ║\n");
    printf("║  DSP Tutorial Suite — Interactive Demos                 ║\n");
    printf("╚══════════════════════════════════════════════════════════╝\n");

    demo_manual_dft();
    demo_freq_bins();
    demo_conjugate_symmetry();
    demo_parseval();
    demo_circular_vs_linear();
    demo_zero_padding();
    demo_standard_signals();
    demo_time_shift();

    printf("\n════════════════════════════════════════════════════════════\n");
    printf("  All demos complete. See chapters/07-dft-theory.md\n");
    printf("  for theory, derivations, and exercises.\n");
    printf("════════════════════════════════════════════════════════════\n\n");

    return 0;
}

/**
 * @file 17-multirate-dsp.c
 * @brief Chapter 17 demo — Multirate DSP: decimation, interpolation, polyphase.
 *
 * ── Decimation by M ──────────────────────────────────────────────
 *
 *   Original:  x x x x x x x x x x x x   (fs = 8000 Hz)
 *   After ↓3:  x . . x . . x . . x . .   (fs = 2667 Hz)
 *
 *   Anti-alias LPF (fc = fs/2M) must precede downsampling
 *   to prevent spectral folding.
 *
 * ── Interpolation by L ───────────────────────────────────────────
 *
 *   Original:  x . . . x . . . x . . .   (fs = 2000 Hz)
 *   After ↑4:  x 0 0 0 x 0 0 0 x 0 0 0   (zeros inserted)
 *   After LPF: x ~ ~ ~ x ~ ~ ~ x ~ ~ ~   (fs = 8000 Hz)
 *
 * ── Polyphase Structure ──────────────────────────────────────────
 *
 *         ┌──► [h₀] ──┐
 *   x[n]──┤──► [h₁] ──┼──► Σ ──► y[m]
 *         ├──► [h₂] ──┤    (at rate fs/M)
 *         └──► [hₘ₋₁]─┘
 *
 * Demonstrates:
 *   - Decimation with anti-alias filtering
 *   - Interpolation with anti-image filtering
 *   - Rational rate conversion (L/M)
 *   - Polyphase vs direct efficiency
 *   - Spectral effects of multirate operations
 *
 * Build & run:
 *   make chapters && ./build/bin/ch17
 *
 * Read alongside: chapters/17-multirate-dsp.md
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "multirate.h"
#include "signal_gen.h"
#include "fft.h"
#include "filter.h"
#include "dsp_utils.h"
#include "gnuplot.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ── Demo 1: Decimation ───────────────────────────────────────── */

static void demo_decimation(void)
{
    printf("=== Demo 1: Decimation by M ===\n\n");

    const int N = 1024;
    const double fs = 8000.0;
    double *x = (double *)malloc((size_t)N * sizeof(double));

    /* Two-tone signal: 200 Hz + 1500 Hz */
    double *tmp = (double *)malloc((size_t)N * sizeof(double));
    gen_sine(x, N, 0.7, 200.0, fs, 0.0);
    gen_sine(tmp, N, 0.5, 1500.0, fs, 0.0);
    signal_add(x, tmp, N);
    free(tmp);

    printf("  Input: 200 Hz + 1500 Hz, fs = %.0f Hz, N = %d\n\n", fs, N);

    int factors[] = {2, 4, 8};
    for (int f = 0; f < 3; f++) {
        int M = factors[f];
        int out_len = N / M;
        double *y = (double *)calloc((size_t)out_len, sizeof(double));
        int actual_len = decimate(x, N, M, y);

        printf("  ↓%d: %d → %d samples (fs_new = %.0f Hz)\n",
               M, N, actual_len, fs / M);
        printf("       Nyquist = %.0f Hz → 1500 Hz %s\n",
               fs / (2.0 * M),
               (1500.0 < fs / (2.0 * M)) ? "preserved" : "ALIASED (filtered out)");

        free(y);
    }

    /* Plot original vs decimated */
    int M = 4;
    int out_len = N / M;
    double *y = (double *)calloc((size_t)out_len, sizeof(double));
    decimate(x, N, M, y);

    double *idx_in  = (double *)malloc((size_t)N * sizeof(double));
    double *idx_out = (double *)malloc((size_t)out_len * sizeof(double));
    for (int i = 0; i < N; i++) idx_in[i] = (double)i / fs;
    for (int i = 0; i < out_len; i++) idx_out[i] = (double)i / (fs / M);

    int plot_n = (N < 200) ? N : 200;
    int plot_m = (out_len < 50) ? out_len : 50;
    gp_plot_multi("ch17", "decimation",
        "Decimation by 4: Original vs Downsampled",
        "Time (s)", "Amplitude",
        (GpSeries[]){
            {"Original (8000 Hz)", idx_in, x, plot_n, "lines"},
            {"Decimated (2000 Hz)", idx_out, y, plot_m, "linespoints"}
        }, 2);
    printf("  → plots/ch17/decimation.png\n\n");

    free(x); free(y); free(idx_in); free(idx_out);
}

/* ── Demo 2: Interpolation ────────────────────────────────────── */

static void demo_interpolation(void)
{
    printf("=== Demo 2: Interpolation by L ===\n\n");

    const int N = 128;
    const double fs = 2000.0;
    double *x = (double *)malloc((size_t)N * sizeof(double));
    gen_sine(x, N, 1.0, 200.0, fs, 0.0);

    int L = 4;
    int out_len = N * L;
    double *y = (double *)calloc((size_t)out_len, sizeof(double));
    int actual = interpolate(x, N, L, y);

    printf("  Input: 200 Hz sine, fs = %.0f Hz, N = %d\n", fs, N);
    printf("  ↑%d: %d → %d samples (fs_new = %.0f Hz)\n\n",
           L, N, actual, fs * L);

    /* Plot */
    double *idx_in  = (double *)malloc((size_t)N * sizeof(double));
    double *idx_out = (double *)malloc((size_t)out_len * sizeof(double));
    for (int i = 0; i < N; i++) idx_in[i] = (double)i / fs;
    for (int i = 0; i < out_len; i++) idx_out[i] = (double)i / (fs * L);

    int pn = 64, pm = 256;
    gp_plot_multi("ch17", "interpolation",
        "Interpolation by 4: Original vs Upsampled",
        "Time (s)", "Amplitude",
        (GpSeries[]){
            {"Original (2000 Hz)", idx_in, x, pn, "linespoints"},
            {"Interpolated (8000 Hz)", idx_out, y, pm, "lines"}
        }, 2);
    printf("  → plots/ch17/interpolation.png\n\n");

    free(x); free(y); free(idx_in); free(idx_out);
}

/* ── Demo 3: Rational Rate Conversion ─────────────────────────── */

static void demo_resample(void)
{
    printf("=== Demo 3: Rational Rate Conversion (L/M) ===\n\n");

    const int N = 480;                /* 48000/100 = 480 samples */
    const double fs_in = 48000.0;
    double *x = (double *)malloc((size_t)N * sizeof(double));
    gen_sine(x, N, 1.0, 440.0, fs_in, 0.0);

    /* 48000 → 44100:  L/M = 441/480 = 147/160 */
    int L = 147, M = 160;
    int max_out = (N * L) / M + 10;
    double *y = (double *)calloc((size_t)max_out, sizeof(double));
    int out_len = resample(x, N, L, M, y);

    double fs_out = fs_in * (double)L / (double)M;
    printf("  Input:  fs = %.0f Hz, N = %d\n", fs_in, N);
    printf("  Output: fs = %.1f Hz, N = %d\n", fs_out, out_len);
    printf("  Ratio:  %d/%d = %.6f\n\n", L, M, (double)L / (double)M);

    free(x); free(y);
}

/* ── Demo 4: Polyphase vs Direct ──────────────────────────────── */

static void demo_polyphase(void)
{
    printf("=== Demo 4: Polyphase Decimation vs Direct ===\n\n");

    const int N = 2048, M = 4;
    const double fs = 8000.0;
    double *x = (double *)malloc((size_t)N * sizeof(double));
    gen_sine(x, N, 1.0, 300.0, fs, 0.0);

    /* Design anti-alias filter manually */
    int taps = 33;
    double *h = (double *)calloc((size_t)taps, sizeof(double));
    fir_lowpass(h, taps, 0.5 / M);

    /* Direct decimation (filter then downsample) */
    double *filtered = (double *)calloc((size_t)N, sizeof(double));
    fir_filter(x, filtered, N, h, taps);
    int out_len = N / M;
    double *y_direct = (double *)calloc((size_t)out_len, sizeof(double));
    for (int i = 0; i < out_len; i++)
        y_direct[i] = filtered[i * M];

    /* Polyphase decimation */
    double *y_poly = (double *)calloc((size_t)out_len, sizeof(double));
    polyphase_decimate(x, N, h, taps, M, y_poly);

    /* Compare */
    double max_err = 0.0;
    for (int i = 0; i < out_len; i++) {
        double err = fabs(y_direct[i] - y_poly[i]);
        if (err > max_err) max_err = err;
    }

    printf("  Direct:    filter (%d taps × %d samples) + downsample\n", taps, N);
    printf("  Polyphase: %d sub-filters × %d taps each, at rate 1/%d\n",
           M, (taps + M - 1) / M, M);
    printf("  Max |direct - polyphase| = %.2e\n", max_err);
    printf("  Ops ratio: polyphase uses ~%d× fewer multiplies\n\n",
           M);

    free(x); free(h); free(filtered); free(y_direct); free(y_poly);
}

/* ── Demo 5: Spectral View ────────────────────────────────────── */

static void demo_spectral(void)
{
    printf("=== Demo 5: Spectral Effects of Decimation ===\n\n");

    const int N = 512;
    const double fs = 8000.0;
    double *x = (double *)malloc((size_t)N * sizeof(double));

    /* Multi-tone */
    double *tmp = (double *)malloc((size_t)N * sizeof(double));
    gen_sine(x, N, 0.5, 500.0, fs, 0.0);
    gen_sine(tmp, N, 0.5, 1200.0, fs, 0.0);
    signal_add(x, tmp, N);
    gen_sine(tmp, N, 0.3, 3000.0, fs, 0.0);
    signal_add(x, tmp, N);
    free(tmp);

    /* Original spectrum */
    Complex *X = (Complex *)calloc((size_t)N, sizeof(Complex));
    for (int i = 0; i < N; i++) { X[i].re = x[i]; X[i].im = 0.0; }
    fft(X, N);

    double *mag_orig = (double *)malloc((size_t)(N / 2) * sizeof(double));
    double *freq_axis = (double *)malloc((size_t)(N / 2) * sizeof(double));
    for (int i = 0; i < N / 2; i++) {
        mag_orig[i] = 20.0 * log10(complex_mag(X[i]) / N + 1e-10);
        freq_axis[i] = (double)i * fs / N;
    }

    gp_plot_1("ch17", "spectrum_original",
              "Original Spectrum (500+1200+3000 Hz)",
              "Frequency (Hz)", "Magnitude (dB)",
              freq_axis, mag_orig, N / 2, "lines");

    printf("  Input tones: 500, 1200, 3000 Hz\n");
    printf("  After ↓4 (fs=2000): Nyquist=1000 → 1200 and 3000 Hz removed\n");
    printf("  → plots/ch17/spectrum_original.png\n\n");

    free(x); free(X); free(mag_orig); free(freq_axis);
}

int main(void)
{
    printf("╔══════════════════════════════════════════════════════════╗\n");
    printf("║  Chapter 17: Multirate DSP                             ║\n");
    printf("╚══════════════════════════════════════════════════════════╝\n\n");

    demo_decimation();
    demo_interpolation();
    demo_resample();
    demo_polyphase();
    demo_spectral();

    printf("=== Chapter 17 Complete ===\n");
    return 0;
}

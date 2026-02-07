/**
 * @file 22-advanced-fir.c
 * @brief Chapter 22 demo — Advanced FIR Design: Parks-McClellan (Remez).
 *
 * ── Equiripple vs Windowed FIR ───────────────────────────────────
 *
 *   Window method:         Remez (equiripple):
 *
 *   Passband ─────         Passband ≈≈≈≈≈ (equal ripple)
 *        ↑                      ↑
 *        │ smooth roll-off      │ sharp transition
 *        │                      │
 *   Stopband ┈┈┈┈┈─        Stopband ≈≈≈≈≈ (equal ripple)
 *
 *   For the same filter order, Remez achieves a narrower
 *   transition band (sharper cutoff).
 *
 * ── Alternation Theorem ──────────────────────────────────────────
 *
 *   An optimal Chebyshev FIR of order N has at least N+2
 *   equiripple extrema in the error function.
 *
 *   Error: ···•···•···•···•···•···•···  (alternating ±δ)
 *
 * Demonstrates:
 *   - Remez lowpass design
 *   - Window method vs Remez comparison
 *   - Bandpass equiripple filter
 *   - Transition band trade-off (order vs width)
 *   - Passband/stopband ripple control
 *
 * Build & run:
 *   make chapters && ./build/bin/ch22
 *
 * Read alongside: chapters/22-advanced-fir.md
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "remez.h"
#include "filter.h"
#include "fft.h"
#include "dsp_utils.h"
#include "gnuplot.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ── Compute frequency response magnitude ─────────────────────── */

static void freq_response_db(const double *h, int taps, int nfft,
                             double *mag_db, double *freq_axis, double fs)
{
    Complex *H = (Complex *)calloc((size_t)nfft, sizeof(Complex));
    for (int i = 0; i < taps && i < nfft; i++) {
        H[i].re = h[i];
        H[i].im = 0.0;
    }
    fft(H, nfft);
    for (int i = 0; i < nfft / 2; i++) {
        mag_db[i] = 20.0 * log10(complex_mag(H[i]) + 1e-12);
        freq_axis[i] = (double)i * fs / nfft;
    }
    free(H);
}

/* ── Demo 1: Remez Lowpass Design ─────────────────────────────── */

static void demo_remez_lowpass(void)
{
    printf("=== Demo 1: Remez Lowpass Design ===\n\n");

    int taps = 51;
    double *h = (double *)calloc((size_t)taps, sizeof(double));

    /* Normalised: passband edge 0.2, stopband edge 0.3 */
    int ret = remez_lowpass(h, taps, 0.2, 0.3, 1.0, 1.0);

    printf("  %d-tap equiripple lowpass\n", taps);
    printf("  Passband: 0 — 0.2·fs/2  Stopband: 0.3·fs/2 — fs/2\n");
    printf("  Remez converged: %s\n\n", ret == 0 ? "YES" : "NO");

    /* Plot frequency response */
    int nfft = 1024;
    double *mag = (double *)malloc((size_t)(nfft / 2) * sizeof(double));
    double *freq = (double *)malloc((size_t)(nfft / 2) * sizeof(double));
    freq_response_db(h, taps, nfft, mag, freq, 2.0); /* fs=2 → normalised */

    gp_plot_1("ch22", "remez_lowpass",
        "51-tap Remez Equiripple Lowpass (fpass=0.2, fstop=0.3)",
        "Normalised Frequency (×π)", "Magnitude (dB)",
        freq, mag, nfft / 2, "lines");
    printf("  → plots/ch22/remez_lowpass.png\n\n");

    free(h); free(mag); free(freq);
}

/* ── Demo 2: Window vs Remez Comparison ───────────────────────── */

static void demo_window_vs_remez(void)
{
    printf("=== Demo 2: Window Method vs Remez ===\n\n");

    int taps = 41;
    double cutoff = 0.25; /* normalised cutoff */

    /* Window method (Hamming) — using existing fir_lowpass */
    double *h_win = (double *)calloc((size_t)taps, sizeof(double));
    fir_lowpass(h_win, taps, cutoff);

    /* Remez method */
    double *h_remez = (double *)calloc((size_t)taps, sizeof(double));
    remez_lowpass(h_remez, taps, cutoff, cutoff + 0.1, 1.0, 1.0);

    /* Compute both frequency responses */
    int nfft = 1024;
    double *mag_win  = (double *)malloc((size_t)(nfft / 2) * sizeof(double));
    double *mag_remez = (double *)malloc((size_t)(nfft / 2) * sizeof(double));
    double *freq     = (double *)malloc((size_t)(nfft / 2) * sizeof(double));

    freq_response_db(h_win, taps, nfft, mag_win, freq, 2.0);
    freq_response_db(h_remez, taps, nfft, mag_remez, freq, 2.0);

    /* Find stopband attenuation */
    int stop_bin = (int)(0.35 * nfft / 2);
    double max_stop_win = -200.0, max_stop_remez = -200.0;
    for (int i = stop_bin; i < nfft / 2; i++) {
        if (mag_win[i] > max_stop_win) max_stop_win = mag_win[i];
        if (mag_remez[i] > max_stop_remez) max_stop_remez = mag_remez[i];
    }

    printf("  Both: %d taps, normalised cutoff = %.2f\n", taps, cutoff);
    printf("  Stopband attenuation:\n");
    printf("    Window (Hamming): %.1f dB\n", max_stop_win);
    printf("    Remez:            %.1f dB\n", max_stop_remez);
    printf("  Remez should have sharper transition band.\n\n");

    gp_plot_multi("ch22", "window_vs_remez",
        "41-tap Filter: Window Method vs Remez",
        "Normalised Frequency (×π)", "Magnitude (dB)",
        (GpSeries[]){
            {"Window (Hamming)",   freq, mag_win, nfft / 2, "lines"},
            {"Remez (equiripple)", freq, mag_remez, nfft / 2, "lines"}
        }, 2);
    printf("  → plots/ch22/window_vs_remez.png\n\n");

    free(h_win); free(h_remez);
    free(mag_win); free(mag_remez); free(freq);
}

/* ── Demo 3: Bandpass Equiripple ──────────────────────────────── */

static void demo_bandpass(void)
{
    printf("=== Demo 3: Bandpass Equiripple Filter ===\n\n");

    int taps = 61;
    double *h = (double *)calloc((size_t)taps, sizeof(double));

    int ret = remez_bandpass(h, taps, 0.1, 0.2, 0.35, 0.45);

    printf("  %d-tap equiripple bandpass\n", taps);
    printf("  Stopband1: 0—0.1, Passband: 0.2—0.35, Stopband2: 0.45—0.5\n");
    printf("  Converged: %s\n\n", ret == 0 ? "YES" : "NO");

    int nfft = 1024;
    double *mag = (double *)malloc((size_t)(nfft / 2) * sizeof(double));
    double *freq = (double *)malloc((size_t)(nfft / 2) * sizeof(double));
    freq_response_db(h, taps, nfft, mag, freq, 2.0);

    gp_plot_1("ch22", "remez_bandpass",
        "61-tap Remez Bandpass (pass: 0.2—0.35)",
        "Normalised Frequency (×π)", "Magnitude (dB)",
        freq, mag, nfft / 2, "lines");
    printf("  → plots/ch22/remez_bandpass.png\n\n");

    free(h); free(mag); free(freq);
}

/* ── Demo 4: Order vs Transition Width ────────────────────────── */

static void demo_order_tradeoff(void)
{
    printf("=== Demo 4: Filter Order vs Transition Width ===\n\n");

    double fpass = 0.2;
    int orders[] = {21, 41, 61, 81};
    double trans[] = {0.15, 0.10, 0.07, 0.05};

    printf("  All: passband edge = %.2f\n", fpass);
    printf("  ┌───────┬──────────────┬───────────────┐\n");
    printf("  │ Taps  │ Δf (trans)   │ Stopband (dB) │\n");
    printf("  ├───────┼──────────────┼───────────────┤\n");

    int nfft = 1024;
    double *mag = (double *)malloc((size_t)(nfft / 2) * sizeof(double));
    double *freq = (double *)malloc((size_t)(nfft / 2) * sizeof(double));

    for (int i = 0; i < 4; i++) {
        double *h = (double *)calloc((size_t)orders[i], sizeof(double));
        remez_lowpass(h, orders[i], fpass, fpass + trans[i], 1.0, 1.0);

        freq_response_db(h, orders[i], nfft, mag, freq, 2.0);

        /* Find max stopband level */
        int stop_bin = (int)((fpass + trans[i]) * nfft / 2);
        double max_stop = -200.0;
        for (int k = stop_bin; k < nfft / 2; k++) {
            if (mag[k] > max_stop) max_stop = mag[k];
        }

        printf("  │  %3d  │    %.2f       │    %+.1f      │\n",
               orders[i], trans[i], max_stop);
        free(h);
    }
    printf("  └───────┴──────────────┴───────────────┘\n");
    printf("  More taps → narrower transition → better attenuation\n\n");

    free(mag); free(freq);
}

/* ── Demo 5: Weighted Ripple Control ──────────────────────────── */

static void demo_weighted(void)
{
    printf("=== Demo 5: Weighted Ripple Control ===\n\n");

    int taps = 51;
    int nfft = 1024;
    double *mag1 = (double *)malloc((size_t)(nfft / 2) * sizeof(double));
    double *mag2 = (double *)malloc((size_t)(nfft / 2) * sizeof(double));
    double *freq = (double *)malloc((size_t)(nfft / 2) * sizeof(double));

    /* Equal weights */
    double *h1 = (double *)calloc((size_t)taps, sizeof(double));
    remez_lowpass(h1, taps, 0.2, 0.3, 1.0, 1.0);
    freq_response_db(h1, taps, nfft, mag1, freq, 2.0);

    /* Heavy stopband weight (want better stopband rejection) */
    double *h2 = (double *)calloc((size_t)taps, sizeof(double));
    remez_lowpass(h2, taps, 0.2, 0.3, 1.0, 10.0);
    freq_response_db(h2, taps, nfft, mag2, freq, 2.0);

    printf("  Both: %d taps, fpass=0.2, fstop=0.3\n", taps);
    printf("  Equal weight (1:1): passband ripple and stopband attenuation balanced\n");
    printf("  Heavy stopband (1:10): more stopband rejection, larger passband ripple\n\n");

    gp_plot_multi("ch22", "weighted_ripple",
        "Weight Trade-off: Equal (1:1) vs Heavy Stopband (1:10)",
        "Normalised Frequency (×π)", "Magnitude (dB)",
        (GpSeries[]){
            {"Equal weight (1:1)", freq, mag1, nfft / 2, "lines"},
            {"Stopband ×10",       freq, mag2, nfft / 2, "lines"}
        }, 2);
    printf("  → plots/ch22/weighted_ripple.png\n\n");

    free(h1); free(h2); free(mag1); free(mag2); free(freq);
}

int main(void)
{
    printf("╔══════════════════════════════════════════════════════════╗\n");
    printf("║  Chapter 22: Advanced FIR Design (Remez / Parks-McClellan)║\n");
    printf("╚══════════════════════════════════════════════════════════╝\n\n");

    demo_remez_lowpass();
    demo_window_vs_remez();
    demo_bandpass();
    demo_order_tradeoff();
    demo_weighted();

    printf("=== Chapter 22 Complete ===\n");
    return 0;
}

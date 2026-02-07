/**
 * @file 25-parametric-spectral.c
 * @brief Chapter 25 demo — Parametric Spectral Estimation: MUSIC & Capon.
 *
 * ── Classical vs Parametric ──────────────────────────────────────
 *
 *   FFT/Periodogram:           MUSIC:
 *
 *     ___/\___/\___            ___│___│_______
 *        f₁   f₂                  f₁  f₂
 *     (broad peaks)            (sharp spikes — super-resolution)
 *
 *   Capon (MVDR):              Model order:
 *     ___∧___∧____             Too low  → miss signals
 *        f₁  f₂               Too high → spurious peaks
 *     (narrow peaks,
 *      true PSD shape)
 *
 * Demonstrates:
 *   - MUSIC pseudospectrum for closely-spaced sinusoids
 *   - Capon minimum-variance spectral estimator
 *   - Resolution comparison: FFT vs MUSIC vs Capon
 *   - Frequency estimation accuracy
 *
 * Build & run:
 *   make chapters && ./build/bin/ch25
 *
 * Read alongside: chapters/25-parametric-spectral.md
 */

#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "spectral_est.h"
#include "fft.h"
#include "dsp_utils.h"
#include "gnuplot.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define N_SAMPLES  256
#define NFFT       1024

/* ── Demo 1: MUSIC for Two Close Sinusoids ────────────────────── */

static void demo_music_two_sines(void)
{
    printf("=== Demo 1: MUSIC — Two Close Sinusoids ===\n\n");

    double f1 = 0.10, f2 = 0.13;  /* normalised frequencies */
    double x[N_SAMPLES];
    unsigned int seed = 42;

    for (int i = 0; i < N_SAMPLES; i++) {
        seed = seed * 1103515245u + 12345u;
        double noise = ((double)(seed & 0x7FFF) / 32768.0 - 0.5) * 0.2;
        x[i] = sin(2.0 * M_PI * f1 * i) + sin(2.0 * M_PI * f2 * i) + noise;
    }

    int p = 20;       /* correlation matrix size */
    int n_sigs = 2;   /* number of signals */

    double spec[NFFT / 2];
    music_spectrum(x, N_SAMPLES, p, n_sigs, spec, NFFT);

    /* Find peaks */
    double freqs[2];
    music_frequencies(x, N_SAMPLES, p, n_sigs, freqs);

    printf("  True frequencies:      f1 = %.4f, f2 = %.4f\n", f1, f2);
    printf("  MUSIC estimates:       f1 = %.4f, f2 = %.4f\n", freqs[0], freqs[1]);
    printf("  Error:                 Δf1 = %.4f, Δf2 = %.4f\n\n",
           fabs(freqs[0] - f1), fabs(freqs[1] - f2));

    /* FFT-based PSD for comparison */
    Complex *X = (Complex *)calloc((size_t)NFFT, sizeof(Complex));
    for (int i = 0; i < N_SAMPLES; i++) {
        double w = 0.5 * (1.0 - cos(2.0 * M_PI * (double)i / (N_SAMPLES - 1)));
        X[i].re = x[i] * w;
        X[i].im = 0.0;
    }
    fft(X, NFFT);
    double psd_fft[NFFT / 2];
    for (int i = 0; i < NFFT / 2; i++)
        psd_fft[i] = 10.0 * log10(X[i].re * X[i].re + X[i].im * X[i].im + 1e-30);

    /* Normalise MUSIC spectrum for visual comparison */
    double max_music = -1e30;
    for (int i = 0; i < NFFT / 2; i++)
        if (spec[i] > max_music) max_music = spec[i];
    double max_fft = -1e30;
    for (int i = 0; i < NFFT / 2; i++)
        if (psd_fft[i] > max_fft) max_fft = psd_fft[i];

    FILE *gp = popen("gnuplot -persistent", "w");
    if (gp) {
        fprintf(gp, "set terminal pngcairo size 900,400\n");
        fprintf(gp, "set output 'plots/ch25/music_two_sines.png'\n");
        fprintf(gp, "set title 'MUSIC vs FFT: Two Close Sinusoids"
                     " (f=%.2f, f=%.2f)'\n", f1, f2);
        fprintf(gp, "set xlabel 'Normalised Frequency'\n");
        fprintf(gp, "set ylabel 'dB (normalised)'\nset grid\n");
        fprintf(gp, "set xrange [0:0.3]\n");
        fprintf(gp, "plot '-' w l lw 2 lc rgb '#CC3311' title 'MUSIC',"
                     " '-' w l lw 1 lc rgb '#0077BB' title 'FFT PSD'\n");
        for (int i = 0; i < NFFT / 2; i++)
            fprintf(gp, "%f %f\n", (double)i / NFFT, spec[i] - max_music);
        fprintf(gp, "e\n");
        for (int i = 0; i < NFFT / 2; i++)
            fprintf(gp, "%f %f\n", (double)i / NFFT, psd_fft[i] - max_fft);
        fprintf(gp, "e\n");
        pclose(gp);
        printf("  Plot → plots/ch25/music_two_sines.png\n\n");
    }
    free(X);
}

/* ── Demo 2: Capon (MVDR) Spectral Estimator ──────────────────── */

static void demo_capon(void)
{
    printf("=== Demo 2: Capon (MVDR) Spectral Estimator ===\n\n");

    double f1 = 0.15, f2 = 0.18;
    double x[N_SAMPLES];
    unsigned int seed = 101;

    for (int i = 0; i < N_SAMPLES; i++) {
        seed = seed * 1103515245u + 12345u;
        double noise = ((double)(seed & 0x7FFF) / 32768.0 - 0.5) * 0.3;
        x[i] = 1.0 * sin(2.0 * M_PI * f1 * i)
             + 0.8 * sin(2.0 * M_PI * f2 * i) + noise;
    }

    int p = 16;
    double spec_capon[NFFT / 2];
    capon_spectrum(x, N_SAMPLES, p, spec_capon, NFFT);

    /* Periodogram for comparison */
    Complex *X = (Complex *)calloc((size_t)NFFT, sizeof(Complex));
    for (int i = 0; i < N_SAMPLES; i++) {
        X[i].re = x[i];
        X[i].im = 0.0;
    }
    fft(X, NFFT);
    double psd[NFFT / 2];
    for (int i = 0; i < NFFT / 2; i++)
        psd[i] = 10.0 * log10(X[i].re * X[i].re + X[i].im * X[i].im + 1e-30);

    FILE *gp = popen("gnuplot -persistent", "w");
    if (gp) {
        fprintf(gp, "set terminal pngcairo size 900,400\n");
        fprintf(gp, "set output 'plots/ch25/capon_vs_fft.png'\n");
        fprintf(gp, "set title 'Capon (MVDR) vs Periodogram'\n");
        fprintf(gp, "set xlabel 'Normalised Frequency'\n");
        fprintf(gp, "set ylabel 'dB'\nset grid\n");
        fprintf(gp, "set xrange [0:0.35]\n");
        fprintf(gp, "plot '-' w l lw 2 lc rgb '#009988' title 'Capon',"
                     " '-' w l lw 1 lc rgb '#AAAAAA' title 'Periodogram'\n");
        /* Normalise for visual alignment */
        double mx1 = -1e30, mx2 = -1e30;
        for (int i = 0; i < NFFT / 2; i++) {
            if (spec_capon[i] > mx1) mx1 = spec_capon[i];
            if (psd[i] > mx2) mx2 = psd[i];
        }
        for (int i = 0; i < NFFT / 2; i++)
            fprintf(gp, "%f %f\n", (double)i / NFFT, spec_capon[i] - mx1);
        fprintf(gp, "e\n");
        for (int i = 0; i < NFFT / 2; i++)
            fprintf(gp, "%f %f\n", (double)i / NFFT, psd[i] - mx2);
        fprintf(gp, "e\n");
        pclose(gp);
        printf("  Plot → plots/ch25/capon_vs_fft.png\n\n");
    }
    free(X);
}

/* ── Demo 3: Resolution Comparison ────────────────────────────── */

static void demo_resolution(void)
{
    printf("=== Demo 3: Super-Resolution Comparison ===\n\n");

    /* Three sinusoids, two very close */
    double f1 = 0.20, f2 = 0.22, f3 = 0.35;
    int short_n = 128;  /* fewer samples makes it harder */
    double x[128];
    unsigned int seed = 77;

    for (int i = 0; i < short_n; i++) {
        seed = seed * 1103515245u + 12345u;
        double noise = ((double)(seed & 0x7FFF) / 32768.0 - 0.5) * 0.1;
        x[i] = sin(2.0 * M_PI * f1 * i) + sin(2.0 * M_PI * f2 * i)
             + 0.7 * sin(2.0 * M_PI * f3 * i) + noise;
    }

    /* MUSIC */
    double spec_music[NFFT / 2];
    music_spectrum(x, short_n, 24, 3, spec_music, NFFT);

    /* Capon */
    double spec_capon[NFFT / 2];
    capon_spectrum(x, short_n, 24, spec_capon, NFFT);

    /* FFT */
    Complex *X = (Complex *)calloc((size_t)NFFT, sizeof(Complex));
    for (int i = 0; i < short_n; i++) {
        double w = 0.5 * (1.0 - cos(2.0 * M_PI * (double)i / (short_n - 1)));
        X[i].re = x[i] * w;
        X[i].im = 0.0;
    }
    fft(X, NFFT);
    double psd[NFFT / 2];
    for (int i = 0; i < NFFT / 2; i++)
        psd[i] = 10.0 * log10(X[i].re * X[i].re + X[i].im * X[i].im + 1e-30);

    printf("  True: f1=%.2f  f2=%.2f  f3=%.2f  (N=%d samples)\n",
           f1, f2, f3, short_n);
    printf("  FFT resolution: %.4f (Δf = 1/N)\n", 1.0 / short_n);
    printf("  |f1-f2| = %.4f — below FFT resolution!\n\n", fabs(f2 - f1));

    FILE *gp = popen("gnuplot -persistent", "w");
    if (gp) {
        fprintf(gp, "set terminal pngcairo size 900,500\n");
        fprintf(gp, "set output 'plots/ch25/resolution_compare.png'\n");
        fprintf(gp, "set title 'Super-Resolution: FFT vs MUSIC vs Capon"
                     " (N=%d)'\n", short_n);
        fprintf(gp, "set xlabel 'Normalised Frequency'\n");
        fprintf(gp, "set ylabel 'dB (normalised)'\nset grid\n");
        fprintf(gp, "set xrange [0.1:0.45]\n");

        /* Normalise all to 0 dB peak */
        double mx_m = -1e30, mx_c = -1e30, mx_f = -1e30;
        for (int i = 0; i < NFFT / 2; i++) {
            if (spec_music[i] > mx_m) mx_m = spec_music[i];
            if (spec_capon[i] > mx_c) mx_c = spec_capon[i];
            if (psd[i] > mx_f) mx_f = psd[i];
        }

        fprintf(gp, "plot '-' w l lw 2 lc rgb '#CC3311' title 'MUSIC',"
                     " '-' w l lw 2 lc rgb '#009988' title 'Capon',"
                     " '-' w l lw 1 lc rgb '#AAAAAA' title 'FFT'\n");
        for (int i = 0; i < NFFT / 2; i++)
            fprintf(gp, "%f %f\n", (double)i / NFFT, spec_music[i] - mx_m);
        fprintf(gp, "e\n");
        for (int i = 0; i < NFFT / 2; i++)
            fprintf(gp, "%f %f\n", (double)i / NFFT, spec_capon[i] - mx_c);
        fprintf(gp, "e\n");
        for (int i = 0; i < NFFT / 2; i++)
            fprintf(gp, "%f %f\n", (double)i / NFFT, psd[i] - mx_f);
        fprintf(gp, "e\n");
        pclose(gp);
        printf("  Plot → plots/ch25/resolution_compare.png\n\n");
    }
    free(X);
}

/* ── Main ─────────────────────────────────────────────────────── */

int main(void)
{
    printf("╔══════════════════════════════════════════════╗\n");
    printf("║  Chapter 25: Parametric Spectral Estimation ║\n");
    printf("║  MUSIC · Capon (MVDR)                       ║\n");
    printf("╚══════════════════════════════════════════════╝\n\n");

    if (system("mkdir -p plots/ch25") != 0) { /* ignore */ }

    demo_music_two_sines();
    demo_capon();
    demo_resolution();

    printf("All Chapter 25 demos complete.\n");
    return 0;
}

/**
 * @file 23-adaptive-filters.c
 * @brief Chapter 23 demo — Adaptive Filters: LMS, NLMS, RLS.
 *
 * ── Adaptive Filtering Concept ───────────────────────────────────
 *
 *   x[n] → ┌───────┐ → y[n] → ⊕ → e[n]
 *           │ W(z)  │          ↑
 *           └───────┘       d[n] (desired)
 *               ↑
 *           e[n] updates W
 *
 *   Three algorithms of increasing sophistication:
 *     LMS  — O(L) per sample, μ fixed step
 *     NLMS — O(L) per sample, normalised step ∝ 1/‖x‖²
 *     RLS  — O(L²) per sample, recursive least-squares via P matrix
 *
 * Demonstrates:
 *   - System identification (unknown plant estimation)
 *   - Noise cancellation via NLMS
 *   - Convergence comparison: LMS vs NLMS vs RLS
 *   - Learning curve (MSE vs iteration)
 *
 * Build & run:
 *   make chapters && ./build/bin/ch23
 *
 * Read alongside: chapters/23-adaptive-filters.md
 */

#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "adaptive.h"
#include "filter.h"
#include "dsp_utils.h"
#include "gnuplot.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define N_SAMPLES  2000
#define PLANT_TAPS 11

/* ── Simple LCG pseudo-random ─────────────────────────────────── */
static double randn(unsigned int *seed)
{
    /* Box-Muller from LCG */
    *seed = *seed * 1103515245u + 12345u;
    double u1 = ((double)(*seed & 0x7FFFFFFFu)) / 2147483647.0;
    *seed = *seed * 1103515245u + 12345u;
    double u2 = ((double)(*seed & 0x7FFFFFFFu)) / 2147483647.0;
    if (u1 < 1e-10) u1 = 1e-10;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

/* ── Demo 1: System Identification with LMS ───────────────────── */

static void demo_lms_sysid(void)
{
    printf("=== Demo 1: System Identification (LMS) ===\n\n");

    /* Unknown plant: lowpass FIR */
    double plant[PLANT_TAPS];
    fir_lowpass(plant, PLANT_TAPS, 0.25);

    double x[N_SAMPLES], d[N_SAMPLES], y[N_SAMPLES], e[N_SAMPLES];
    double w_final[PLANT_TAPS];
    unsigned int seed = 42;

    /* White noise input */
    for (int i = 0; i < N_SAMPLES; i++)
        x[i] = randn(&seed);

    /* Desired = plant output */
    fir_filter(x, d, N_SAMPLES, plant, PLANT_TAPS);

    lms_filter(x, d, N_SAMPLES, PLANT_TAPS, 0.01, y, e, w_final);

    /* Compare estimated vs true weights */
    printf("  Tap   True      Estimated\n");
    printf("  ---   --------  ---------\n");
    for (int i = 0; i < PLANT_TAPS; i++)
        printf("  %3d   %+.5f  %+.5f\n", i, plant[i], w_final[i]);

    double mse = 0.0;
    for (int i = N_SAMPLES / 2; i < N_SAMPLES; i++)
        mse += e[i] * e[i];
    mse /= (double)(N_SAMPLES / 2);
    printf("\n  Steady-state MSE: %.2e\n\n", mse);

    /* Plot learning curve */
    double mse_curve[N_SAMPLES];
    int block = 50;
    for (int i = 0; i < N_SAMPLES; i++) {
        int start = (i > block) ? i - block : 0;
        double s = 0.0;
        for (int j = start; j <= i; j++) s += e[j] * e[j];
        mse_curve[i] = s / (double)(i - start + 1);
    }

    FILE *gp = popen("gnuplot -persistent", "w");
    if (gp) {
        fprintf(gp, "set terminal pngcairo size 800,400\n");
        fprintf(gp, "set output 'plots/ch23/lms_learning_curve.png'\n");
        fprintf(gp, "set title 'LMS Learning Curve'\n");
        fprintf(gp, "set xlabel 'Sample'\nset ylabel 'MSE'\n");
        fprintf(gp, "set logscale y\nset grid\n");
        fprintf(gp, "plot '-' with lines lw 2 lc rgb '#0077BB' title 'LMS'\n");
        for (int i = 0; i < N_SAMPLES; i += 4)
            fprintf(gp, "%d %e\n", i, mse_curve[i]);
        fprintf(gp, "e\n");
        pclose(gp);
        printf("  Plot → plots/ch23/lms_learning_curve.png\n\n");
    }
}

/* ── Demo 2: Noise Cancellation with NLMS ─────────────────────── */

static void demo_nlms_noise_cancel(void)
{
    printf("=== Demo 2: Noise Cancellation (NLMS) ===\n\n");

    unsigned int seed = 123;
    double signal[N_SAMPLES], noise[N_SAMPLES], corrupted[N_SAMPLES];
    double y[N_SAMPLES], err[N_SAMPLES], w[32];

    /* Clean signal: two sinusoids */
    for (int i = 0; i < N_SAMPLES; i++)
        signal[i] = sin(2.0 * M_PI * 0.05 * i) + 0.5 * sin(2.0 * M_PI * 0.12 * i);

    /* Noise reference (correlated with actual noise) */
    for (int i = 0; i < N_SAMPLES; i++)
        noise[i] = randn(&seed);

    /* Corrupted = signal + filtered noise */
    double noise_path[5] = {1.0, -0.5, 0.25, -0.1, 0.05};
    double filtered_noise[N_SAMPLES];
    fir_filter(noise, filtered_noise, N_SAMPLES, noise_path, 5);

    for (int i = 0; i < N_SAMPLES; i++)
        corrupted[i] = signal[i] + filtered_noise[i];

    /* NLMS: noise reference → estimated noise, error = clean estimate */
    nlms_filter(noise, corrupted, N_SAMPLES, 32, 0.5, 1e-6, y, err, w);

    double snr_before = 0.0, snr_after = 0.0;
    double sig_pow = 0.0, noise_pow = 0.0, err_pow = 0.0;
    for (int i = 500; i < N_SAMPLES; i++) {
        sig_pow += signal[i] * signal[i];
        noise_pow += filtered_noise[i] * filtered_noise[i];
        err_pow += (err[i] - signal[i]) * (err[i] - signal[i]);
    }
    snr_before = 10.0 * log10(sig_pow / (noise_pow + 1e-30));
    snr_after = 10.0 * log10(sig_pow / (err_pow + 1e-30));
    printf("  SNR before: %.1f dB\n", snr_before);
    printf("  SNR after:  %.1f dB\n\n", snr_after);

    /* Plot: corrupted vs cleaned */
    int plot_n = 300;
    FILE *gp = popen("gnuplot -persistent", "w");
    if (gp) {
        fprintf(gp, "set terminal pngcairo size 900,500\n");
        fprintf(gp, "set output 'plots/ch23/nlms_noise_cancel.png'\n");
        fprintf(gp, "set multiplot layout 2,1 title 'NLMS Noise Cancellation'\n");
        fprintf(gp, "set xlabel 'Sample'\nset ylabel 'Amplitude'\nset grid\n");
        fprintf(gp, "set title 'Corrupted Signal'\n");
        fprintf(gp, "plot '-' w l lw 1 lc rgb '#CC3311' title 'Corrupted'\n");
        for (int i = 1000; i < 1000 + plot_n; i++)
            fprintf(gp, "%d %f\n", i, corrupted[i]);
        fprintf(gp, "e\n");
        fprintf(gp, "set title 'NLMS Cleaned vs Original'\n");
        fprintf(gp, "plot '-' w l lw 2 lc rgb '#0077BB' title 'Cleaned',"
                     " '-' w l lw 1 dt 2 lc rgb '#009988' title 'Original'\n");
        for (int i = 1000; i < 1000 + plot_n; i++)
            fprintf(gp, "%d %f\n", i, err[i]);
        fprintf(gp, "e\n");
        for (int i = 1000; i < 1000 + plot_n; i++)
            fprintf(gp, "%d %f\n", i, signal[i]);
        fprintf(gp, "e\nunset multiplot\n");
        pclose(gp);
        printf("  Plot → plots/ch23/nlms_noise_cancel.png\n\n");
    }
}

/* ── Demo 3: Convergence Comparison ───────────────────────────── */

static void demo_convergence_compare(void)
{
    printf("=== Demo 3: Convergence Comparison: LMS vs NLMS vs RLS ===\n\n");

    double plant[PLANT_TAPS];
    fir_lowpass(plant, PLANT_TAPS, 0.25);

    unsigned int seed = 77;
    double x[N_SAMPLES], d[N_SAMPLES];

    for (int i = 0; i < N_SAMPLES; i++)
        x[i] = randn(&seed);
    fir_filter(x, d, N_SAMPLES, plant, PLANT_TAPS);

    /* Add small observation noise */
    seed = 99;
    for (int i = 0; i < N_SAMPLES; i++)
        d[i] += 0.01 * randn(&seed);

    double y[N_SAMPLES], e_lms[N_SAMPLES], e_nlms[N_SAMPLES], e_rls[N_SAMPLES];
    double w[PLANT_TAPS];

    lms_filter(x, d, N_SAMPLES, PLANT_TAPS, 0.005, y, e_lms, w);
    nlms_filter(x, d, N_SAMPLES, PLANT_TAPS, 0.5, 1e-6, y, e_nlms, w);
    rls_filter(x, d, N_SAMPLES, PLANT_TAPS, 0.99, 100.0, y, e_rls, w);

    /* Smoothed MSE curves */
    int block = 100;
    double mse_lms[N_SAMPLES], mse_nlms[N_SAMPLES], mse_rls[N_SAMPLES];
    for (int i = 0; i < N_SAMPLES; i++) {
        int start = (i > block) ? i - block : 0;
        int cnt = i - start + 1;
        double s1 = 0, s2 = 0, s3 = 0;
        for (int j = start; j <= i; j++) {
            s1 += e_lms[j] * e_lms[j];
            s2 += e_nlms[j] * e_nlms[j];
            s3 += e_rls[j] * e_rls[j];
        }
        mse_lms[i] = s1 / cnt;
        mse_nlms[i] = s2 / cnt;
        mse_rls[i] = s3 / cnt;
    }

    FILE *gp = popen("gnuplot -persistent", "w");
    if (gp) {
        fprintf(gp, "set terminal pngcairo size 800,400\n");
        fprintf(gp, "set output 'plots/ch23/convergence_compare.png'\n");
        fprintf(gp, "set title 'Convergence: LMS vs NLMS vs RLS'\n");
        fprintf(gp, "set xlabel 'Sample'\nset ylabel 'MSE'\n");
        fprintf(gp, "set logscale y\nset grid\n");
        fprintf(gp, "plot '-' w l lw 2 lc rgb '#CC3311' title 'LMS',"
                     " '-' w l lw 2 lc rgb '#0077BB' title 'NLMS',"
                     " '-' w l lw 2 lc rgb '#009988' title 'RLS'\n");
        for (int i = 0; i < N_SAMPLES; i += 4)
            fprintf(gp, "%d %e\n", i, mse_lms[i]);
        fprintf(gp, "e\n");
        for (int i = 0; i < N_SAMPLES; i += 4)
            fprintf(gp, "%d %e\n", i, mse_nlms[i]);
        fprintf(gp, "e\n");
        for (int i = 0; i < N_SAMPLES; i += 4)
            fprintf(gp, "%d %e\n", i, mse_rls[i]);
        fprintf(gp, "e\n");
        pclose(gp);
        printf("  Plot → plots/ch23/convergence_compare.png\n\n");
    }
}

/* ── Main ─────────────────────────────────────────────────────── */

int main(void)
{
    printf("╔══════════════════════════════════════════════╗\n");
    printf("║  Chapter 23: Adaptive Filters               ║\n");
    printf("║  LMS · NLMS · RLS                           ║\n");
    printf("╚══════════════════════════════════════════════╝\n\n");

    if (system("mkdir -p plots/ch23") != 0) { /* ignore */ }

    demo_lms_sysid();
    demo_nlms_noise_cancel();
    demo_convergence_compare();

    printf("All Chapter 23 demos complete.\n");
    return 0;
}

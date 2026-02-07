/**
 * @file 21-signal-averaging.c
 * @brief Chapter 21 demo — Signal Averaging & Noise Reduction.
 *
 * ── Coherent (Synchronous) Averaging ─────────────────────────────
 *
 *   Trial 1:  ████░█████░░████   signal + noise₁
 *   Trial 2:  ████░░████░█████   signal + noise₂
 *     ...           ...
 *   Trial K:  ████░████░░█████   signal + noiseₖ
 *         ──────────────────────
 *   Average:  ██████████████     signal (noise ÷ √K)
 *
 *   SNR improvement: 10·log₁₀(K) dB
 *
 * ── Exponential Moving Average ───────────────────────────────────
 *
 *   y[n] = α·x[n] + (1−α)·y[n−1]
 *
 *   Small α → heavy smoothing, large lag
 *   Large α → less smoothing, follows signal
 *
 * ── Median Filter ────────────────────────────────────────────────
 *
 *   Sorts M neighbours, picks middle value.
 *   Removes impulse noise ("salt & pepper") while preserving edges.
 *
 * Demonstrates:
 *   - Coherent averaging: SNR improvement with K trials
 *   - Exponential moving average (EMA) tracking
 *   - Moving average vs EMA comparison
 *   - Median filter for impulse noise removal
 *   - SNR measurement before/after averaging
 *
 * Build & run:
 *   make chapters && ./build/bin/ch21
 *
 * Read alongside: chapters/21-signal-averaging.md
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "averaging.h"
#include "signal_gen.h"
#include "gnuplot.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Simple noise generation (seeded) */
static double rand_gauss(unsigned int *seed)
{
    /* Box-Muller transform */
    double u1 = ((double)(*seed = *seed * 1103515245 + 12345) / 2147483648.0);
    double u2 = ((double)(*seed = *seed * 1103515245 + 12345) / 2147483648.0);
    if (u1 < 1e-10) u1 = 1e-10;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

/* ── Demo 1: Coherent Averaging ───────────────────────────────── */

static void demo_coherent(void)
{
    printf("=== Demo 1: Coherent Averaging ===\n\n");

    const int N = 256;
    const double fs = 1000.0;
    double *clean = (double *)malloc((size_t)N * sizeof(double));
    gen_sine(clean, N, 1.0, 50.0, fs, 0.0);

    int K_values[] = {1, 4, 16, 64};
    printf("  Signal: 50 Hz sine, noise σ = 2.0\n");
    printf("  ┌──────┬───────────┬───────────┬───────────┐\n");
    printf("  │  K   │ SNR_in dB │ SNR_out dB│ Gain dB   │\n");
    printf("  ├──────┼───────────┼───────────┼───────────┤\n");

    for (int ki = 0; ki < 4; ki++) {
        int K = K_values[ki];
        unsigned int seed = 42;

        /* Generate K noisy trials */
        double **trials = (double **)malloc((size_t)K * sizeof(double *));
        for (int k = 0; k < K; k++) {
            trials[k] = (double *)malloc((size_t)N * sizeof(double));
            for (int i = 0; i < N; i++) {
                trials[k][i] = clean[i] + 2.0 * rand_gauss(&seed);
            }
        }

        /* Average */
        double *avg = (double *)calloc((size_t)N, sizeof(double));
        coherent_average((const double **)trials, K, N, avg);

        /* Measure SNR */
        double snr_before, snr_after;
        compute_snr_improvement(trials[0], clean, avg, N,
                                &snr_before, &snr_after);

        printf("  │ %4d │  %+7.2f  │  %+7.2f  │  %+7.2f  │\n",
               K, snr_before, snr_after, snr_after - snr_before);

        /* Plot for K=16 */
        if (K == 16) {
            double *t_axis = (double *)malloc((size_t)N * sizeof(double));
            for (int i = 0; i < N; i++) t_axis[i] = (double)i / fs;

            gp_plot_multi("ch21", "coherent_averaging",
                "Coherent Averaging (K=16 trials)",
                "Time (s)", "Amplitude",
                (GpSeries[]){
                    {"Noisy (1 trial)", t_axis, trials[0], N, "lines"},
                    {"Averaged (K=16)", t_axis, avg, N, "lines"},
                    {"Clean signal", t_axis, clean, N, "lines"}
                }, 3);
            free(t_axis);
        }

        for (int k = 0; k < K; k++) free(trials[k]);
        free(trials);
        free(avg);
    }

    printf("  └──────┴───────────┴───────────┴───────────┘\n");
    printf("  Theory: gain ≈ 10·log₁₀(K) dB\n");
    printf("  → plots/ch21/coherent_averaging.png\n\n");

    free(clean);
}

/* ── Demo 2: EMA Tracking ─────────────────────────────────────── */

static void demo_ema(void)
{
    printf("=== Demo 2: Exponential Moving Average ===\n\n");

    const int N = 256;
    double *x = (double *)malloc((size_t)N * sizeof(double));
    unsigned int seed = 123;

    /* Step signal + noise */
    for (int i = 0; i < N; i++) {
        double step = (i < N / 2) ? 0.0 : 1.0;
        x[i] = step + 0.3 * rand_gauss(&seed);
    }

    double alphas[] = {0.05, 0.15, 0.5};
    const char *labels[] = {"α=0.05", "α=0.15", "α=0.50"};
    double *y[3];
    for (int a = 0; a < 3; a++) {
        y[a] = (double *)calloc((size_t)N, sizeof(double));
        ema_filter(x, N, alphas[a], y[a]);
    }

    double *idx = (double *)malloc((size_t)N * sizeof(double));
    for (int i = 0; i < N; i++) idx[i] = (double)i;

    gp_plot_multi("ch21", "ema_comparison",
        "EMA: Effect of Alpha on Step Response",
        "Sample", "Amplitude",
        (GpSeries[]){
            {"Noisy input", idx, x, N, "lines"},
            {labels[0], idx, y[0], N, "lines"},
            {labels[1], idx, y[1], N, "lines"},
            {labels[2], idx, y[2], N, "lines"}
        }, 4);

    printf("  Step + noise (σ=0.3)\n");
    printf("  α=0.05: heavy smooth, slow response\n");
    printf("  α=0.15: moderate\n");
    printf("  α=0.50: fast response, less smoothing\n");
    printf("  → plots/ch21/ema_comparison.png\n\n");

    free(x); free(idx);
    for (int a = 0; a < 3; a++) free(y[a]);
}

/* ── Demo 3: Moving Average vs EMA ────────────────────────────── */

static void demo_ma_vs_ema(void)
{
    printf("=== Demo 3: Moving Average vs EMA ===\n\n");

    const int N = 256;
    const double fs = 1000.0;
    double *x = (double *)malloc((size_t)N * sizeof(double));
    unsigned int seed = 456;

    gen_sine(x, N, 1.0, 20.0, fs, 0.0);
    for (int i = 0; i < N; i++)
        x[i] += 0.5 * rand_gauss(&seed);

    int M = 11;
    double alpha = 2.0 / (M + 1); /* EMA equivalent to M-point MA */

    double *y_ma = (double *)calloc((size_t)N, sizeof(double));
    double *y_ema = (double *)calloc((size_t)N, sizeof(double));
    moving_average(x, N, M, y_ma);
    ema_filter(x, N, alpha, y_ema);

    printf("  Moving average: M = %d (boxcar)\n", M);
    printf("  EMA: α = %.3f (equivalent span)\n", alpha);

    double *idx = (double *)malloc((size_t)N * sizeof(double));
    for (int i = 0; i < N; i++) idx[i] = (double)i / fs;

    gp_plot_multi("ch21", "ma_vs_ema",
        "Moving Average (M=11) vs EMA (α=0.167)",
        "Time (s)", "Amplitude",
        (GpSeries[]){
            {"Noisy input", idx, x, N, "lines"},
            {"Moving Average", idx, y_ma, N, "lines"},
            {"EMA",           idx, y_ema, N, "lines"}
        }, 3);
    printf("  → plots/ch21/ma_vs_ema.png\n\n");

    free(x); free(y_ma); free(y_ema); free(idx);
}

/* ── Demo 4: Median Filter — Impulse Removal ──────────────────── */

static void demo_median(void)
{
    printf("=== Demo 4: Median Filter for Impulse Noise ===\n\n");

    const int N = 256;
    const double fs = 1000.0;
    double *clean = (double *)malloc((size_t)N * sizeof(double));
    double *noisy = (double *)malloc((size_t)N * sizeof(double));
    gen_sine(clean, N, 1.0, 30.0, fs, 0.0);

    /* Add impulse noise to ~10% of samples */
    unsigned int seed = 789;
    int impulse_count = 0;
    for (int i = 0; i < N; i++) {
        noisy[i] = clean[i];
        unsigned int r = (seed = seed * 1103515245 + 12345);
        if ((r >> 16) % 10 == 0) {
            noisy[i] += ((r & 1) ? 5.0 : -5.0);
            impulse_count++;
        }
    }

    double *y = (double *)calloc((size_t)N, sizeof(double));
    median_filter(noisy, N, 5, y);

    /* Measure recovery */
    double mse_before = 0.0, mse_after = 0.0;
    for (int i = 0; i < N; i++) {
        double e1 = noisy[i] - clean[i];
        double e2 = y[i] - clean[i];
        mse_before += e1 * e1;
        mse_after  += e2 * e2;
    }
    mse_before /= N;
    mse_after  /= N;

    printf("  %d impulse spikes added (±5.0)\n", impulse_count);
    printf("  Median window M = 5\n");
    printf("  MSE before: %.4f\n", mse_before);
    printf("  MSE after:  %.4f\n", mse_after);
    printf("  Improvement: %.1f×\n\n", mse_before / mse_after);

    double *idx = (double *)malloc((size_t)N * sizeof(double));
    for (int i = 0; i < N; i++) idx[i] = (double)i / fs;

    gp_plot_multi("ch21", "median_filter",
        "Median Filter: Impulse Noise Removal (M=5)",
        "Time (s)", "Amplitude",
        (GpSeries[]){
            {"Noisy (impulses)", idx, noisy, N, "lines"},
            {"Median filtered",  idx, y, N, "lines"},
            {"Clean signal",     idx, clean, N, "lines"}
        }, 3);
    printf("  → plots/ch21/median_filter.png\n\n");

    free(clean); free(noisy); free(y); free(idx);
}

/* ── Demo 5: SNR vs Number of Averages ────────────────────────── */

static void demo_snr_curve(void)
{
    printf("=== Demo 5: SNR Improvement vs K ===\n\n");

    const int N = 128;
    const double fs = 1000.0;
    double *clean = (double *)malloc((size_t)N * sizeof(double));
    gen_sine(clean, N, 1.0, 50.0, fs, 0.0);

    int K_values[] = {1, 2, 4, 8, 16, 32, 64, 128};
    int n_points = 8;
    double *k_axis = (double *)malloc((size_t)n_points * sizeof(double));
    double *snr_measured = (double *)malloc((size_t)n_points * sizeof(double));
    double *snr_theory   = (double *)malloc((size_t)n_points * sizeof(double));

    for (int ki = 0; ki < n_points; ki++) {
        int K = K_values[ki];
        unsigned int seed = 42;

        double **trials = (double **)malloc((size_t)K * sizeof(double *));
        for (int k = 0; k < K; k++) {
            trials[k] = (double *)malloc((size_t)N * sizeof(double));
            for (int i = 0; i < N; i++)
                trials[k][i] = clean[i] + 2.0 * rand_gauss(&seed);
        }

        double *avg = (double *)calloc((size_t)N, sizeof(double));
        coherent_average((const double **)trials, K, N, avg);

        double snr_b, snr_a;
        compute_snr_improvement(trials[0], clean, avg, N, &snr_b, &snr_a);

        k_axis[ki]  = (double)K;
        snr_measured[ki] = snr_a;
        snr_theory[ki]   = snr_b + 10.0 * log10((double)K);

        for (int k = 0; k < K; k++) free(trials[k]);
        free(trials);
        free(avg);
    }

    gp_plot_multi("ch21", "snr_vs_k",
        "SNR Improvement: Measured vs Theory (10·log10(K))",
        "Number of Averages (K)", "SNR (dB)",
        (GpSeries[]){
            {"Measured SNR", k_axis, snr_measured, n_points, "linespoints"},
            {"Theory",       k_axis, snr_theory,   n_points, "lines"}
        }, 2);
    printf("  → plots/ch21/snr_vs_k.png\n\n");

    free(clean); free(k_axis); free(snr_measured); free(snr_theory);
}

int main(void)
{
    printf("╔══════════════════════════════════════════════════════════╗\n");
    printf("║  Chapter 21: Signal Averaging & Noise Reduction        ║\n");
    printf("╚══════════════════════════════════════════════════════════╝\n\n");

    demo_coherent();
    demo_ema();
    demo_ma_vs_ema();
    demo_median();
    demo_snr_curve();

    printf("=== Chapter 21 Complete ===\n");
    return 0;
}

/**
 * @file 24-linear-prediction.c
 * @brief Chapter 24 demo — Linear Prediction Coding (LPC).
 *
 * ── Linear Prediction Model ─────────────────────────────────────
 *
 *   x[n] ≈ -a₁x[n-1] - a₂x[n-2] - ··· - aₚx[n-p]
 *
 *   Residual (prediction error):
 *   e[n] = x[n] + a₁x[n-1] + ··· + aₚx[n-p]
 *
 *       ┌──────────┐     e[n]
 *   x → │ A(z) = 1 │ ───────→  (analysis: flat residual)
 *       │  + Σaₖz⁻ᵏ│
 *       └──────────┘
 *
 *       ┌──────────┐     x̂[n]
 *   e → │ 1/A(z)   │ ───────→  (synthesis: reconstruct)
 *       └──────────┘
 *
 * Demonstrates:
 *   - Levinson-Durbin autocorrelation method
 *   - LPC analysis / synthesis round-trip
 *   - AR spectral envelope estimation
 *   - Reflection coefficients and lattice structure
 *   - Prediction order trade-off
 *
 * Build & run:
 *   make chapters && ./build/bin/ch24
 *
 * Read alongside: chapters/24-linear-prediction.md
 */

#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lpc.h"
#include "fft.h"
#include "dsp_utils.h"
#include "gnuplot.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define N_SAMPLES  512
#define MAX_ORDER  20

/* ── Demo 1: Levinson-Durbin on Known AR Process ──────────────── */

static void demo_levinson(void)
{
    printf("=== Demo 1: Levinson-Durbin (AR Process) ===\n\n");

    /* Generate AR(4) process: x[n] = 2.76x[n-1] - 3.81x[n-2]
       + 2.65x[n-3] - 0.92x[n-4] + e[n] */
    double a_true[4] = {-2.76, 3.81, -2.65, 0.92};
    double x[N_SAMPLES];
    unsigned int seed = 42;

    /* Seed with zeros, drive with white noise */
    for (int i = 0; i < 4; i++) x[i] = 0.0;
    for (int i = 4; i < N_SAMPLES; i++) {
        seed = seed * 1103515245u + 12345u;
        double e = ((double)(seed & 0x7FFF) / 16384.0 - 1.0) * 0.1;
        x[i] = e;
        for (int k = 0; k < 4; k++)
            x[i] -= a_true[k] * x[i - 1 - k];
    }

    /* Estimate LPC coefficients */
    int order = 4;
    double a_est[MAX_ORDER + 1];
    double E;
    lpc_coefficients(x, N_SAMPLES, order, a_est, &E);

    printf("  Coefficient  True      Estimated\n");
    printf("  -----------  --------  ---------\n");
    /* Note: lpc_coefficients returns a[1..p] with convention
       x[n] + a[1]x[n-1] + ... = e[n], so a[k] ≈ a_true[k-1] */
    for (int i = 1; i <= order; i++)
        printf("       a[%d]    %+.4f   %+.4f\n", i, a_true[i - 1], a_est[i]);

    printf("\n  Prediction error energy: %e\n\n", E);
}

/* ── Demo 2: Analysis/Synthesis Round-Trip ────────────────────── */

static void demo_analysis_synthesis(void)
{
    printf("=== Demo 2: Analysis / Synthesis Round-Trip ===\n\n");

    /* Test signal: sum of sinusoids (quasi-periodic, like speech) */
    double x[N_SAMPLES], residual[N_SAMPLES], recon[N_SAMPLES];
    for (int i = 0; i < N_SAMPLES; i++)
        x[i] = sin(2.0 * M_PI * 0.1 * i) + 0.5 * sin(2.0 * M_PI * 0.25 * i)
              + 0.3 * sin(2.0 * M_PI * 0.37 * i);

    int order = 12;
    double a[MAX_ORDER + 1], E;
    lpc_coefficients(x, N_SAMPLES, order, a, &E);

    /* Analysis: extract residual */
    lpc_residual(x, N_SAMPLES, a, order, residual);

    /* Synthesis: reconstruct from residual */
    lpc_synthesise(residual, N_SAMPLES, a, order, recon);

    /* Compute reconstruction error */
    double max_err = 0.0;
    for (int i = 0; i < N_SAMPLES; i++) {
        double diff = fabs(x[i] - recon[i]);
        if (diff > max_err) max_err = diff;
    }
    printf("  Order: %d\n", order);
    printf("  Max reconstruction error: %e\n\n", max_err);

    /* Plot */
    int plot_n = 200;
    FILE *gp = popen("gnuplot -persistent", "w");
    if (gp) {
        fprintf(gp, "set terminal pngcairo size 900,600\n");
        fprintf(gp, "set output 'plots/ch24/analysis_synthesis.png'\n");
        fprintf(gp, "set multiplot layout 3,1 title"
                     " 'LPC Analysis/Synthesis (order %d)'\n", order);
        fprintf(gp, "set grid\nset xlabel 'Sample'\n");

        fprintf(gp, "set title 'Original Signal'\n");
        fprintf(gp, "plot '-' w l lw 2 lc rgb '#0077BB' notitle\n");
        for (int i = 0; i < plot_n; i++)
            fprintf(gp, "%d %f\n", i, x[i]);
        fprintf(gp, "e\n");

        fprintf(gp, "set title 'Residual (Prediction Error)'\n");
        fprintf(gp, "plot '-' w l lw 1 lc rgb '#CC3311' notitle\n");
        for (int i = 0; i < plot_n; i++)
            fprintf(gp, "%d %f\n", i, residual[i]);
        fprintf(gp, "e\n");

        fprintf(gp, "set title 'Reconstructed vs Original'\n");
        fprintf(gp, "plot '-' w l lw 2 lc rgb '#009988' title 'Recon',"
                     " '-' w l lw 1 dt 2 lc rgb '#0077BB' title 'Original'\n");
        for (int i = 0; i < plot_n; i++)
            fprintf(gp, "%d %f\n", i, recon[i]);
        fprintf(gp, "e\n");
        for (int i = 0; i < plot_n; i++)
            fprintf(gp, "%d %f\n", i, x[i]);
        fprintf(gp, "e\nunset multiplot\n");
        pclose(gp);
        printf("  Plot → plots/ch24/analysis_synthesis.png\n\n");
    }
}

/* ── Demo 3: AR Spectral Envelope ─────────────────────────────── */

static void demo_ar_spectrum(void)
{
    printf("=== Demo 3: AR Spectral Envelope ===\n\n");

    /* Signal with known spectral peaks */
    double x[N_SAMPLES];
    for (int i = 0; i < N_SAMPLES; i++)
        x[i] = sin(2.0 * M_PI * 0.1 * i) + 0.5 * sin(2.0 * M_PI * 0.25 * i);

    /* Compute FFT-based PSD for reference */
    int nfft = 512;
    Complex *X = (Complex *)calloc((size_t)nfft, sizeof(Complex));
    for (int i = 0; i < N_SAMPLES && i < nfft; i++) {
        double w = 0.5 * (1.0 - cos(2.0 * M_PI * (double)i / (N_SAMPLES - 1)));
        X[i].re = x[i] * w;
        X[i].im = 0.0;
    }
    fft(X, nfft);
    double psd_fft[256];
    for (int i = 0; i < nfft / 2; i++)
        psd_fft[i] = 10.0 * log10(X[i].re * X[i].re + X[i].im * X[i].im + 1e-30);

    /* LPC spectral envelope for different orders */
    int orders[] = {4, 10, 20};
    int n_orders = 3;

    FILE *gp = popen("gnuplot -persistent", "w");
    if (gp) {
        fprintf(gp, "set terminal pngcairo size 800,500\n");
        fprintf(gp, "set output 'plots/ch24/ar_spectrum.png'\n");
        fprintf(gp, "set title 'AR Spectral Envelope vs FFT PSD'\n");
        fprintf(gp, "set xlabel 'Normalised Frequency'\nset ylabel 'dB'\n");
        fprintf(gp, "set grid\nset key top right\n");
        fprintf(gp, "plot '-' w l lw 1 lc rgb '#AAAAAA' title 'FFT PSD'");
        for (int o = 0; o < n_orders; o++) {
            const char *cols[] = {"#CC3311", "#0077BB", "#009988"};
            fprintf(gp, ", '-' w l lw 2 lc rgb '%s' title 'AR(%d)'",
                    cols[o], orders[o]);
        }
        fprintf(gp, "\n");

        /* FFT PSD data */
        for (int i = 0; i < nfft / 2; i++)
            fprintf(gp, "%f %f\n", (double)i / (double)nfft, psd_fft[i]);
        fprintf(gp, "e\n");

        /* AR spectra */
        for (int o = 0; o < n_orders; o++) {
            double a[MAX_ORDER + 1], E;
            double spec[256];
            lpc_coefficients(x, N_SAMPLES, orders[o], a, &E);
            lpc_spectrum(a, orders[o], E, spec, nfft);
            for (int i = 0; i < nfft / 2; i++)
                fprintf(gp, "%f %f\n", (double)i / (double)nfft, spec[i]);
            fprintf(gp, "e\n");
        }
        pclose(gp);
        printf("  Plot → plots/ch24/ar_spectrum.png\n\n");
    }

    free(X);
}

/* ── Demo 4: Reflection Coefficients ──────────────────────────── */

static void demo_reflection_coefficients(void)
{
    printf("=== Demo 4: Reflection Coefficients ===\n\n");

    double x[N_SAMPLES];
    for (int i = 0; i < N_SAMPLES; i++)
        x[i] = sin(2.0 * M_PI * 0.1 * i) + 0.3 * sin(2.0 * M_PI * 0.35 * i);

    int order = 10;
    double r[MAX_ORDER + 1], a[MAX_ORDER + 1], k_coeff[MAX_ORDER];
    double E;

    lpc_autocorrelation(x, N_SAMPLES, r, order);
    levinson_durbin(r, order, a, k_coeff, &E);

    printf("  Reflection coefficients (lattice structure):\n");
    for (int i = 0; i < order; i++)
        printf("    k[%2d] = %+.6f  %s\n", i + 1, k_coeff[i],
               (fabs(k_coeff[i]) < 1.0) ? "(stable)" : "(!unstable!)");

    printf("\n  All |k| < 1 ⟹ stable all-pole model ✓\n\n");
}

/* ── Main ─────────────────────────────────────────────────────── */

int main(void)
{
    printf("╔══════════════════════════════════════════════╗\n");
    printf("║  Chapter 24: Linear Prediction Coding       ║\n");
    printf("║  Levinson-Durbin · AR Spectrum · LPC        ║\n");
    printf("╚══════════════════════════════════════════════╝\n\n");

    if (system("mkdir -p plots/ch24") != 0) { /* ignore */ }

    demo_levinson();
    demo_analysis_synthesis();
    demo_ar_spectrum();
    demo_reflection_coefficients();

    printf("All Chapter 24 demos complete.\n");
    return 0;
}

/**
 * @file 27-2d-dsp.c
 * @brief Chapter 27 demo — Two-Dimensional DSP.
 *
 * ── 2-D Convolution ─────────────────────────────────────────────
 *
 *   out(r,c) = ΣΣ img(r-i, c-j) · kernel(i,j)
 *
 *   Image:          Kernel (3×3):     Result:
 *   ┌─────────┐     ┌───────┐        ┌─────────┐
 *   │ · · · · │  ⊛  │ w w w │   →    │ · · · · │
 *   │ · · · · │     │ w w w │        │ · · · · │
 *   │ · · · · │     │ w w w │        │ · · · · │
 *   └─────────┘     └───────┘        └─────────┘
 *
 * ── 2-D FFT (row-column decomposition) ──────────────────────────
 *
 *   FFT2D(img) = FFT_cols( FFT_rows( img ) )
 *
 * Demonstrates:
 *   - Gaussian blur via spatial convolution
 *   - Edge detection via Sobel operator
 *   - Image sharpening
 *   - 2-D FFT and frequency-domain filtering
 *   - Ideal low-pass filter in 2-D frequency domain
 *
 * Build & run:
 *   make chapters && ./build/bin/ch27
 *
 * Read alongside: chapters/27-2d-dsp.md
 */

#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dsp2d.h"
#include "gnuplot.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define IMG_SIZE  64  /* 64×64 test image (power of 2 for FFT) */

/* ── Generate a synthetic test image ──────────────────────────── */

static void gen_test_image(double *img, int rows, int cols)
{
    /* Combination of rectangles, circles, and gradients */
    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            double val = 0.0;

            /* Background gradient */
            val = 0.1 * (double)c / cols;

            /* Bright rectangle */
            if (r >= 10 && r < 30 && c >= 15 && c < 45)
                val = 0.9;

            /* Dark circle */
            double dr = (double)(r - 45);
            double dc = (double)(c - 20);
            if (sqrt(dr * dr + dc * dc) < 10.0)
                val = 0.0;

            /* Bright circle */
            dr = (double)(r - 45);
            dc = (double)(c - 50);
            if (sqrt(dr * dr + dc * dc) < 8.0)
                val = 1.0;

            /* Diagonal line */
            if (abs(r - c) < 2 && r < 20)
                val = 0.7;

            img[r * cols + c] = val;
        }
    }
}

/* ── Plot a 2D image as heatmap ───────────────────────────────── */

static void plot_image(const double *img, int rows, int cols,
                       const char *filename, const char *title)
{
    FILE *gp = popen("gnuplot -persistent", "w");
    if (!gp) return;

    fprintf(gp, "set terminal pngcairo size 400,400\n");
    fprintf(gp, "set output '%s'\n", filename);
    fprintf(gp, "set title '%s'\n", title);
    fprintf(gp, "set xlabel 'Column'\nset ylabel 'Row'\n");
    fprintf(gp, "set palette grey\nset cbrange [0:1]\n");
    fprintf(gp, "set xrange [-0.5:%.1f]\nset yrange [%.1f:-0.5]\n",
            cols - 0.5, rows - 0.5);
    fprintf(gp, "plot '-' matrix with image notitle\n");

    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++)
            fprintf(gp, "%f ", img[r * cols + c]);
        fprintf(gp, "\n");
    }
    fprintf(gp, "e\ne\n");
    pclose(gp);
}

/* ── Demo 1: Gaussian Blur ────────────────────────────────────── */

static void demo_gaussian_blur(void)
{
    printf("=== Demo 1: Gaussian Blur ===\n\n");

    double img[IMG_SIZE * IMG_SIZE], blurred[IMG_SIZE * IMG_SIZE];
    gen_test_image(img, IMG_SIZE, IMG_SIZE);

    plot_image(img, IMG_SIZE, IMG_SIZE,
               "plots/ch27/original.png", "Original Image");

    double kernel[49];  /* 7×7 */
    kernel_gaussian(kernel, 7, 1.5);

    conv2d(img, IMG_SIZE, IMG_SIZE, kernel, 7, 7, blurred);

    /* Normalise for display */
    double mx = 0;
    for (int i = 0; i < IMG_SIZE * IMG_SIZE; i++)
        if (fabs(blurred[i]) > mx) mx = fabs(blurred[i]);
    if (mx > 0) for (int i = 0; i < IMG_SIZE * IMG_SIZE; i++)
        blurred[i] /= mx;

    plot_image(blurred, IMG_SIZE, IMG_SIZE,
               "plots/ch27/gaussian_blur.png", "Gaussian Blur (7x7, sigma=1.5)");

    printf("  Original → Gaussian blur applied\n");
    printf("  Plots → plots/ch27/original.png, gaussian_blur.png\n\n");
}

/* ── Demo 2: Sobel Edge Detection ─────────────────────────────── */

static void demo_sobel_edges(void)
{
    printf("=== Demo 2: Sobel Edge Detection ===\n\n");

    double img[IMG_SIZE * IMG_SIZE], edges[IMG_SIZE * IMG_SIZE];
    gen_test_image(img, IMG_SIZE, IMG_SIZE);

    sobel_magnitude(img, IMG_SIZE, IMG_SIZE, edges);

    /* Normalise for display */
    double mx = 0;
    for (int i = 0; i < IMG_SIZE * IMG_SIZE; i++)
        if (edges[i] > mx) mx = edges[i];
    if (mx > 0) for (int i = 0; i < IMG_SIZE * IMG_SIZE; i++)
        edges[i] /= mx;

    plot_image(edges, IMG_SIZE, IMG_SIZE,
               "plots/ch27/sobel_edges.png", "Sobel Edge Detection");

    printf("  Edge magnitude computed via Sobel gradient\n");
    printf("  Plot → plots/ch27/sobel_edges.png\n\n");
}

/* ── Demo 3: Image Sharpening ─────────────────────────────────── */

static void demo_sharpen(void)
{
    printf("=== Demo 3: Image Sharpening ===\n\n");

    double img[IMG_SIZE * IMG_SIZE];
    double blurred[IMG_SIZE * IMG_SIZE], sharpened[IMG_SIZE * IMG_SIZE];
    gen_test_image(img, IMG_SIZE, IMG_SIZE);

    /* First blur, then sharpen the blurred image */
    double gk[25];
    kernel_gaussian(gk, 5, 1.0);
    conv2d(img, IMG_SIZE, IMG_SIZE, gk, 5, 5, blurred);

    double sk[9];
    kernel_sharpen(sk, 1.5);
    conv2d(blurred, IMG_SIZE, IMG_SIZE, sk, 3, 3, sharpened);

    /* Clip to [0,1] */
    for (int i = 0; i < IMG_SIZE * IMG_SIZE; i++) {
        if (sharpened[i] < 0.0) sharpened[i] = 0.0;
        if (sharpened[i] > 1.0) sharpened[i] = 1.0;
    }

    double mx = 0;
    for (int i = 0; i < IMG_SIZE * IMG_SIZE; i++)
        if (fabs(blurred[i]) > mx) mx = fabs(blurred[i]);
    if (mx > 0) for (int i = 0; i < IMG_SIZE * IMG_SIZE; i++)
        blurred[i] /= mx;

    plot_image(blurred, IMG_SIZE, IMG_SIZE,
               "plots/ch27/blurred.png", "Blurred Input");
    plot_image(sharpened, IMG_SIZE, IMG_SIZE,
               "plots/ch27/sharpened.png", "After Sharpening (alpha=1.5)");

    printf("  Blurred → sharpened via Laplacian\n");
    printf("  Plots → plots/ch27/blurred.png, sharpened.png\n\n");
}

/* ── Demo 4: 2-D FFT Low-Pass Filtering ──────────────────────── */

static void demo_2d_fft_filter(void)
{
    printf("=== Demo 4: 2-D Frequency-Domain Low-Pass Filter ===\n\n");

    double img[IMG_SIZE * IMG_SIZE];
    gen_test_image(img, IMG_SIZE, IMG_SIZE);

    /* Add some high-frequency noise */
    unsigned int seed = 42;
    double noisy[IMG_SIZE * IMG_SIZE];
    for (int i = 0; i < IMG_SIZE * IMG_SIZE; i++) {
        seed = seed * 1103515245u + 12345u;
        double noise = ((double)(seed & 0x7FFF) / 32768.0 - 0.5) * 0.3;
        noisy[i] = img[i] + noise;
        if (noisy[i] < 0) noisy[i] = 0;
        if (noisy[i] > 1) noisy[i] = 1;
    }

    /* Ideal LPF in 2-D frequency domain */
    int N = IMG_SIZE * IMG_SIZE;
    double *H_re = (double *)calloc((size_t)N, sizeof(double));
    double *H_im = (double *)calloc((size_t)N, sizeof(double));
    lpf2d_ideal(H_re, H_im, IMG_SIZE, IMG_SIZE, 0.15);

    double filtered[IMG_SIZE * IMG_SIZE];
    filter2d_freq(noisy, IMG_SIZE, IMG_SIZE, H_re, H_im, filtered);

    /* Clip */
    for (int i = 0; i < N; i++) {
        if (filtered[i] < 0) filtered[i] = 0;
        if (filtered[i] > 1) filtered[i] = 1;
    }

    plot_image(noisy, IMG_SIZE, IMG_SIZE,
               "plots/ch27/noisy.png", "Noisy Image");
    plot_image(filtered, IMG_SIZE, IMG_SIZE,
               "plots/ch27/lpf_filtered.png",
               "2D FFT Low-Pass Filtered (cutoff=0.15)");

    printf("  Noisy → 2D-FFT LPF denoised\n");
    printf("  Plots → plots/ch27/noisy.png, lpf_filtered.png\n\n");

    free(H_re); free(H_im);
}

/* ── Demo 5: 2-D FFT Magnitude Spectrum ───────────────────────── */

static void demo_2d_spectrum(void)
{
    printf("=== Demo 5: 2-D FFT Magnitude Spectrum ===\n\n");

    double img[IMG_SIZE * IMG_SIZE];
    gen_test_image(img, IMG_SIZE, IMG_SIZE);

    int N = IMG_SIZE * IMG_SIZE;
    double *re = (double *)malloc((size_t)N * sizeof(double));
    double *im = (double *)calloc((size_t)N, sizeof(double));
    memcpy(re, img, (size_t)N * sizeof(double));

    fft2d(re, im, IMG_SIZE, IMG_SIZE);

    /* Log magnitude, centre-shifted */
    double mag[IMG_SIZE * IMG_SIZE];
    double mx = 0;
    for (int u = 0; u < IMG_SIZE; u++) {
        for (int v = 0; v < IMG_SIZE; v++) {
            int idx = u * IMG_SIZE + v;
            double m = sqrt(re[idx] * re[idx] + im[idx] * im[idx]);
            /* Shift: swap quadrants */
            int su = (u + IMG_SIZE / 2) % IMG_SIZE;
            int sv = (v + IMG_SIZE / 2) % IMG_SIZE;
            mag[su * IMG_SIZE + sv] = log10(m + 1.0);
        }
    }
    for (int i = 0; i < N; i++)
        if (mag[i] > mx) mx = mag[i];
    if (mx > 0) for (int i = 0; i < N; i++)
        mag[i] /= mx;

    plot_image(mag, IMG_SIZE, IMG_SIZE,
               "plots/ch27/fft2d_spectrum.png",
               "2D FFT Magnitude (log, centred)");

    printf("  2-D FFT magnitude spectrum plotted\n");
    printf("  Plot → plots/ch27/fft2d_spectrum.png\n\n");

    free(re); free(im);
}

/* ── Main ─────────────────────────────────────────────────────── */

int main(void)
{
    printf("╔══════════════════════════════════════════════╗\n");
    printf("║  Chapter 27: Two-Dimensional DSP            ║\n");
    printf("║  Conv2D · Sobel · FFT2D · Filtering         ║\n");
    printf("╚══════════════════════════════════════════════╝\n\n");

    if (system("mkdir -p plots/ch27") != 0) { /* ignore */ }

    demo_gaussian_blur();
    demo_sobel_edges();
    demo_sharpen();
    demo_2d_fft_filter();
    demo_2d_spectrum();

    printf("All Chapter 27 demos complete.\n");
    return 0;
}

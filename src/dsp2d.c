/**
 * @file dsp2d.c
 * @brief 2-D DSP — convolution, FFT (row-column), image kernels.
 *
 * ── 2-D FFT Strategy ────────────────────────────────────────────
 *   1. Apply 1-D FFT to each row
 *   2. Apply 1-D FFT to each column
 *   Complexity: O(MN·log(MN))   for M×N image
 *
 * ── Spatial Convolution ──────────────────────────────────────────
 *   out(r,c) = ΣΣ img(r-i, c-j) · kernel(i,j)
 *   Zero-padded boundary conditions.
 */

#include "dsp2d.h"
#include "dsp_utils.h"  /* Complex */
#include "fft.h"        /* fft, ifft */
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ================================================================== */
/*  2-D Spatial Convolution                                            */
/* ================================================================== */

void conv2d(const double *img, int rows, int cols,
            const double *kernel, int krows, int kcols,
            double *out)
{
    int kr2 = krows / 2;
    int kc2 = kcols / 2;

    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            double sum = 0.0;
            for (int ki = 0; ki < krows; ki++) {
                for (int kj = 0; kj < kcols; kj++) {
                    int ri = r + ki - kr2;
                    int ci = c + kj - kc2;
                    if (ri >= 0 && ri < rows && ci >= 0 && ci < cols)
                        sum += img[ri * cols + ci] * kernel[ki * kcols + kj];
                }
            }
            out[r * cols + c] = sum;
        }
    }
}

/* ================================================================== */
/*  Standard Kernels                                                   */
/* ================================================================== */

void kernel_gaussian(double *kernel, int ksize, double sigma)
{
    int half = ksize / 2;
    double sum = 0.0;
    double s2 = 2.0 * sigma * sigma;

    for (int i = 0; i < ksize; i++) {
        for (int j = 0; j < ksize; j++) {
            double di = (double)(i - half);
            double dj = (double)(j - half);
            double val = exp(-(di * di + dj * dj) / s2);
            kernel[i * ksize + j] = val;
            sum += val;
        }
    }
    /* Normalise */
    for (int i = 0; i < ksize * ksize; i++)
        kernel[i] /= sum;
}

void kernel_sobel(double *gx, double *gy)
{
    /* Sobel X (horizontal gradient) */
    static const double sx[9] = {
        -1, 0, 1,
        -2, 0, 2,
        -1, 0, 1
    };
    /* Sobel Y (vertical gradient) */
    static const double sy[9] = {
        -1, -2, -1,
         0,  0,  0,
         1,  2,  1
    };
    memcpy(gx, sx, 9 * sizeof(double));
    memcpy(gy, sy, 9 * sizeof(double));
}

void kernel_log(double *kernel, int ksize, double sigma)
{
    int half = ksize / 2;
    double s2 = sigma * sigma;
    double s4 = s2 * s2;
    double sum = 0.0;

    for (int i = 0; i < ksize; i++) {
        for (int j = 0; j < ksize; j++) {
            double di = (double)(i - half);
            double dj = (double)(j - half);
            double r2 = di * di + dj * dj;
            double val = -(1.0 / (M_PI * s4)) *
                         (1.0 - r2 / (2.0 * s2)) *
                         exp(-r2 / (2.0 * s2));
            kernel[i * ksize + j] = val;
            sum += val;
        }
    }
    /* Zero-mean normalisation */
    double mean = sum / (double)(ksize * ksize);
    for (int i = 0; i < ksize * ksize; i++)
        kernel[i] -= mean;
}

void kernel_sharpen(double *kernel, double alpha)
{
    /* Centre = 1 + 4α, edges = -α, corners = 0 */
    double k[9] = {
        0.0,    -alpha,       0.0,
        -alpha, 1.0 + 4*alpha, -alpha,
        0.0,    -alpha,       0.0
    };
    memcpy(kernel, k, 9 * sizeof(double));
}

/* ================================================================== */
/*  Sobel Gradient Magnitude                                           */
/* ================================================================== */

void sobel_magnitude(const double *img, int rows, int cols, double *mag)
{
    double gx[9], gy[9];
    kernel_sobel(gx, gy);

    double *Gx = (double *)malloc((size_t)(rows * cols) * sizeof(double));
    double *Gy = (double *)malloc((size_t)(rows * cols) * sizeof(double));

    conv2d(img, rows, cols, gx, 3, 3, Gx);
    conv2d(img, rows, cols, gy, 3, 3, Gy);

    for (int i = 0; i < rows * cols; i++)
        mag[i] = sqrt(Gx[i] * Gx[i] + Gy[i] * Gy[i]);

    free(Gx); free(Gy);
}

/* ================================================================== */
/*  2-D FFT (row-column decomposition)                                 */
/* ================================================================== */

void fft2d(double *data_re, double *data_im, int rows, int cols)
{
    Complex *tmp = (Complex *)malloc((size_t)(rows > cols ? rows : cols) * sizeof(Complex));

    /* FFT along each row */
    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            tmp[c].re = data_re[r * cols + c];
            tmp[c].im = data_im[r * cols + c];
        }
        fft(tmp, cols);
        for (int c = 0; c < cols; c++) {
            data_re[r * cols + c] = tmp[c].re;
            data_im[r * cols + c] = tmp[c].im;
        }
    }

    /* FFT along each column */
    for (int c = 0; c < cols; c++) {
        for (int r = 0; r < rows; r++) {
            tmp[r].re = data_re[r * cols + c];
            tmp[r].im = data_im[r * cols + c];
        }
        fft(tmp, rows);
        for (int r = 0; r < rows; r++) {
            data_re[r * cols + c] = tmp[r].re;
            data_im[r * cols + c] = tmp[r].im;
        }
    }

    free(tmp);
}

void ifft2d(double *data_re, double *data_im, int rows, int cols)
{
    Complex *tmp = (Complex *)malloc((size_t)(rows > cols ? rows : cols) * sizeof(Complex));

    /* IFFT along each row */
    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            tmp[c].re = data_re[r * cols + c];
            tmp[c].im = data_im[r * cols + c];
        }
        ifft(tmp, cols);
        for (int c = 0; c < cols; c++) {
            data_re[r * cols + c] = tmp[c].re;
            data_im[r * cols + c] = tmp[c].im;
        }
    }

    /* IFFT along each column */
    for (int c = 0; c < cols; c++) {
        for (int r = 0; r < rows; r++) {
            tmp[r].re = data_re[r * cols + c];
            tmp[r].im = data_im[r * cols + c];
        }
        ifft(tmp, rows);
        for (int r = 0; r < rows; r++) {
            data_re[r * cols + c] = tmp[r].re;
            data_im[r * cols + c] = tmp[r].im;
        }
    }

    free(tmp);
}

/* ================================================================== */
/*  Frequency-Domain 2-D Filtering                                     */
/* ================================================================== */

void filter2d_freq(const double *img, int rows, int cols,
                   const double *H_re, const double *H_im,
                   double *out)
{
    int N = rows * cols;
    double *re = (double *)malloc((size_t)N * sizeof(double));
    double *im = (double *)calloc((size_t)N, sizeof(double));

    memcpy(re, img, (size_t)N * sizeof(double));

    fft2d(re, im, rows, cols);

    /* Multiply: X·H */
    for (int i = 0; i < N; i++) {
        double a = re[i], b = im[i];
        double c = H_re[i], d = H_im[i];
        re[i] = a * c - b * d;
        im[i] = a * d + b * c;
    }

    ifft2d(re, im, rows, cols);

    for (int i = 0; i < N; i++)
        out[i] = re[i];

    free(re); free(im);
}

void lpf2d_ideal(double *H_re, double *H_im, int rows, int cols,
                 double cutoff)
{
    for (int u = 0; u < rows; u++) {
        for (int v = 0; v < cols; v++) {
            /* Centre-shifted frequency */
            double fu = (double)u / (double)rows;
            double fv = (double)v / (double)cols;
            if (fu > 0.5) fu -= 1.0;
            if (fv > 0.5) fv -= 1.0;
            double dist = sqrt(fu * fu + fv * fv);

            H_re[u * cols + v] = (dist <= cutoff) ? 1.0 : 0.0;
            H_im[u * cols + v] = 0.0;
        }
    }
}

/**
 * @file dsp2d.h
 * @brief 2-D DSP — convolution, FFT, image kernels.
 *
 * Chapter 27: Two-Dimensional DSP
 *
 * All images represented as row-major double arrays: pixel(r,c) = img[r*cols+c].
 * Greyscale only (single channel).
 */

#ifndef DSP2D_H
#define DSP2D_H

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------ */
/*  2-D Convolution                                                    */
/* ------------------------------------------------------------------ */

/**
 * 2-D linear convolution (direct, zero-padded boundary).
 *
 * @param img       Input image (rows × cols)
 * @param rows      Image height
 * @param cols      Image width
 * @param kernel    Convolution kernel (krows × kcols)
 * @param krows     Kernel height
 * @param kcols     Kernel width
 * @param out       Output image (rows × cols)
 */
void conv2d(const double *img, int rows, int cols,
            const double *kernel, int krows, int kcols,
            double *out);

/* ------------------------------------------------------------------ */
/*  Standard Kernels                                                   */
/* ------------------------------------------------------------------ */

/**
 * Generate a Gaussian blur kernel (ksize × ksize).
 * @param kernel    Output (ksize × ksize)
 * @param ksize     Kernel size (odd, e.g. 3, 5, 7)
 * @param sigma     Standard deviation
 */
void kernel_gaussian(double *kernel, int ksize, double sigma);

/**
 * Sobel gradient kernels (3×3).
 * @param gx   Output Sobel-X kernel (3×3)
 * @param gy   Output Sobel-Y kernel (3×3)
 */
void kernel_sobel(double *gx, double *gy);

/**
 * Laplacian of Gaussian (LoG) kernel.
 * @param kernel    Output (ksize × ksize)
 * @param ksize     Kernel size (odd)
 * @param sigma     Standard deviation
 */
void kernel_log(double *kernel, int ksize, double sigma);

/**
 * Sharpening kernel (3×3): identity + α·Laplacian.
 * @param kernel    Output (3×3)
 * @param alpha     Sharpening strength
 */
void kernel_sharpen(double *kernel, double alpha);

/* ------------------------------------------------------------------ */
/*  Edge Detection via Sobel Gradient                                  */
/* ------------------------------------------------------------------ */

/**
 * Compute gradient magnitude via Sobel operator.
 * @param img   Input image (rows × cols)
 * @param rows  Height
 * @param cols  Width
 * @param mag   Output gradient magnitude (rows × cols)
 */
void sobel_magnitude(const double *img, int rows, int cols, double *mag);

/* ------------------------------------------------------------------ */
/*  2-D FFT (row-column decomposition)                                 */
/* ------------------------------------------------------------------ */

/**
 * Forward 2-D FFT (in-place, rows × cols must both be powers of 2).
 * Data layout: row-major Complex array [rows * cols].
 *
 * @param data  Complex image (rows × cols), in-place
 * @param rows  Height (power of 2)
 * @param cols  Width (power of 2)
 */
void fft2d(double *data_re, double *data_im, int rows, int cols);

/**
 * Inverse 2-D FFT (in-place).
 * @param data_re  Real part (rows × cols)
 * @param data_im  Imag part (rows × cols)
 * @param rows     Height (power of 2)
 * @param cols     Width (power of 2)
 */
void ifft2d(double *data_re, double *data_im, int rows, int cols);

/**
 * 2-D frequency-domain filtering.
 * Steps: FFT2D → multiply by H → IFFT2D.
 *
 * @param img       Input image (rows × cols, real)
 * @param rows      Height (power of 2)
 * @param cols      Width (power of 2)
 * @param H_re      Filter freq response real (rows × cols)
 * @param H_im      Filter freq response imag (rows × cols)
 * @param out       Output image (rows × cols, real)
 */
void filter2d_freq(const double *img, int rows, int cols,
                   const double *H_re, const double *H_im,
                   double *out);

/**
 * Generate ideal low-pass 2-D filter in frequency domain.
 * H(u,v) = 1 if sqrt(u² + v²) <= cutoff else 0.
 *
 * @param H_re      Output real part (rows × cols)
 * @param H_im      Output imag part (rows × cols, zeroed)
 * @param rows      Height
 * @param cols      Width
 * @param cutoff    Cutoff (normalised, 0..0.5)
 */
void lpf2d_ideal(double *H_re, double *H_im, int rows, int cols,
                 double cutoff);

#ifdef __cplusplus
}
#endif

#endif /* DSP2D_H */

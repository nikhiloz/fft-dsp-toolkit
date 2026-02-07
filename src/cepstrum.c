/**
 * @file cepstrum.c
 * @brief Cepstral analysis — real/complex cepstrum, liftering, MFCC pipeline.
 *
 * ── Real Cepstrum ────────────────────────────────────────────────
 *   c[n] = IFFT{ log|FFT{x}| }
 *   Symmetric: c[n] = c[-n].  Separates convolution in quefrency domain.
 *
 * ── MFCC Pipeline ────────────────────────────────────────────────
 *   frame → Hamming → FFT → |X|² → Mel filterbank → log → DCT-II
 */

#include "cepstrum.h"
#include "dsp_utils.h"   /* Complex */
#include "fft.h"         /* fft, ifft */
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ================================================================== */
/*  Real Cepstrum                                                      */
/* ================================================================== */

void cepstrum_real(const double *x, int n, double *c, int nfft)
{
    Complex *X = (Complex *)calloc((size_t)nfft, sizeof(Complex));
    for (int i = 0; i < n && i < nfft; i++) {
        X[i].re = x[i];
        X[i].im = 0.0;
    }

    fft(X, nfft);

    /* log|X[k]| */
    for (int k = 0; k < nfft; k++) {
        double mag = sqrt(X[k].re * X[k].re + X[k].im * X[k].im);
        X[k].re = log(mag + 1e-30);
        X[k].im = 0.0;
    }

    ifft(X, nfft);

    for (int i = 0; i < nfft; i++)
        c[i] = X[i].re;

    free(X);
}

/* ================================================================== */
/*  Complex Cepstrum                                                   */
/* ================================================================== */

void cepstrum_complex(const double *x, int n, double *c, int nfft)
{
    Complex *X = (Complex *)calloc((size_t)nfft, sizeof(Complex));
    for (int i = 0; i < n && i < nfft; i++) {
        X[i].re = x[i];
        X[i].im = 0.0;
    }

    fft(X, nfft);

    /* log(X[k]) = log|X[k]| + j·phase_unwrapped(X[k]) */
    /* Simple phase unwrapping */
    double prev_phase = 0.0;
    double unwrapped = 0.0;
    for (int k = 0; k < nfft; k++) {
        double mag = sqrt(X[k].re * X[k].re + X[k].im * X[k].im);
        double phase = atan2(X[k].im, X[k].re);

        if (k > 0) {
            double dp = phase - prev_phase;
            /* Wrap to [-π, π] */
            while (dp > M_PI) dp -= 2.0 * M_PI;
            while (dp < -M_PI) dp += 2.0 * M_PI;
            unwrapped += dp;
        } else {
            unwrapped = phase;
        }
        prev_phase = phase;

        X[k].re = log(mag + 1e-30);
        X[k].im = unwrapped;
    }

    ifft(X, nfft);

    for (int i = 0; i < nfft; i++)
        c[i] = X[i].re;

    free(X);
}

/* ================================================================== */
/*  Liftering (spectral envelope extraction)                           */
/* ================================================================== */

void cepstrum_lifter(const double *c, int nfft, int L, double *envelope)
{
    Complex *X = (Complex *)calloc((size_t)nfft, sizeof(Complex));

    /* Low-time lifter: keep c[0..L-1] and mirror c[nfft-L+1..nfft-1] */
    for (int i = 0; i < L && i < nfft; i++) {
        X[i].re = c[i];
        X[i].im = 0.0;
    }
    /* Mirror for real-valued result */
    for (int i = 1; i < L && (nfft - i) >= L; i++) {
        X[nfft - i].re = c[nfft - i];
        X[nfft - i].im = 0.0;
    }

    fft(X, nfft);

    /* Output: log magnitude envelope in dB (positive frequencies) */
    int half = nfft / 2;
    for (int k = 0; k < half; k++) {
        double mag = sqrt(X[k].re * X[k].re + X[k].im * X[k].im);
        envelope[k] = 20.0 * log10(mag + 1e-30);
    }

    free(X);
}

/* ================================================================== */
/*  Mel Scale Conversions                                              */
/* ================================================================== */

double hz_to_mel(double f_hz)
{
    return 2595.0 * log10(1.0 + f_hz / 700.0);
}

double mel_to_hz(double mel)
{
    return 700.0 * (pow(10.0, mel / 2595.0) - 1.0);
}

/* ================================================================== */
/*  Mel Filterbank                                                     */
/* ================================================================== */

void mel_filterbank(const double *power_spec, int nfft, double fs,
                    int n_filters, double *fbank,
                    double f_low, double f_high)
{
    int half = nfft / 2;

    /* Compute n_filters + 2 equally-spaced points on Mel scale */
    double mel_low = hz_to_mel(f_low);
    double mel_high = hz_to_mel(f_high);
    int n_pts = n_filters + 2;
    double *mel_pts = (double *)malloc((size_t)n_pts * sizeof(double));
    double *hz_pts = (double *)malloc((size_t)n_pts * sizeof(double));
    int *bin_pts = (int *)malloc((size_t)n_pts * sizeof(int));

    for (int i = 0; i < n_pts; i++) {
        mel_pts[i] = mel_low + (mel_high - mel_low) * (double)i / (double)(n_pts - 1);
        hz_pts[i] = mel_to_hz(mel_pts[i]);
        bin_pts[i] = (int)floor(hz_pts[i] * (double)nfft / fs);
        if (bin_pts[i] >= half) bin_pts[i] = half - 1;
        if (bin_pts[i] < 0) bin_pts[i] = 0;
    }

    /* Apply triangular filters */
    for (int m = 0; m < n_filters; m++) {
        int f_start = bin_pts[m];
        int f_center = bin_pts[m + 1];
        int f_end = bin_pts[m + 2];

        double sum = 0.0;
        for (int k = f_start; k <= f_end && k < half; k++) {
            double weight = 0.0;
            if (k >= f_start && k <= f_center && f_center > f_start)
                weight = (double)(k - f_start) / (double)(f_center - f_start);
            else if (k > f_center && k <= f_end && f_end > f_center)
                weight = (double)(f_end - k) / (double)(f_end - f_center);
            sum += power_spec[k] * weight;
        }
        fbank[m] = sum;
    }

    free(mel_pts); free(hz_pts); free(bin_pts);
}

/* ================================================================== */
/*  DCT Type-II                                                        */
/* ================================================================== */

void dct_ii(const double *x, double *y, int n)
{
    for (int k = 0; k < n; k++) {
        double sum = 0.0;
        for (int i = 0; i < n; i++)
            sum += x[i] * cos(M_PI * (2.0 * (double)i + 1.0) * (double)k
                              / (2.0 * (double)n));
        y[k] = sum;
    }
}

/* ================================================================== */
/*  MFCC Computation                                                   */
/* ================================================================== */

void compute_mfcc(const double *frame, int frame_len, int nfft,
                  double fs, int n_filters, int n_mfcc, double *mfcc)
{
    /* 1. Apply Hamming window */
    double *windowed = (double *)calloc((size_t)nfft, sizeof(double));
    for (int i = 0; i < frame_len && i < nfft; i++) {
        double w = 0.54 - 0.46 * cos(2.0 * M_PI * (double)i / (double)(frame_len - 1));
        windowed[i] = frame[i] * w;
    }

    /* 2. FFT → power spectrum */
    Complex *X = (Complex *)calloc((size_t)nfft, sizeof(Complex));
    for (int i = 0; i < nfft; i++) {
        X[i].re = windowed[i];
        X[i].im = 0.0;
    }
    fft(X, nfft);

    int half = nfft / 2;
    double *power_spec = (double *)malloc((size_t)half * sizeof(double));
    for (int k = 0; k < half; k++)
        power_spec[k] = X[k].re * X[k].re + X[k].im * X[k].im;

    /* 3. Mel filterbank */
    double *fbank = (double *)malloc((size_t)n_filters * sizeof(double));
    mel_filterbank(power_spec, nfft, fs, n_filters, fbank, 0.0, fs / 2.0);

    /* 4. Log filterbank energies */
    double *log_fbank = (double *)malloc((size_t)n_filters * sizeof(double));
    for (int m = 0; m < n_filters; m++)
        log_fbank[m] = log(fbank[m] + 1e-30);

    /* 5. DCT-II → MFCCs */
    double *dct_out = (double *)malloc((size_t)n_filters * sizeof(double));
    dct_ii(log_fbank, dct_out, n_filters);

    for (int i = 0; i < n_mfcc && i < n_filters; i++)
        mfcc[i] = dct_out[i];

    free(windowed); free(X); free(power_spec);
    free(fbank); free(log_fbank); free(dct_out);
}

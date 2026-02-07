/**
 * @file test_phase6.c
 * @brief Unit tests for Phase 6 modules: adaptive, lpc, spectral_est,
 *        cepstrum, dsp2d.
 *
 * Tests:
 *   1.  LMS converges to known plant weights
 *   2.  NLMS steady-state MSE < threshold
 *   3.  RLS converges faster than LMS
 *   4.  Levinson-Durbin recovers AR coefficients
 *   5.  LPC analysis/synthesis round-trip
 *   6.  LPC spectral envelope has peaks at signal frequencies
 *   7.  Reflection coefficients all |k| < 1 for stable signal
 *   8.  MUSIC finds two sinusoidal frequencies
 *   9.  Capon spectrum has peaks at signal frequencies
 *  10.  Correlation matrix is symmetric
 *  11.  Real cepstrum of impulse is near-zero
 *  12.  Cepstral pitch detection matches known period
 *  13.  Mel scale round-trip: Hz → Mel → Hz
 *  14.  MFCC output has valid values (no NaN/Inf)
 *  15.  DCT-II of constant is non-zero only at index 0
 *  16.  Conv2D with identity kernel preserves image
 *  17.  Gaussian kernel sums to 1.0
 *  18.  Sobel of flat image returns zeros
 *  19.  2D FFT → IFFT round-trip preserves data
 *
 * Run: make test
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "test_framework.h"
#include "adaptive.h"
#include "lpc.h"
#include "spectral_est.h"
#include "cepstrum.h"
#include "dsp2d.h"
#include "filter.h"
#include "fft.h"
#include "dsp_utils.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Simple pseudo-random from LCG */
static double test_randn(unsigned int *seed)
{
    *seed = *seed * 1103515245u + 12345u;
    double u1 = ((double)(*seed & 0x7FFFFFFFu)) / 2147483647.0;
    *seed = *seed * 1103515245u + 12345u;
    double u2 = ((double)(*seed & 0x7FFFFFFFu)) / 2147483647.0;
    if (u1 < 1e-10) u1 = 1e-10;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

int main(void)
{
    TEST_SUITE("Phase 6: Adaptive, LPC, Spectral Est, Cepstrum, 2D DSP");

    /* ── Test 1: LMS converges to known plant ────────────── */
    TEST_CASE_BEGIN("LMS converges to known plant weights");
    {
        int taps = 7;
        double plant[7];
        fir_lowpass(plant, taps, 0.2);

        int N = 4000;
        double *x = (double *)malloc((size_t)N * sizeof(double));
        double *d = (double *)malloc((size_t)N * sizeof(double));
        double *y = (double *)malloc((size_t)N * sizeof(double));
        double *e = (double *)malloc((size_t)N * sizeof(double));
        double w[7];
        unsigned int seed = 42;
        for (int i = 0; i < N; i++) x[i] = test_randn(&seed);
        fir_filter(x, d, N, plant, taps);

        lms_filter(x, d, N, taps, 0.01, y, e, w);

        double max_diff = 0;
        for (int i = 0; i < taps; i++) {
            double diff = fabs(w[i] - plant[i]);
            if (diff > max_diff) max_diff = diff;
        }
        if (max_diff < 0.05) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("LMS weight error > 0.05"); }

        free(x); free(d); free(y); free(e);
    }

    /* ── Test 2: NLMS steady-state MSE ────────────────────── */
    TEST_CASE_BEGIN("NLMS steady-state MSE < 1e-3");
    {
        int taps = 7;
        double plant[7];
        fir_lowpass(plant, taps, 0.2);

        int N = 4000;
        double *x = (double *)malloc((size_t)N * sizeof(double));
        double *d = (double *)malloc((size_t)N * sizeof(double));
        double *y = (double *)malloc((size_t)N * sizeof(double));
        double *e = (double *)malloc((size_t)N * sizeof(double));
        double w[7];
        unsigned int seed = 55;
        for (int i = 0; i < N; i++) x[i] = test_randn(&seed);
        fir_filter(x, d, N, plant, taps);

        nlms_filter(x, d, N, taps, 0.5, 1e-6, y, e, w);

        double mse = 0;
        for (int i = N / 2; i < N; i++) mse += e[i] * e[i];
        mse /= (double)(N / 2);
        if (mse < 1e-3) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("NLMS MSE too high"); }

        free(x); free(d); free(y); free(e);
    }

    /* ── Test 3: RLS converges faster than LMS ────────────── */
    TEST_CASE_BEGIN("RLS converges faster than LMS");
    {
        int taps = 5, N = 1000;
        double plant[5];
        fir_lowpass(plant, taps, 0.25);
        double *x = (double *)malloc((size_t)N * sizeof(double));
        double *d = (double *)malloc((size_t)N * sizeof(double));
        double *y = (double *)malloc((size_t)N * sizeof(double));
        double *e_lms = (double *)malloc((size_t)N * sizeof(double));
        double *e_rls = (double *)malloc((size_t)N * sizeof(double));
        double w[5];
        unsigned int seed = 99;
        for (int i = 0; i < N; i++) x[i] = test_randn(&seed);
        fir_filter(x, d, N, plant, taps);

        lms_filter(x, d, N, taps, 0.01, y, e_lms, w);
        rls_filter(x, d, N, taps, 0.99, 100.0, y, e_rls, w);

        /* Compare MSE in first 200 samples */
        double mse_lms = 0, mse_rls = 0;
        int check = 200;
        for (int i = 50; i < check; i++) {
            mse_lms += e_lms[i] * e_lms[i];
            mse_rls += e_rls[i] * e_rls[i];
        }
        if (mse_rls < mse_lms) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("RLS not faster than LMS"); }

        free(x); free(d); free(y); free(e_lms); free(e_rls);
    }

    /* ── Test 4: Levinson-Durbin recovers AR coefficients ── */
    TEST_CASE_BEGIN("Levinson-Durbin recovers AR(2) coefficients");
    {
        /* AR(2): x[n] = 1.5x[n-1] - 0.7x[n-2] + e[n] */
        int N = 4000;
        double *x = (double *)malloc((size_t)N * sizeof(double));
        unsigned int seed = 42;
        x[0] = 0; x[1] = 0;
        for (int i = 2; i < N; i++) {
            x[i] = 1.5 * x[i - 1] - 0.7 * x[i - 2] + test_randn(&seed) * 0.1;
        }

        double a[2], E;
        lpc_coefficients(x, N, 2, a, &E);

        /* a[0] should be ≈ -1.5, a[1] should be ≈ 0.7 (0-based) */
        int ok = (fabs(a[0] - (-1.5)) < 0.1) && (fabs(a[1] - 0.7) < 0.1);
        if (ok) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("LPC coefficients mismatch"); }

        free(x);
    }

    /* ── Test 5: LPC analysis/synthesis round-trip ─────────── */
    TEST_CASE_BEGIN("LPC analysis/synthesis round-trip");
    {
        int N = 256, order = 10;
        double *x = (double *)malloc((size_t)N * sizeof(double));
        for (int i = 0; i < N; i++)
            x[i] = sin(2.0 * M_PI * 0.1 * i) + 0.5 * sin(2.0 * M_PI * 0.3 * i);

        double a[11], E;
        lpc_coefficients(x, N, order, a, &E);

        double *residual = (double *)malloc((size_t)N * sizeof(double));
        double *recon = (double *)malloc((size_t)N * sizeof(double));
        lpc_residual(x, N, a, order, residual);
        lpc_synthesise(residual, N, a, order, recon);

        double max_err = 0;
        for (int i = 0; i < N; i++) {
            double diff = fabs(x[i] - recon[i]);
            if (diff > max_err) max_err = diff;
        }
        if (max_err < 1e-10) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("Round-trip error too large"); }

        free(x); free(residual); free(recon);
    }

    /* ── Test 6: LPC spectral envelope peaks ──────────────── */
    TEST_CASE_BEGIN("LPC spectral envelope has peaks at signal freq");
    {
        int N = 512, order = 12, nfft = 512;
        double *x = (double *)malloc((size_t)N * sizeof(double));
        for (int i = 0; i < N; i++)
            x[i] = sin(2.0 * M_PI * 0.1 * i);

        double a[13], E;
        lpc_coefficients(x, N, order, a, &E);
        double *spec = (double *)malloc((size_t)(nfft / 2) * sizeof(double));
        lpc_spectrum(a, order, E, spec, nfft);

        /* Find max in spec — should be near bin nfft*0.1 = 51 */
        int peak_bin = 0;
        double peak_val = -1e30;
        for (int i = 1; i < nfft / 2; i++)
            if (spec[i] > peak_val) { peak_val = spec[i]; peak_bin = i; }

        int expected = (int)(0.1 * nfft);
        if (abs(peak_bin - expected) < 5) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("AR peak not near expected freq"); }

        free(x); free(spec);
    }

    /* ── Test 7: Reflection coefficients all |k| < 1 ──────── */
    TEST_CASE_BEGIN("Reflection coefficients |k| < 1 (stable)");
    {
        int N = 512, order = 8;
        double *x = (double *)malloc((size_t)N * sizeof(double));
        for (int i = 0; i < N; i++)
            x[i] = sin(2.0 * M_PI * 0.15 * i) + 0.3 * sin(2.0 * M_PI * 0.35 * i);

        double r[9], a[9], k_coeff[8], E;
        lpc_autocorrelation(x, N, r, order);
        levinson_durbin(r, order, a, k_coeff, &E);

        int all_stable = 1;
        for (int i = 0; i < order; i++)
            if (fabs(k_coeff[i]) >= 1.0) all_stable = 0;

        if (all_stable) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("Some |k| >= 1"); }

        free(x);
    }

    /* ── Test 8: MUSIC finds two sinusoidal frequencies ───── */
    TEST_CASE_BEGIN("MUSIC finds two sinusoidal frequencies");
    {
        double f1 = 0.10, f2 = 0.25;
        int N = 256;
        double *x = (double *)malloc((size_t)N * sizeof(double));
        unsigned int seed = 42;
        for (int i = 0; i < N; i++) {
            x[i] = sin(2.0 * M_PI * f1 * i) + sin(2.0 * M_PI * f2 * i)
                  + test_randn(&seed) * 0.05;
        }

        double freqs[2];
        music_frequencies(x, N, 20, 2, freqs);

        int ok = (fabs(freqs[0] - f1) < 0.02) && (fabs(freqs[1] - f2) < 0.02);
        if (ok) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("MUSIC freq estimates inaccurate"); }

        free(x);
    }

    /* ── Test 9: Capon has peaks at signal frequencies ────── */
    TEST_CASE_BEGIN("Capon spectrum has peaks at signal frequencies");
    {
        double f1 = 0.15;
        int N = 256;
        double *x = (double *)malloc((size_t)N * sizeof(double));
        unsigned int seed_c = 55;
        for (int i = 0; i < N; i++)
            x[i] = sin(2.0 * M_PI * f1 * i) + test_randn(&seed_c) * 0.01;

        int nfft = 512;
        double *spec = (double *)malloc((size_t)(nfft / 2) * sizeof(double));
        capon_spectrum(x, N, 16, spec, nfft);

        /* Find peak */
        int peak_bin = 0;
        double peak_val = -1e30;
        for (int i = 1; i < nfft / 2; i++)
            if (spec[i] > peak_val) { peak_val = spec[i]; peak_bin = i; }

        double f_est = (double)peak_bin / (double)nfft;
        if (fabs(f_est - f1) < 0.02) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("Capon peak not near expected freq"); }

        free(x); free(spec);
    }

    /* ── Test 10: Correlation matrix is symmetric ─────────── */
    TEST_CASE_BEGIN("Correlation matrix is symmetric");
    {
        int N = 128, p = 8;
        double *x = (double *)malloc((size_t)N * sizeof(double));
        unsigned int seed = 33;
        for (int i = 0; i < N; i++) x[i] = test_randn(&seed);

        double *R = (double *)malloc((size_t)(p * p) * sizeof(double));
        correlation_matrix(x, N, p, R);

        double max_asym = 0;
        for (int i = 0; i < p; i++)
            for (int j = i + 1; j < p; j++) {
                double d = fabs(R[i * p + j] - R[j * p + i]);
                if (d > max_asym) max_asym = d;
            }

        if (max_asym < 1e-14) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("R not symmetric"); }

        free(x); free(R);
    }

    /* ── Test 11: Real cepstrum of impulse ────────────────── */
    TEST_CASE_BEGIN("Real cepstrum of impulse is near-zero (>0)");
    {
        int nfft = 256;
        double *x = (double *)calloc((size_t)nfft, sizeof(double));
        x[0] = 1.0;  /* impulse → flat spectrum → log(1)=0 → IFFT=0 */

        double *c = (double *)malloc((size_t)nfft * sizeof(double));
        cepstrum_real(x, nfft, c, nfft);

        /* c[0] should be ~0 (log|1|=0), rest should be ~0 */
        double max_val = 0;
        for (int i = 1; i < nfft; i++)
            if (fabs(c[i]) > max_val) max_val = fabs(c[i]);

        if (max_val < 1e-10) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("Cepstrum of impulse not zero"); }

        free(x); free(c);
    }

    /* ── Test 12: Cepstral pitch detection ────────────────── */
    TEST_CASE_BEGIN("Cepstral pitch detection matches known period");
    {
        int nfft = 512, period = 50;
        double *x = (double *)calloc((size_t)nfft, sizeof(double));
        /* Impulse train with known period */
        for (int i = 0; i < nfft; i += period) x[i] = 1.0;

        double *c = (double *)malloc((size_t)nfft * sizeof(double));
        cepstrum_real(x, nfft, c, nfft);

        /* Find peak near expected period */
        int search_start = 20, search_end = 100;
        int peak_q = 0;
        double peak_val = 0;
        for (int q = search_start; q < search_end; q++)
            if (c[q] > peak_val) { peak_val = c[q]; peak_q = q; }

        if (abs(peak_q - period) <= 2) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("Cepstral pitch mismatch"); }

        free(x); free(c);
    }

    /* ── Test 13: Mel scale round-trip ────────────────────── */
    TEST_CASE_BEGIN("Mel scale round-trip: Hz → Mel → Hz");
    {
        double freqs[] = {100, 500, 1000, 2000, 4000};
        double max_err = 0;
        for (int i = 0; i < 5; i++) {
            double rt = mel_to_hz(hz_to_mel(freqs[i]));
            double err = fabs(rt - freqs[i]);
            if (err > max_err) max_err = err;
        }
        if (max_err < 1e-8) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("Mel round-trip error"); }
    }

    /* ── Test 14: MFCC output valid ───────────────────────── */
    TEST_CASE_BEGIN("MFCC output has valid values (no NaN/Inf)");
    {
        int frame_len = 256, nfft = 256;
        double *frame = (double *)malloc((size_t)frame_len * sizeof(double));
        for (int i = 0; i < frame_len; i++)
            frame[i] = sin(2.0 * M_PI * 0.1 * i);

        double mfcc[13];
        compute_mfcc(frame, frame_len, nfft, 8000.0, 26, 13, mfcc);

        int valid = 1;
        for (int i = 0; i < 13; i++)
            if (mfcc[i] != mfcc[i] || mfcc[i] == 1.0/0.0 || mfcc[i] == -1.0/0.0)
                valid = 0;

        if (valid) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("MFCC contains NaN or Inf"); }

        free(frame);
    }

    /* ── Test 15: DCT-II of constant ──────────────────────── */
    TEST_CASE_BEGIN("DCT-II of constant: non-zero only at k=0");
    {
        int N = 16;
        double *x = (double *)malloc((size_t)N * sizeof(double));
        double *y = (double *)malloc((size_t)N * sizeof(double));
        for (int i = 0; i < N; i++) x[i] = 3.0;

        dct_ii(x, y, N);

        /* y[0] should be 3.0 * N (or scaled), y[k>0] ≈ 0 */
        int ok = (fabs(y[0]) > 1e-6);
        for (int i = 1; i < N; i++)
            if (fabs(y[i]) > 1e-10) ok = 0;

        if (ok) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("DCT of constant not concentrated at k=0"); }

        free(x); free(y);
    }

    /* ── Test 16: Conv2D with identity kernel ─────────────── */
    TEST_CASE_BEGIN("Conv2D with identity kernel preserves image");
    {
        int rows = 8, cols = 8;
        double img[64], out[64];
        for (int i = 0; i < 64; i++) img[i] = (double)i / 64.0;

        double kernel[9] = {0,0,0, 0,1,0, 0,0,0};  /* 3×3 identity */
        conv2d(img, rows, cols, kernel, 3, 3, out);

        /* Interior pixels should match exactly */
        double max_err = 0;
        for (int r = 1; r < rows - 1; r++)
            for (int c = 1; c < cols - 1; c++) {
                double d = fabs(out[r * cols + c] - img[r * cols + c]);
                if (d > max_err) max_err = d;
            }
        if (max_err < 1e-14) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("Identity kernel changed interior"); }
    }

    /* ── Test 17: Gaussian kernel sums to 1.0 ─────────────── */
    TEST_CASE_BEGIN("Gaussian kernel sums to 1.0");
    {
        double kernel[49];
        kernel_gaussian(kernel, 7, 1.5);
        double sum = 0;
        for (int i = 0; i < 49; i++) sum += kernel[i];
        if (fabs(sum - 1.0) < 1e-10) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("Gaussian kernel sum != 1"); }
    }

    /* ── Test 18: Sobel of flat image returns zeros ────────── */
    TEST_CASE_BEGIN("Sobel of flat image returns zeros");
    {
        int rows = 16, cols = 16;
        double img[256], mag[256];
        for (int i = 0; i < 256; i++) img[i] = 0.5;

        sobel_magnitude(img, rows, cols, mag);

        /* Check only interior pixels (boundary has zero-pad effect) */
        double max_mag = 0;
        for (int r = 1; r < rows - 1; r++)
            for (int c = 1; c < cols - 1; c++) {
                double v = fabs(mag[r * cols + c]);
                if (v > max_mag) max_mag = v;
            }

        if (max_mag < 1e-10) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("Sobel of flat not zero"); }
    }

    /* ── Test 19: 2D FFT → IFFT round-trip ────────────────── */
    TEST_CASE_BEGIN("2D FFT/IFFT round-trip preserves data");
    {
        int rows = 8, cols = 8, N = 64;
        double re[64], im[64], orig[64];
        unsigned int seed = 77;
        for (int i = 0; i < N; i++) {
            re[i] = test_randn(&seed);
            orig[i] = re[i];
            im[i] = 0.0;
        }

        fft2d(re, im, rows, cols);
        ifft2d(re, im, rows, cols);

        double max_err = 0;
        for (int i = 0; i < N; i++) {
            double d = fabs(re[i] - orig[i]);
            if (d > max_err) max_err = d;
        }
        if (max_err < 1e-10) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("2D FFT round-trip error"); }
    }

    /* ── Summary ──────────────────────────────────────────── */
    printf("\n  ────────────────────────────\n");
    printf("  Results: %d/%d passed", test_passed, test_count);
    if (test_failed > 0)
        printf(", %d FAILED", test_failed);
    printf("\n\n");

    return test_failed > 0 ? 1 : 0;
}

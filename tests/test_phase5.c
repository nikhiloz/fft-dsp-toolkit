/**
 * @file test_phase5.c
 * @brief Unit tests for Phase 5 modules: multirate, hilbert, averaging, remez.
 *
 * Tests:
 *   1.  Decimate preserves low-frequency content
 *   2.  Interpolate→Decimate round-trip
 *   3.  Polyphase matches direct decimation
 *   4.  Hilbert FIR is antisymmetric & odd-only nonzero
 *   5.  Analytic signal envelope ≈ 1.0 for pure sine
 *   6.  Analytic signal FFT method matches FIR (steady state)
 *   7.  Envelope detects AM modulation
 *   8.  Instantaneous frequency of pure sine
 *   9.  Coherent averaging reduces noise
 *  10.  Moving average of constant = constant
 *  11.  Median filter removes single impulse
 *  12.  EMA tracks step input
 *  13.  Remez lowpass passband gain ≈ 0 dB
 *  14.  Remez lowpass stopband attenuation
 *  15.  Remez bandpass passband gain ≈ 0 dB
 *
 * Run: make test
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "test_framework.h"
#include "multirate.h"
#include "hilbert.h"
#include "averaging.h"
#include "remez.h"
#include "fft.h"
#include "filter.h"
#include "dsp_utils.h"
#include "signal_gen.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void)
{
    TEST_SUITE("Phase 5: Multirate, Hilbert, Averaging, Remez");

    /* ── Test 1: Decimate preserves low frequency ────────── */
    TEST_CASE_BEGIN("Decimate preserves low-frequency content");
    {
        int N = 512, M = 4;
        double fs = 8000.0;
        double *x = (double *)calloc((size_t)N, sizeof(double));
        gen_sine(x, N, 1.0, 100.0, fs, 0.0);  /* 100 Hz, well below 1000 Hz Nyquist */

        double *y = (double *)calloc((size_t)(N / M), sizeof(double));
        int out_len = decimate(x, N, M, y);

        /* Check: output energy should be significant (signal not destroyed) */
        double energy_in = 0.0, energy_out = 0.0;
        int skip = 10;
        for (int i = skip; i < N - skip; i++)
            energy_in += x[i] * x[i];
        for (int i = skip; i < out_len - skip; i++)
            energy_out += y[i] * y[i];
        /* Scale: energy_out at 1/M rate should be ~energy_in/M */
        double ratio = (energy_out * M) / (energy_in + 1e-20);

        if (ratio > 0.5 && ratio < 2.0) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("Decimated energy ratio out of range"); }
        free(x); free(y);
    }

    /* ── Test 2: Interpolate→Decimate round-trip ─────────── */
    TEST_CASE_BEGIN("Interpolate then decimate round-trip");
    {
        int N = 128, L = 4;
        double *x = (double *)calloc((size_t)N, sizeof(double));
        gen_sine(x, N, 1.0, 0.05, 1.0, 0.0);  /* very low freq */

        double *up = (double *)calloc((size_t)(N * L), sizeof(double));
        interpolate(x, N, L, up);

        double *down = (double *)calloc((size_t)N, sizeof(double));
        decimate(up, N * L, L, down);

        /* Should roughly recover original (with filter delay shift) */
        double energy_orig = 0.0, energy_out = 0.0;
        for (int i = 0; i < N; i++) {
            energy_orig += x[i] * x[i];
            energy_out  += down[i] * down[i];
        }

        double ratio = energy_out / (energy_orig + 1e-20);
        if (ratio > 0.3 && ratio < 3.0) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("Energy ratio out of range"); }
        free(x); free(up); free(down);
    }

    /* ── Test 3: Polyphase matches direct decimation ─────── */
    TEST_CASE_BEGIN("Polyphase matches direct decimation");
    {
        int N = 256, M = 4, taps = 33;
        double *x = (double *)calloc((size_t)N, sizeof(double));
        gen_sine(x, N, 1.0, 0.05, 1.0, 0.0);

        double *h = (double *)calloc((size_t)taps, sizeof(double));
        fir_lowpass(h, taps, 0.5 / M);

        /* Direct: filter then pick every M-th */
        double *filtered = (double *)calloc((size_t)N, sizeof(double));
        fir_filter(x, filtered, N, h, taps);
        int out_len = N / M;
        double *y_direct = (double *)calloc((size_t)out_len, sizeof(double));
        for (int i = 0; i < out_len; i++)
            y_direct[i] = filtered[i * M];

        /* Polyphase */
        double *y_poly = (double *)calloc((size_t)out_len, sizeof(double));
        polyphase_decimate(x, N, h, taps, M, y_poly);

        double max_err = 0.0;
        for (int i = 0; i < out_len; i++) {
            double err = fabs(y_direct[i] - y_poly[i]);
            if (err > max_err) max_err = err;
        }

        if (max_err < 1e-10) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("Polyphase != direct"); }

        free(x); free(h); free(filtered); free(y_direct); free(y_poly);
    }

    /* ── Test 4: Hilbert FIR symmetry ────────────────────── */
    TEST_CASE_BEGIN("Hilbert FIR antisymmetric, even coeffs zero");
    {
        int taps = 31;
        double *h = (double *)calloc((size_t)taps, sizeof(double));
        hilbert_design(h, taps);

        int mid = taps / 2;
        int ok = 1;
        /* h[mid] should be 0 (center of antisymmetric) */
        if (fabs(h[mid]) > 1e-12) ok = 0;
        /* Even-indexed offsets from centre should be 0 */
        for (int i = 0; i < mid && ok; i++) {
            int offset = mid - i;
            if (offset % 2 == 0 && fabs(h[i]) > 1e-12) ok = 0;
        }
        /* Antisymmetric: h[mid+k] = -h[mid-k] */
        for (int k = 1; k <= mid && ok; k++) {
            if (fabs(h[mid + k] + h[mid - k]) > 1e-12) ok = 0;
        }

        if (ok) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("Hilbert FIR symmetry violated"); }
        free(h);
    }

    /* ── Test 5: Analytic signal envelope ≈ 1.0 ──────────── */
    TEST_CASE_BEGIN("Analytic signal envelope ≈ 1.0 for pure sine");
    {
        int N = 256;
        double *x = (double *)calloc((size_t)N, sizeof(double));
        /* Use bin-aligned frequency: f = 8/256 = 0.03125 normalised */
        for (int i = 0; i < N; i++)
            x[i] = sin(2.0 * M_PI * 8.0 * (double)i / (double)N);

        Complex *z = (Complex *)calloc((size_t)N, sizeof(Complex));
        analytic_signal_fft(x, N, z);

        double max_err = 0.0;
        int skip = 2;
        for (int i = skip; i < N - skip; i++) {
            double env = complex_mag(z[i]);
            double err = fabs(env - 1.0);
            if (err > max_err) max_err = err;
        }

        if (max_err < 0.01) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("Envelope deviates from 1.0"); }
        free(x); free(z);
    }

    /* ── Test 6: FFT analytic matches FIR (steady state) ─── */
    TEST_CASE_BEGIN("Analytic signal: FFT vs FIR method");
    {
        int N = 256;
        double *x = (double *)calloc((size_t)N, sizeof(double));
        gen_sine(x, N, 1.0, 0.1, 1.0, 0.0);

        Complex *z_fir = (Complex *)calloc((size_t)N, sizeof(Complex));
        Complex *z_fft = (Complex *)calloc((size_t)N, sizeof(Complex));
        analytic_signal(x, N, z_fir, 51);
        analytic_signal_fft(x, N, z_fft);

        double max_env_diff = 0.0;
        int skip = 30;
        for (int i = skip; i < N - skip; i++) {
            double e1 = complex_mag(z_fir[i]);
            double e2 = complex_mag(z_fft[i]);
            double d = fabs(e1 - e2);
            if (d > max_env_diff) max_env_diff = d;
        }

        if (max_env_diff < 0.1) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("FIR vs FFT envelope differ too much"); }
        free(x); free(z_fir); free(z_fft);
    }

    /* ── Test 7: Envelope detects AM modulation ──────────── */
    TEST_CASE_BEGIN("Envelope detects AM modulation depth");
    {
        int N = 512;
        double fs = 4000.0;
        double *x = (double *)calloc((size_t)N, sizeof(double));
        for (int i = 0; i < N; i++) {
            double t = (double)i / fs;
            x[i] = (1.0 + 0.5 * sin(2.0 * M_PI * 10.0 * t))
                  * cos(2.0 * M_PI * 500.0 * t);
        }

        double *env = (double *)calloc((size_t)N, sizeof(double));
        envelope(x, N, env, 0);

        /* Check envelope range: should be ~[0.5, 1.5] */
        double env_min = 1e10, env_max = -1e10;
        int skip = 20;
        for (int i = skip; i < N - skip; i++) {
            if (env[i] < env_min) env_min = env[i];
            if (env[i] > env_max) env_max = env[i];
        }

        int ok = (env_min > 0.3 && env_min < 0.7 &&
                  env_max > 1.3 && env_max < 1.7);
        if (ok) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("AM envelope range incorrect"); }
        free(x); free(env);
    }

    /* ── Test 8: Inst frequency of pure sine ─────────────── */
    TEST_CASE_BEGIN("Instantaneous frequency of pure sine");
    {
        int N = 256;
        double f_norm = 0.1; /* normalised frequency */
        double *x = (double *)calloc((size_t)N, sizeof(double));
        gen_sine(x, N, 1.0, f_norm, 1.0, 0.0);

        double *freq = (double *)calloc((size_t)N, sizeof(double));
        inst_frequency(x, N, freq, 0);

        double max_err = 0.0;
        int skip = 10;
        for (int i = skip; i < N - skip; i++) {
            double err = fabs(freq[i] - f_norm);
            if (err > max_err) max_err = err;
        }

        if (max_err < 0.02) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("IF of pure sine too far from expected"); }
        free(x); free(freq);
    }

    /* ── Test 9: Coherent averaging reduces noise ────────── */
    TEST_CASE_BEGIN("Coherent averaging SNR improvement");
    {
        int N = 256, K = 64;
        double *clean = (double *)calloc((size_t)N, sizeof(double));
        gen_sine(clean, N, 1.0, 0.05, 1.0, 0.0);

        unsigned int seed = 42;
        double **trials = (double **)malloc((size_t)K * sizeof(double *));
        for (int k = 0; k < K; k++) {
            trials[k] = (double *)malloc((size_t)N * sizeof(double));
            for (int i = 0; i < N; i++) {
                seed = seed * 1103515245 + 12345;
                double u1 = (double)seed / 4294967296.0;
                seed = seed * 1103515245 + 12345;
                double u2 = (double)seed / 4294967296.0;
                if (u1 < 1e-10) u1 = 1e-10;
                double noise = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
                trials[k][i] = clean[i] + 2.0 * noise;
            }
        }

        double *avg = (double *)calloc((size_t)N, sizeof(double));
        coherent_average((const double **)trials, K, N, avg);

        double snr_before, snr_after;
        compute_snr_improvement(trials[0], clean, avg, N,
                                &snr_before, &snr_after);

        /* With K=64, expect ~18 dB improvement (10*log10(64)) */
        double gain = snr_after - snr_before;
        if (gain > 10.0) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("Averaging gain < 10 dB with K=64"); }

        for (int k = 0; k < K; k++) free(trials[k]);
        free(trials); free(avg); free(clean);
    }

    /* ── Test 10: Moving average of constant ─────────────── */
    TEST_CASE_BEGIN("Moving average of constant = constant");
    {
        int N = 100, M = 7;
        double *x = (double *)malloc((size_t)N * sizeof(double));
        for (int i = 0; i < N; i++) x[i] = 3.14;

        double *y = (double *)calloc((size_t)N, sizeof(double));
        moving_average(x, N, M, y);

        int ok = 1;
        for (int i = M; i < N; i++) {
            if (fabs(y[i] - 3.14) > 1e-10) { ok = 0; break; }
        }

        if (ok) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("MA of constant should be constant"); }
        free(x); free(y);
    }

    /* ── Test 11: Median filter removes impulse ──────────── */
    TEST_CASE_BEGIN("Median filter removes single impulse");
    {
        int N = 64, M = 5;
        double *x = (double *)calloc((size_t)N, sizeof(double));
        x[32] = 100.0;  /* single impulse */

        double *y = (double *)calloc((size_t)N, sizeof(double));
        median_filter(x, N, M, y);

        /* Output at position 32 should be 0 (median of 4 zeros + 1 spike) */
        if (fabs(y[32]) < 1e-10) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("Median didn't remove impulse"); }
        free(x); free(y);
    }

    /* ── Test 12: EMA tracks step input ──────────────────── */
    TEST_CASE_BEGIN("EMA converges to step value");
    {
        int N = 200;
        double alpha = 0.1;
        double *x = (double *)calloc((size_t)N, sizeof(double));
        for (int i = 0; i < N; i++) x[i] = 1.0;  /* step at n=0 */

        double *y = (double *)calloc((size_t)N, sizeof(double));
        ema_filter(x, N, alpha, y);

        /* Should converge close to 1.0 by end */
        if (fabs(y[N - 1] - 1.0) < 0.001) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("EMA didn't converge to 1.0"); }
        free(x); free(y);
    }

    /* ── Test 13: Remez lowpass passband gain ─────────────── */
    TEST_CASE_BEGIN("Remez lowpass passband gain ≈ 0 dB");
    {
        int taps = 51;
        double *h = (double *)calloc((size_t)taps, sizeof(double));
        int ret = remez_lowpass(h, taps, 0.2, 0.3, 1.0, 1.0);

        if (ret != 0) {
            TEST_FAIL_STMT("Remez failed to converge");
        } else {
            /* Check DC gain (sum of coefficients) */
            double sum = 0.0;
            for (int i = 0; i < taps; i++) sum += h[i];
            double dc_db = 20.0 * log10(fabs(sum) + 1e-20);

            if (fabs(dc_db) < 1.0) { TEST_PASS_STMT; }
            else { TEST_FAIL_STMT("DC gain not near 0 dB"); }
        }
        free(h);
    }

    /* ── Test 14: Remez lowpass stopband attenuation ─────── */
    TEST_CASE_BEGIN("Remez lowpass stopband attenuation");
    {
        int taps = 51, nfft = 1024;
        double *h = (double *)calloc((size_t)taps, sizeof(double));
        remez_lowpass(h, taps, 0.2, 0.3, 1.0, 1.0);

        Complex *H = (Complex *)calloc((size_t)nfft, sizeof(Complex));
        for (int i = 0; i < taps; i++) { H[i].re = h[i]; H[i].im = 0.0; }
        fft(H, nfft);

        /* Stopband: f=0.3 to f=0.5 (normalised, Nyquist=0.5)
         * bin k = f * nfft for unit sample rate analysis */
        int stop_start = (int)(0.3 * (double)nfft);
        int stop_end   = nfft / 2;
        double max_stop = -200.0;
        for (int i = stop_start; i <= stop_end; i++) {
            double mag_db = 20.0 * log10(complex_mag(H[i]) + 1e-20);
            if (mag_db > max_stop) max_stop = mag_db;
        }

        /* Expect at least 15 dB stopband attenuation */
        if (max_stop < -15.0) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("Stopband attenuation < 15 dB"); }
        free(h); free(H);
    }

    /* ── Test 15: Remez bandpass passband ─────────────────── */
    TEST_CASE_BEGIN("Remez bandpass passband gain ≈ 0 dB");
    {
        int taps = 61, nfft = 1024;
        double *h = (double *)calloc((size_t)taps, sizeof(double));
        int ret = remez_bandpass(h, taps, 0.1, 0.2, 0.35, 0.45);

        if (ret != 0) {
            TEST_FAIL_STMT("Remez bandpass failed to converge");
        } else {
            Complex *H = (Complex *)calloc((size_t)nfft, sizeof(Complex));
            for (int i = 0; i < taps; i++) { H[i].re = h[i]; H[i].im = 0.0; }
            fft(H, nfft);

            /* Passband centre: f = 0.275 normalised → bin = 0.275 * nfft */
            int pass_bin = (int)(0.275 * (double)nfft);
            double gain_db = 20.0 * log10(complex_mag(H[pass_bin]) + 1e-20);

            if (fabs(gain_db) < 3.0) { TEST_PASS_STMT; }
            else { TEST_FAIL_STMT("Passband gain not near 0 dB"); }
            free(H);
        }
        free(h);
    }

    printf("\n=== Test Summary ===\n");
    printf("Total: %d, Passed: %d, Failed: %d\n",
           test_count, test_passed, test_failed);
    printf("Pass Rate: %.1f%%\n",
           test_count > 0 ? (100.0 * test_passed / test_count) : 0.0);
    return (test_failed == 0) ? 0 : 1;
}

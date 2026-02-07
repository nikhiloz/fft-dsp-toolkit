/**
 * @file 26-cepstrum-mfcc.c
 * @brief Chapter 26 demo — Cepstrum Analysis & MFCCs.
 *
 * ── Cepstral Domain ─────────────────────────────────────────────
 *
 *   Time domain:    x[n] = s[n] * h[n]   (convolution)
 *   Log spectrum:   log|X| = log|S| + log|H|  (addition)
 *   Cepstrum:       c[n] = cₛ[n] + cₕ[n]     (addition in quefrency)
 *
 *      low quefrency ← spectral envelope (vocal tract)
 *     high quefrency ← pitch (glottal excitation)
 *
 *           cepstrum
 *    ┌──────────────────────┐
 *    │ ↑                    │
 *    │ │ env  pitch         │
 *    │ │  ↓    ↓            │
 *    │ ·····│·····│·····    │
 *    └──────────────────────┘
 *      0   L-1   quefrency →
 *
 * ── MFCC Pipeline ────────────────────────────────────────────────
 *
 *   frame → [Hamming] → [FFT] → [|·|²] → [Mel Fbank] → [log] → [DCT]
 *
 * Demonstrates:
 *   - Real cepstrum computation
 *   - Spectral envelope via liftering
 *   - Mel filterbank visualisation
 *   - MFCC extraction
 *   - Pitch detection via cepstral peak
 *
 * Build & run:
 *   make chapters && ./build/bin/ch26
 *
 * Read alongside: chapters/26-cepstrum-mfcc.md
 */

#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "cepstrum.h"
#include "fft.h"
#include "dsp_utils.h"
#include "gnuplot.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define FRAME_LEN 512
#define NFFT      512

/* ── Demo 1: Real Cepstrum and Liftering ──────────────────────── */

static void demo_cepstrum_lifter(void)
{
    printf("=== Demo 1: Real Cepstrum & Spectral Envelope ===\n\n");

    /* Simulate a "speech-like" signal: impulse train through resonant filter */
    double x[FRAME_LEN];
    memset(x, 0, sizeof(x));
    int pitch_period = 80;  /* ~100 Hz at 8000 Hz sample rate */
    for (int i = 0; i < FRAME_LEN; i += pitch_period)
        x[i] = 1.0;

    /* Simple resonance: IIR y[n] = x[n] + 1.5 y[n-1] - 0.85 y[n-2] */
    double y[FRAME_LEN];
    y[0] = x[0]; y[1] = x[1] + 1.5 * y[0];
    for (int i = 2; i < FRAME_LEN; i++)
        y[i] = x[i] + 1.5 * y[i - 1] - 0.85 * y[i - 2];

    /* Normalise */
    double mx = 0;
    for (int i = 0; i < FRAME_LEN; i++)
        if (fabs(y[i]) > mx) mx = fabs(y[i]);
    for (int i = 0; i < FRAME_LEN; i++) y[i] /= mx;

    /* Compute real cepstrum */
    double cepstrum[NFFT];
    cepstrum_real(y, FRAME_LEN, cepstrum, NFFT);

    /* Detect pitch from cepstral peak */
    int min_q = 20, max_q = 200;
    double max_val = 0.0;
    int pitch_q = 0;
    for (int q = min_q; q < max_q && q < NFFT / 2; q++) {
        if (cepstrum[q] > max_val) {
            max_val = cepstrum[q];
            pitch_q = q;
        }
    }
    printf("  Pitch period (cepstral peak): %d samples\n", pitch_q);
    printf("  Expected pitch period:        %d samples\n", pitch_period);
    printf("  At fs=8000 Hz → f0 ≈ %d Hz\n\n", 8000 / pitch_q);

    /* Lifter: extract spectral envelope */
    int L = 30;
    double envelope[NFFT / 2];
    cepstrum_lifter(cepstrum, NFFT, L, envelope);

    /* Original log spectrum */
    Complex *X = (Complex *)calloc((size_t)NFFT, sizeof(Complex));
    for (int i = 0; i < FRAME_LEN; i++) {
        X[i].re = y[i];
        X[i].im = 0.0;
    }
    fft(X, NFFT);
    double log_spec[NFFT / 2];
    for (int i = 0; i < NFFT / 2; i++)
        log_spec[i] = 20.0 * log10(complex_mag(X[i]) + 1e-30);

    FILE *gp = popen("gnuplot -persistent", "w");
    if (gp) {
        fprintf(gp, "set terminal pngcairo size 900,600\n");
        fprintf(gp, "set output 'plots/ch26/cepstrum_lifter.png'\n");
        fprintf(gp, "set multiplot layout 2,1 title"
                     " 'Cepstral Analysis & Spectral Envelope'\n");
        fprintf(gp, "set grid\n");

        fprintf(gp, "set title 'Real Cepstrum'\n");
        fprintf(gp, "set xlabel 'Quefrency (samples)'\nset ylabel 'Amplitude'\n");
        fprintf(gp, "set xrange [0:200]\n");
        fprintf(gp, "plot '-' w l lw 1 lc rgb '#0077BB' notitle\n");
        for (int i = 0; i < 200; i++)
            fprintf(gp, "%d %f\n", i, cepstrum[i]);
        fprintf(gp, "e\n");

        fprintf(gp, "set title 'Log Spectrum + Envelope (L=%d)'\n", L);
        fprintf(gp, "set xlabel 'Frequency Bin'\nset ylabel 'dB'\n");
        fprintf(gp, "set xrange [0:256]\nset auto y\n");
        fprintf(gp, "plot '-' w l lw 1 lc rgb '#AAAAAA' title 'Log|X|',"
                     " '-' w l lw 3 lc rgb '#CC3311' title 'Envelope'\n");
        for (int i = 0; i < NFFT / 2; i++)
            fprintf(gp, "%d %f\n", i, log_spec[i]);
        fprintf(gp, "e\n");
        for (int i = 0; i < NFFT / 2; i++)
            fprintf(gp, "%d %f\n", i, envelope[i]);
        fprintf(gp, "e\nunset multiplot\n");
        pclose(gp);
        printf("  Plot → plots/ch26/cepstrum_lifter.png\n\n");
    }
    free(X);
}

/* ── Demo 2: Mel Filterbank ───────────────────────────────────── */

static void demo_mel_filterbank(void)
{
    printf("=== Demo 2: Mel Filterbank ===\n\n");

    double fs = 8000.0;
    int n_filters = 26;

    printf("  Mel scale mapping (selected points):\n");
    double freqs[] = {0, 200, 500, 1000, 2000, 3000, 4000};
    for (int i = 0; i < 7; i++)
        printf("    %5.0f Hz → %6.1f Mel\n", freqs[i], hz_to_mel(freqs[i]));
    printf("\n");

    /* Apply each filter individually */
    FILE *gp = popen("gnuplot -persistent", "w");
    if (gp) {
        fprintf(gp, "set terminal pngcairo size 800,400\n");
        fprintf(gp, "set output 'plots/ch26/mel_filterbank.png'\n");
        fprintf(gp, "set title 'Mel Filterbank (%d filters, fs=%.0f Hz)'\n",
                n_filters, fs);
        fprintf(gp, "set xlabel 'Frequency (Hz)'\nset ylabel 'Weight'\n");
        fprintf(gp, "set grid\n");

        /* Reconstruct triangular filters manually for plotting */
        double mel_low = hz_to_mel(0.0);
        double mel_high = hz_to_mel(fs / 2.0);
        int n_pts = n_filters + 2;
        double *mel_pts = (double *)malloc((size_t)n_pts * sizeof(double));
        double *hz_pts = (double *)malloc((size_t)n_pts * sizeof(double));

        for (int i = 0; i < n_pts; i++) {
            mel_pts[i] = mel_low + (mel_high - mel_low) * i / (n_pts - 1);
            hz_pts[i] = mel_to_hz(mel_pts[i]);
        }

        fprintf(gp, "plot ");
        for (int m = 0; m < n_filters; m++) {
            if (m > 0) fprintf(gp, ", ");
            fprintf(gp, "'-' w l lw 1 notitle");
        }
        fprintf(gp, "\n");

        for (int m = 0; m < n_filters; m++) {
            fprintf(gp, "%f 0\n", hz_pts[m]);
            fprintf(gp, "%f 1\n", hz_pts[m + 1]);
            fprintf(gp, "%f 0\n", hz_pts[m + 2]);
            fprintf(gp, "e\n");
        }
        pclose(gp);
        printf("  Plot → plots/ch26/mel_filterbank.png\n\n");

        free(mel_pts); free(hz_pts);
    }
}

/* ── Demo 3: MFCC Extraction ─────────────────────────────────── */

static void demo_mfcc(void)
{
    printf("=== Demo 3: MFCC Extraction ===\n\n");

    double fs = 8000.0;
    int n_filters = 26;
    int n_mfcc = 13;

    /* Generate two different "signals" */
    double frame_a[FRAME_LEN], frame_b[FRAME_LEN];
    for (int i = 0; i < FRAME_LEN; i++) {
        /* Signal A: low frequency */
        frame_a[i] = sin(2.0 * M_PI * 200.0 * i / fs)
                   + 0.5 * sin(2.0 * M_PI * 400.0 * i / fs);
        /* Signal B: higher frequency */
        frame_b[i] = sin(2.0 * M_PI * 1000.0 * i / fs)
                   + 0.7 * sin(2.0 * M_PI * 2000.0 * i / fs);
    }

    double mfcc_a[13], mfcc_b[13];
    compute_mfcc(frame_a, FRAME_LEN, NFFT, fs, n_filters, n_mfcc, mfcc_a);
    compute_mfcc(frame_b, FRAME_LEN, NFFT, fs, n_filters, n_mfcc, mfcc_b);

    printf("  MFCC coefficients:\n");
    printf("  Index   Signal A     Signal B\n");
    printf("  -----   ----------   ----------\n");
    for (int i = 0; i < n_mfcc; i++)
        printf("    %2d    %+9.3f    %+9.3f\n", i, mfcc_a[i], mfcc_b[i]);

    /* Euclidean distance */
    double dist = 0.0;
    for (int i = 0; i < n_mfcc; i++)
        dist += (mfcc_a[i] - mfcc_b[i]) * (mfcc_a[i] - mfcc_b[i]);
    dist = sqrt(dist);
    printf("\n  MFCC Euclidean distance: %.3f\n\n", dist);

    /* Bar chart */
    FILE *gp = popen("gnuplot -persistent", "w");
    if (gp) {
        fprintf(gp, "set terminal pngcairo size 800,400\n");
        fprintf(gp, "set output 'plots/ch26/mfcc_comparison.png'\n");
        fprintf(gp, "set title 'MFCC Comparison (13 coefficients)'\n");
        fprintf(gp, "set xlabel 'MFCC Index'\nset ylabel 'Value'\n");
        fprintf(gp, "set style data histogram\nset style histogram cluster gap 1\n");
        fprintf(gp, "set style fill solid 0.7\nset grid ytics\n");
        fprintf(gp, "plot '-' using 2:xtic(1) title 'Signal A' lc rgb '#0077BB',"
                     " '' using 2:xtic(1) title 'Signal B' lc rgb '#CC3311'\n");
        for (int i = 0; i < n_mfcc; i++)
            fprintf(gp, "%d %f\n", i, mfcc_a[i]);
        fprintf(gp, "e\n");
        for (int i = 0; i < n_mfcc; i++)
            fprintf(gp, "%d %f\n", i, mfcc_b[i]);
        fprintf(gp, "e\n");
        pclose(gp);
        printf("  Plot → plots/ch26/mfcc_comparison.png\n\n");
    }
}

/* ── Main ─────────────────────────────────────────────────────── */

int main(void)
{
    printf("╔══════════════════════════════════════════════╗\n");
    printf("║  Chapter 26: Cepstrum Analysis & MFCCs      ║\n");
    printf("║  Cepstrum · Liftering · Mel · MFCC          ║\n");
    printf("╚══════════════════════════════════════════════╝\n\n");

    if (system("mkdir -p plots/ch26") != 0) { /* ignore */ }

    demo_cepstrum_lifter();
    demo_mel_filterbank();
    demo_mfcc();

    printf("All Chapter 26 demos complete.\n");
    return 0;
}

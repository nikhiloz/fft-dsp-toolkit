/**
 * @file 08-putting-it-together.c
 * @brief Chapter 8 demo — End-to-end pipeline combining FFT + filter.
 *
 * Demonstrates:
 *   - Generate a noisy multi-tone signal
 *   - Show the spectrum BEFORE filtering
 *   - Apply a lowpass filter
 *   - Show the spectrum AFTER filtering
 *   - Compare RMS and spectral content
 *
 * Build & run:
 *   make chapters && ./build/bin/ch08
 *
 * Read alongside: chapters/08-putting-it-together.md
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "fft.h"
#include "filter.h"
#include "dsp_utils.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define N     512
#define FS    8000.0
#define TAPS  31

static void print_spectrum(const char *title, const double *sig, int n) {
    double windowed[512];
    Complex spectrum[512];
    double mag[512];

    memcpy(windowed, sig, (size_t)n * sizeof(double));
    apply_window(windowed, n, hann_window);
    fft_real(windowed, spectrum, n);
    fft_magnitude(spectrum, mag, n);

    printf("  %s:\n", title);
    printf("  Frequency   | dB        | Visual\n");
    printf("  ────────────┼───────────┼──────────\n");

    for (int k = 0; k < n / 2; k++) {
        double freq = (double)k * FS / n;
        double db = db_from_magnitude(mag[k] / (n / 2));
        if (db > -50.0) {
            /* Simple bar chart */
            int bar_len = (int)(db + 50) / 2;
            if (bar_len < 0) bar_len = 0;
            if (bar_len > 25) bar_len = 25;

            printf("  %7.1f Hz  | %+6.1f dB | ", freq, db);
            for (int b = 0; b < bar_len; b++) printf("█");
            printf("\n");
        }
    }
    printf("\n");
}

int main(void) {
    printf("=== Chapter 8: Putting It All Together ===\n\n");
    printf("Pipeline: signal → spectrum → filter → spectrum (compare)\n\n");

    /* ── Step 1: Generate signal with desired + unwanted components ─ */
    double clean[N], noisy[N], filtered[N];

    printf("── Step 1: Generate signal ──\n");
    printf("  Desired:  200 Hz (amplitude 1.0) + 500 Hz (amplitude 0.5)\n");
    printf("  Noise:    2800 Hz (0.3) + 3500 Hz (0.2)\n\n");

    for (int i = 0; i < N; i++) {
        double t = (double)i / FS;
        clean[i] = sin(2.0 * M_PI * 200.0 * t)
                 + 0.5 * sin(2.0 * M_PI * 500.0 * t);
        noisy[i] = clean[i]
                 + 0.3 * sin(2.0 * M_PI * 2800.0 * t)
                 + 0.2 * sin(2.0 * M_PI * 3500.0 * t);
    }

    /* ── Step 2: Show spectrum BEFORE filtering ──────────────────── */
    printf("── Step 2: Spectrum before filtering ──\n\n");
    print_spectrum("Noisy signal", noisy, N);

    /* ── Step 3: Design and apply lowpass filter ─────────────────── */
    printf("── Step 3: Apply %d-tap lowpass filter (cutoff 800 Hz) ──\n\n", TAPS);

    double h[TAPS];
    fir_lowpass(h, TAPS, 800.0 / FS);
    fir_filter(noisy, filtered, N, h, TAPS);

    /* ── Step 4: Show spectrum AFTER filtering ────────────────────── */
    printf("── Step 4: Spectrum after filtering ──\n\n");
    print_spectrum("Filtered signal (skip transient)", filtered + TAPS, N - TAPS >= 512 ? 512 : 256);

    /* ── Step 5: Quantitative comparison ─────────────────────────── */
    printf("── Step 5: Comparison ──\n\n");

    double rms_clean = rms(clean + TAPS, N - TAPS);
    double rms_noisy = rms(noisy + TAPS, N - TAPS);
    double rms_filt  = rms(filtered + TAPS, N - TAPS);

    double err_before = 0, err_after = 0;
    for (int i = TAPS; i < N; i++) {
        err_before += (noisy[i] - clean[i]) * (noisy[i] - clean[i]);
        err_after  += (filtered[i] - clean[i]) * (filtered[i] - clean[i]);
    }
    err_before = sqrt(err_before / (N - TAPS));
    err_after  = sqrt(err_after / (N - TAPS));

    printf("  Metric              | Before filter | After filter\n");
    printf("  ────────────────────┼───────────────┼─────────────\n");
    printf("  RMS signal          | %13.4f | %12.4f\n", rms_noisy, rms_filt);
    printf("  RMS error vs clean  | %13.4f | %12.4f\n", err_before, err_after);
    printf("  Clean reference RMS | %13.4f |\n", rms_clean);
    double improvement = 20.0 * log10(err_before / err_after);
    printf("\n  SNR improvement: %.1f dB\n", improvement);
    printf("  The lowpass filter successfully removed high-frequency noise.\n");

    return 0;
}

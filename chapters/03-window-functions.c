/**
 * @file 03-window-functions.c
 * @brief Chapter 3 demo — Compare window functions and their effect on spectra.
 *
 * Demonstrates:
 *   - Hann, Hamming, and Blackman window shapes
 *   - FFT of a windowed vs unwindowed signal
 *   - Side-lobe levels for each window
 *
 * Build & run:
 *   make chapters && ./build/bin/ch03
 *
 * Read alongside: chapters/03-window-functions.md
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "fft.h"
#include "dsp_utils.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define N   256
#define FS  8000.0

/* Run FFT and find peak + highest side-lobe within a range */
static void analyze_window(const char *name, window_fn win, double *base_signal) {
    double signal[N];
    Complex spectrum[N];
    double mag[N];

    memcpy(signal, base_signal, N * sizeof(double));

    if (win) {
        apply_window(signal, N, win);
    }

    fft_real(signal, spectrum, N);
    fft_magnitude(spectrum, mag, N);

    /* Find peak bin */
    int peak_bin = 0;
    double peak_mag = 0;
    for (int k = 1; k < N / 2; k++) {
        if (mag[k] > peak_mag) {
            peak_mag = mag[k];
            peak_bin = k;
        }
    }

    double peak_db = db_from_magnitude(peak_mag);

    /* Find highest side-lobe (bins more than 3 away from peak) */
    double max_sidelobe = 0;
    for (int k = 1; k < N / 2; k++) {
        if (abs(k - peak_bin) > 3 && mag[k] > max_sidelobe) {
            max_sidelobe = mag[k];
        }
    }
    double sidelobe_db = db_from_magnitude(max_sidelobe);

    printf("  %-12s  peak at bin %2d (%6.1f Hz)  |  peak: %+6.1f dB  |  "
           "side-lobe: %+6.1f dB  |  suppression: %.0f dB\n",
           name, peak_bin, (double)peak_bin * FS / N,
           peak_db, sidelobe_db, peak_db - sidelobe_db);
}

int main(void) {
    printf("=== Chapter 3: Window Functions ===\n\n");

    /* ── Part 1: Print window shapes ─────────────────────────────── */
    printf("── Window shapes (N=16, showing w[i]) ──\n\n");
    printf("  i   | Rectangular | Hann    | Hamming | Blackman\n");
    printf("  ────┼────────────┼─────────┼─────────┼─────────\n");
    for (int i = 0; i < 16; i++) {
        printf("  %2d  |   %5.3f    | %5.3f   | %5.3f   | %5.3f\n",
               i,
               1.0,
               hann_window(16, i),
               hamming_window(16, i),
               blackman_window(16, i));
    }

    /* ── Part 2: Spectral leakage comparison ─────────────────────── */
    printf("\n── Spectral leakage comparison ──\n");
    printf("  Signal: 440 Hz sine at fs=%g Hz, N=%d\n", FS, N);
    printf("  440/31.25 = 14.08 → falls BETWEEN bins → leakage expected\n\n");

    double signal[N];
    for (int i = 0; i < N; i++) {
        double t = (double)i / FS;
        signal[i] = sin(2.0 * M_PI * 440.0 * t);
    }

    analyze_window("Rectangular", NULL, signal);
    analyze_window("Hann", hann_window, signal);
    analyze_window("Hamming", hamming_window, signal);
    analyze_window("Blackman", blackman_window, signal);

    /* ── Part 3: A bin-centred frequency (no leakage) ────────────── */
    printf("\n── Control: bin-centred frequency (500 Hz = bin 16.0) ──\n");
    printf("  500/31.25 = 16.0 → falls EXACTLY on bin → no leakage\n\n");

    for (int i = 0; i < N; i++) {
        double t = (double)i / FS;
        signal[i] = sin(2.0 * M_PI * 500.0 * t);
    }

    analyze_window("Rectangular", NULL, signal);
    analyze_window("Hann", hann_window, signal);

    printf("\n  When the frequency falls exactly on a bin, all windows\n");
    printf("  give a clean peak. Windows only matter for non-integer bins.\n");

    return 0;
}

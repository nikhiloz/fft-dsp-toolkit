/**
 * @file 05-spectral-analysis.c
 * @brief Chapter 5 demo — Full spectral analysis pipeline.
 *
 * Demonstrates:
 *   - Signal → window → FFT → magnitude → dB → display
 *   - Frequency resolution and bin interpretation
 *   - Effect of FFT size on spectral detail
 *   - Comparing windowed vs unwindowed spectra
 *
 * Build & run:
 *   make chapters && ./build/bin/ch05
 *
 * Read alongside: chapters/05-spectral-analysis.md
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "fft.h"
#include "dsp_utils.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define FS  8000.0

static void run_analysis(const char *title, const double *raw, int n, window_fn win) {
    double signal[1024];
    Complex spectrum[1024];
    double mag[1024];

    memcpy(signal, raw, (size_t)n * sizeof(double));
    if (win) apply_window(signal, n, win);
    fft_real(signal, spectrum, n);
    fft_magnitude(spectrum, mag, n);

    printf("  %s (N=%d, resolution=%.1f Hz/bin):\n", title, n, FS / n);
    printf("  Frequency   | dB\n");
    printf("  ────────────┼──────\n");
    int printed = 0;
    for (int k = 0; k < n / 2; k++) {
        double freq = (double)k * FS / n;
        double db = db_from_magnitude(mag[k] / (n / 2));
        if (db > -40.0) {
            printf("  %7.1f Hz  | %+6.1f dB", freq, db);
            if (db > -15.0) printf("  ◄── peak");
            printf("\n");
            printed++;
        }
    }
    if (printed == 0) printf("  (no bins above -40 dB threshold)\n");
    printf("\n");
}

int main(void) {
    printf("=== Chapter 5: Spectral Analysis ===\n\n");

    /* ── Generate test signal: 440 Hz + 1000 Hz + 2500 Hz ───────── */
    double signal_256[256], signal_512[512];

    for (int i = 0; i < 512; i++) {
        double t = (double)i / FS;
        double s = 1.0 * sin(2.0 * M_PI * 440.0 * t)
                 + 0.5 * sin(2.0 * M_PI * 1000.0 * t)
                 + 0.3 * sin(2.0 * M_PI * 2500.0 * t);
        if (i < 256) signal_256[i] = s;
        signal_512[i] = s;
    }

    /* ── Part 1: No window vs Hann window ────────────────────────── */
    printf("── Part 1: Windowing effect (N=256) ──\n\n");
    run_analysis("No window (rectangular)", signal_256, 256, NULL);
    run_analysis("Hann window", signal_256, 256, hann_window);

    printf("  With the Hann window, peaks are narrower and side lobes\n");
    printf("  are suppressed — cleaner spectrum at the cost of slightly\n");
    printf("  wider main lobes.\n\n");

    /* ── Part 2: FFT size comparison ─────────────────────────────── */
    printf("── Part 2: N=256 vs N=512 (Hann window) ──\n\n");
    run_analysis("N=256 (31.25 Hz/bin)", signal_256, 256, hann_window);
    run_analysis("N=512 (15.63 Hz/bin)", signal_512, 512, hann_window);

    printf("  Doubling N halves the bin width → peaks are better resolved.\n");
    printf("  Trade-off: need 2x more samples and 2x more computation.\n\n");

    /* ── Part 3: RMS comparison ──────────────────────────────────── */
    printf("── Part 3: Signal statistics ──\n\n");
    printf("  RMS (256 samples): %.4f\n", rms(signal_256, 256));
    printf("  RMS (512 samples): %.4f\n", rms(signal_512, 512));
    printf("  Nyquist frequency: %.0f Hz (fs/2)\n", FS / 2.0);
    printf("  All 3 frequencies (440, 1000, 2500 Hz) are below Nyquist → no aliasing.\n");

    return 0;
}

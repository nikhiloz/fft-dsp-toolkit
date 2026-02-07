/**
 * @file 02-fft-fundamentals.c
 * @brief Chapter 2 demo — Step-by-step FFT on small signals.
 *
 * Demonstrates:
 *   - 8-point FFT on known signals (impulse, DC, sine)
 *   - Frequency bin interpretation
 *   - FFT → IFFT roundtrip
 *   - Dual-tone spectrum analysis (440 Hz + 1000 Hz)
 *
 * Build & run:
 *   make chapters && ./build/bin/ch02
 *
 * Read alongside: chapters/02-fft-fundamentals.md
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "fft.h"
#include "dsp_utils.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static void print_spectrum(const char *title, Complex *x, int n) {
    printf("  %s:\n", title);
    for (int k = 0; k < n; k++) {
        double mag = complex_mag(x[k]);
        printf("    bin[%d]  %+7.3f %+7.3fi   |X| = %.3f\n",
               k, x[k].re, x[k].im, mag);
    }
    printf("\n");
}

int main(void) {
    printf("=== Chapter 2: FFT Fundamentals ===\n\n");

    /* ── Demo 1: Impulse → flat spectrum ─────────────────────────── */
    printf("── Demo 1: Impulse signal (delta function) ──\n");
    Complex impulse[8] = { {1,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0} };
    fft(impulse, 8);
    print_spectrum("FFT of [1,0,0,0,0,0,0,0]", impulse, 8);
    printf("  → All bins have magnitude 1.0 (flat spectrum)\n");
    printf("  → An impulse contains ALL frequencies equally.\n\n");

    /* ── Demo 2: DC signal → energy in bin 0 ─────────────────────── */
    printf("── Demo 2: DC (constant) signal ──\n");
    Complex dc[8];
    for (int i = 0; i < 8; i++) dc[i] = (Complex){1.0, 0.0};
    fft(dc, 8);
    print_spectrum("FFT of [1,1,1,1,1,1,1,1]", dc, 8);
    printf("  → Only bin 0 has energy (magnitude 8 = N).\n");
    printf("  → DC = zero frequency.\n\n");

    /* ── Demo 3: Alternating → Nyquist ───────────────────────────── */
    printf("── Demo 3: Alternating signal [1,-1,1,-1,...] ──\n");
    Complex alt[8];
    for (int i = 0; i < 8; i++) alt[i] = (Complex){ (i % 2 == 0) ? 1.0 : -1.0, 0.0 };
    fft(alt, 8);
    print_spectrum("FFT of [1,-1,1,-1,1,-1,1,-1]", alt, 8);
    printf("  → Only bin N/2 (bin 4) has energy.\n");
    printf("  → This is the Nyquist frequency (highest representable).\n\n");

    /* ── Demo 4: FFT ↔ IFFT roundtrip ────────────────────────────── */
    printf("── Demo 4: FFT then IFFT recovers original ──\n");
    Complex original[8], roundtrip[8];
    for (int i = 0; i < 8; i++) {
        original[i] = (Complex){ sin(2.0 * M_PI * i / 8.0), 0.0 };
        roundtrip[i] = original[i];
    }
    fft(roundtrip, 8);
    printf("  After FFT:\n");
    for (int i = 0; i < 8; i++)
        printf("    X[%d] = %+.3f %+.3fi\n", i, roundtrip[i].re, roundtrip[i].im);
    ifft(roundtrip, 8);
    printf("  After IFFT:\n");
    double max_err = 0;
    for (int i = 0; i < 8; i++) {
        double err = fabs(roundtrip[i].re - original[i].re);
        if (err > max_err) max_err = err;
        printf("    x[%d] = %+.6f  (original: %+.6f,  error: %.1e)\n",
               i, roundtrip[i].re, original[i].re, err);
    }
    printf("  Max roundtrip error: %.1e\n\n", max_err);

    /* ── Demo 5: Real-world spectrum analysis ────────────────────── */
    printf("── Demo 5: 256-point FFT of 440 Hz + 1000 Hz ──\n");
    #define N   256
    #define FS  8000.0
    double signal[N];
    Complex spectrum[N];
    double mag[N];

    for (int i = 0; i < N; i++) {
        double t = (double)i / FS;
        signal[i] = sin(2.0 * M_PI * 440.0 * t) + 0.5 * sin(2.0 * M_PI * 1000.0 * t);
    }
    apply_window(signal, N, hann_window);
    fft_real(signal, spectrum, N);
    fft_magnitude(spectrum, mag, N);

    printf("  Frequency   | Magnitude (dB)\n");
    printf("  ────────────┼───────────────\n");
    for (int k = 0; k < N / 2; k++) {
        double freq = (double)k * FS / N;
        double db = db_from_magnitude(mag[k] / (N / 2));
        if (db > -40.0) {
            printf("  %7.1f Hz  | %+6.1f dB", freq, db);
            if (fabs(freq - 440.0) < FS / N || fabs(freq - 1000.0) < FS / N)
                printf("  ◄── peak!");
            printf("\n");
        }
    }
    printf("  Resolution: %.1f Hz/bin\n", FS / N);

    return 0;
}

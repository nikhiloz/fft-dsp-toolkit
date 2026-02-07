/**
 * @file 04-digital-filters.c
 * @brief Chapter 4 demo — FIR filter design and noise reduction.
 *
 * Demonstrates:
 *   - Moving average filter on a step signal
 *   - Windowed-sinc lowpass filter design
 *   - Filtering a noisy signal and measuring SNR improvement
 *   - Filter coefficient inspection
 *
 * Build & run:
 *   make chapters && ./build/bin/ch04
 *
 * Read alongside: chapters/04-digital-filters.md
 */

#include <stdio.h>
#include <math.h>
#include "filter.h"
#include "dsp_utils.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define N      256
#define FS     8000.0

int main(void) {
    printf("=== Chapter 4: Digital Filters ===\n\n");

    /* ── Part 1: Moving average on a step signal ─────────────────── */
    printf("── Part 1: 5-tap moving average on a step signal ──\n\n");

    double step_in[16], step_out[16];
    double ma_h[5];
    fir_moving_average(ma_h, 5);

    for (int i = 0; i < 16; i++) step_in[i] = (i >= 4) ? 1.0 : 0.0;
    fir_filter(step_in, step_out, 16, ma_h, 5);

    printf("  n  | input | output (5-pt avg) | note\n");
    printf("  ───┼───────┼───────────────────┼──────\n");
    for (int i = 0; i < 16; i++) {
        printf("  %2d |  %.1f  |      %.3f        ", i, step_in[i], step_out[i]);
        if (i >= 4 && i < 9)
            printf("| ← ramp (settling)");
        else if (i >= 9)
            printf("| ← settled to 1.0");
        printf("\n");
    }
    printf("\n  The moving average smoothly ramps from 0 to 1.\n");
    printf("  Settling time = %d samples (= filter order).\n\n", 5);

    /* ── Part 2: Lowpass filter coefficients ──────────────────────── */
    printf("── Part 2: 31-tap lowpass filter coefficients ──\n\n");

    #define TAPS 31
    double h[TAPS];
    double cutoff = 500.0 / FS;  /* 0.0625 normalized */
    fir_lowpass(h, TAPS, cutoff);

    printf("  Cutoff: 500 Hz (normalized: %.4f)\n", cutoff);
    printf("  Taps: %d (centre at tap %d)\n\n", TAPS, TAPS / 2);

    /* Print symmetry */
    printf("  Symmetry check (linear-phase FIR):\n");
    int symmetric = 1;
    for (int i = 0; i < TAPS / 2; i++) {
        if (fabs(h[i] - h[TAPS - 1 - i]) > 1e-12) symmetric = 0;
    }
    printf("  h[i] == h[M-1-i]?  %s\n\n", symmetric ? "YES ✓" : "NO ✗");

    double coeff_sum = 0;
    printf("  Coefficients (all %d):\n", TAPS);
    for (int i = 0; i < TAPS; i++) {
        printf("    h[%2d] = %+.6f", i, h[i]);
        if (i == TAPS / 2) printf("  ← centre (largest)");
        printf("\n");
        coeff_sum += h[i];
    }
    printf("  Sum of coefficients: %.6f  (should be 1.0 for unity DC gain)\n\n", coeff_sum);

    /* ── Part 3: Noise reduction ─────────────────────────────────── */
    printf("── Part 3: Lowpass filtering a noisy signal ──\n\n");

    double clean[N], noisy[N], filtered[N];

    for (int i = 0; i < N; i++) {
        double t = (double)i / FS;
        clean[i] = sin(2.0 * M_PI * 200.0 * t);
        noisy[i] = clean[i]
                  + 0.3 * sin(2.0 * M_PI * 2800.0 * t)
                  + 0.2 * sin(2.0 * M_PI * 3500.0 * t);
    }

    fir_filter(noisy, filtered, N, h, TAPS);

    double rms_clean    = rms(clean, N);
    double rms_noisy    = rms(noisy, N);
    double rms_filt     = rms(filtered + TAPS, N - TAPS);

    double error_sum = 0;
    for (int i = TAPS; i < N; i++) {
        double e = filtered[i] - clean[i];
        error_sum += e * e;
    }
    double rms_error = sqrt(error_sum / (N - TAPS));

    printf("  Signal: 200 Hz sine + noise at 2800 Hz and 3500 Hz\n");
    printf("  Filter: %d-tap lowpass at 500 Hz\n\n", TAPS);
    printf("  Clean RMS:     %.4f\n", rms_clean);
    printf("  Noisy RMS:     %.4f  (noise added %.1f%%)\n",
           rms_noisy, (rms_noisy - rms_clean) / rms_clean * 100);
    printf("  Filtered RMS:  %.4f  (after %d-sample settling)\n", rms_filt, TAPS);
    printf("  Error vs clean: %.4f RMS\n", rms_error);
    printf("\n  The filter removed the high-frequency noise while\n");
    printf("  preserving the 200 Hz signal.\n");

    return 0;
}

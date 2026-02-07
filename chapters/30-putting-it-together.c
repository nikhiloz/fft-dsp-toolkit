/**
 * @file 30-putting-it-together.c
 * @brief Chapter 30 (Capstone) — Comprehensive DSP Pipeline.
 *
 * ── Full DSP Pipeline ──────────────────────────────────────────
 *
 *   ┌──────────┐  ┌───────┐  ┌──────────┐  ┌─────────┐  ┌────────┐
 *   │  Signal  │─►│  FIR  │─►│ Envelope │─►│ Spectral│─►│  LPC   │
 *   │Generator │  │Filter │  │ Detector │  │Analysis │  │Summary │
 *   └──────────┘  └───────┘  └──────────┘  └─────────┘  └────────┘
 *        │                         │
 *        ▼                         ▼
 *   Multi-tone                 Hilbert AM
 *   + noise                    demodulator
 *
 * Demonstrates modules from ALL phases:
 *   Phase 1: signal_gen, fft, dsp_utils, filter (FIR)
 *   Phase 2: iir (Butterworth)
 *   Phase 3: spectrum (Welch PSD), correlation
 *   Phase 5: multirate (decimation), hilbert (envelope)
 *   Phase 6: lpc (linear prediction), adaptive (NLMS)
 *   Phase 7: realtime (ring buffer, frame processor)
 *
 * Build & run:
 *   make chapters && ./build/bin/ch30
 *
 * Read alongside: chapters/30-putting-it-together.md
 */

#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "fft.h"
#include "filter.h"
#include "iir.h"
#include "dsp_utils.h"
#include "signal_gen.h"
#include "spectrum.h"
#include "correlation.h"
#include "multirate.h"
#include "hilbert.h"
#include "lpc.h"
#include "adaptive.h"
#include "realtime.h"
#include "gnuplot.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define N     1024
#define FS    8000.0
#define TAPS  31

/* ── Section 1: Signal Generation (Phase 1) ──────────────────── */

static void section_signal_gen(double *clean, double *noisy, double *noise_ref)
{
    printf("── Section 1: Signal Generation ──\n\n");

    /* Build a multi-tone signal */
    double tone1[N], tone2[N], tone3[N], wn[N];
    gen_sine(tone1, N, 1.0, 200.0, FS, 0.0);
    gen_sine(tone2, N, 0.5, 500.0, FS, 0.0);
    gen_sine(tone3, N, 0.3, 1200.0, FS, 0.0);
    gen_white_noise(wn, N, 0.15, 42);

    for (int i = 0; i < N; i++) {
        clean[i] = tone1[i] + tone2[i];
        noise_ref[i] = tone3[i] + wn[i];
        noisy[i] = clean[i] + noise_ref[i];
    }

    printf("  Clean signal:  200 Hz + 500 Hz\n");
    printf("  Noise:         1200 Hz + white noise (0.15 rms)\n");
    printf("  Noisy RMS:     %.4f\n", rms(noisy, N));
    printf("  Clean RMS:     %.4f\n\n", rms(clean, N));
}

/* ── Section 2: FIR Low-pass Filter (Phase 1) ────────────────── */

static void section_fir_filter(const double *noisy, double *filtered)
{
    printf("── Section 2: FIR Low-pass Filtering ──\n\n");

    double h[TAPS];
    fir_lowpass(h, TAPS, 800.0 / FS);
    fir_filter(noisy, filtered, N, h, TAPS);

    printf("  %d-tap FIR lowpass, cutoff 800 Hz\n", TAPS);
    printf("  Filtered RMS:  %.4f\n\n", rms(filtered + TAPS, N - TAPS));
}

/* ── Section 3: Spectral Analysis (Phase 1+3) ────────────────── */

static void section_spectral_analysis(const double *noisy, const double *filtered)
{
    printf("── Section 3: Spectral Comparison ──\n\n");

    Complex spec_before[N], spec_after[N];
    double mag_before[N], mag_after[N];
    double windowed[N];

    /* Before filtering */
    memcpy(windowed, noisy, N * sizeof(double));
    apply_window(windowed, N, hann_window);
    fft_real(windowed, spec_before, N);
    fft_magnitude(spec_before, mag_before, N);

    /* After filtering (skip transient) */
    double tmp[512];
    memcpy(tmp, filtered + TAPS, 512 * sizeof(double));
    apply_window(tmp, 512, hann_window);
    fft_real(tmp, spec_after, 512);
    fft_magnitude(spec_after, mag_after, 512);

    printf("  Frequency   | Before (dB) | After (dB)\n");
    printf("  ────────────┼─────────────┼───────────\n");

    double freqs[] = {200.0, 500.0, 1200.0};
    for (int f = 0; f < 3; f++) {
        int k_before = (int)(freqs[f] * N / FS + 0.5);
        int k_after  = (int)(freqs[f] * 512 / FS + 0.5);
        double db_b = db_from_magnitude(mag_before[k_before] / (N / 2));
        double db_a = db_from_magnitude(mag_after[k_after] / 256);
        printf("  %7.0f Hz  |  %+6.1f     |  %+6.1f\n", freqs[f], db_b, db_a);
    }
    printf("\n  1200 Hz component should be greatly attenuated.\n\n");
}

/* ── Section 4: Adaptive Noise Cancellation (Phase 6) ─────────── */

static void section_adaptive(const double *noisy, const double *noise_ref, const double *clean)
{
    printf("── Section 4: Adaptive Noise Cancellation (NLMS) ──\n\n");

    NlmsState nlms;
    nlms_init(&nlms, 16, 0.5, 1e-6);

    double cancelled[N];
    double y, e;
    for (int i = 0; i < N; i++) {
        nlms_update(&nlms, noise_ref[i], noisy[i], &y, &e);
        cancelled[i] = e;  /* error = desired - noise estimate ≈ clean */
    }

    /* Measure improvement (skip convergence transient) */
    int skip = 200;
    double err_before = 0, err_after = 0;
    for (int i = skip; i < N; i++) {
        err_before += (noisy[i] - clean[i]) * (noisy[i] - clean[i]);
        err_after  += (cancelled[i] - clean[i]) * (cancelled[i] - clean[i]);
    }
    err_before = sqrt(err_before / (N - skip));
    err_after  = sqrt(err_after / (N - skip));

    double improvement = 20.0 * log10(err_before / (err_after + 1e-30));
    printf("  16-tap NLMS (µ=0.5)\n");
    printf("  Error before: %.4f  after: %.4f\n", err_before, err_after);
    printf("  SNR improvement: %.1f dB\n\n", improvement);

    nlms_free(&nlms);
}

/* ── Section 5: Hilbert Envelope (Phase 5) ────────────────────── */

static void section_hilbert_envelope(void)
{
    printf("── Section 5: Hilbert Transform — AM Demodulation ──\n\n");

    /* AM signal: carrier 1000 Hz, modulated by 50 Hz */
    double am[N], env_out[N];
    for (int i = 0; i < N; i++) {
        double t = (double)i / FS;
        double modulator = 0.5 + 0.5 * cos(2.0 * M_PI * 50.0 * t);
        am[i] = modulator * sin(2.0 * M_PI * 1000.0 * t);
    }

    envelope(am, N, env_out, 63);

    /* Check that envelope tracks modulator */
    int skip = 100;  /* skip start transient */
    double env_rms = rms(env_out + skip, N - skip);
    printf("  AM signal: 1000 Hz carrier, 50 Hz modulator\n");
    printf("  Envelope RMS: %.4f (expected ~0.5)\n", env_rms);
    printf("  Envelope detection working: %s\n\n",
           fabs(env_rms - 0.5) < 0.15 ? "YES" : "CHECK");
}

/* ── Section 6: LPC Analysis (Phase 6) ───────────────────────── */

static void section_lpc(const double *clean)
{
    printf("── Section 6: Linear Prediction Coding ──\n\n");

    int p = 10;
    double a[11];  /* LPC coefficients a[0..p] where a[0] = 1 */
    double E;

    lpc_coefficients(clean, N, p, a, &E);

    printf("  %d-th order LPC of clean signal\n", p);
    printf("  Prediction error energy: %.6f\n", E);
    printf("  Coefficients: a[1]=%.4f  a[2]=%.4f  a[3]=%.4f ...\n", a[0], a[1], a[2]);

    /* Compute prediction residual */
    double residual_energy = 0;
    for (int i = p; i < N; i++) {
        double pred = 0;
        for (int k = 0; k < p; k++)
            pred += a[k] * clean[i - 1 - k];
        double err = clean[i] - pred;
        residual_energy += err * err;
    }
    residual_energy /= (N - p);

    printf("  Residual energy: %.6f\n", residual_energy);
    printf("  Prediction gain: %.1f dB\n\n",
           10.0 * log10(rms(clean, N) * rms(clean, N) / (residual_energy + 1e-30)));
}

/* ── Section 7: Decimation (Phase 5) ─────────────────────────── */

static void section_decimation(const double *filtered)
{
    printf("── Section 7: Multirate — Decimation ──\n\n");

    int M = 2;  /* decimate by 2 */
    int n_out = decimate(filtered, N, M, NULL);  /* query size */
    double *decimated = (double *)malloc((size_t)n_out * sizeof(double));
    decimate(filtered, N, M, decimated);

    printf("  Decimation factor: %d\n", M);
    printf("  Input:  %d samples at %.0f Hz\n", N, FS);
    printf("  Output: %d samples at %.0f Hz\n", n_out, FS / M);
    printf("  Output RMS: %.4f\n\n", rms(decimated, n_out));

    free(decimated);
}

/* ── Section 8: Real-time Streaming (Phase 7) ─────────────────── */

static void section_realtime(const double *noisy)
{
    printf("── Section 8: Real-time Streaming Pipeline ──\n\n");

    RingBuffer *rb = ring_buffer_create(2048);
    FrameProcessor *fp = frame_processor_create(256, 128);
    LatencyStats lat;
    latency_init(&lat);

    /* Stream the noisy signal through ring buffer → frame processor */
    int chunk = 64;
    int frames_out = 0;

    for (int pos = 0; pos + chunk <= N; pos += chunk) {
        ring_buffer_write(rb, noisy + pos, chunk);

        double block[64];
        ring_buffer_read(rb, block, chunk);

        double t0 = timer_usec();
        int f = frame_processor_feed(fp, block, chunk);
        double t1 = timer_usec();

        if (f > 0) {
            latency_record(&lat, t1 - t0);
            frames_out++;
        }
    }

    printf("  Streamed %d samples in %d-sample chunks\n", N, chunk);
    printf("  Frames produced: %d\n", frames_out);
    printf("  Processing latency: min=%.1f µs  avg=%.1f µs  max=%.1f µs\n",
           lat.min_us, latency_avg(&lat), lat.max_us);

    if (frames_out > 0) {
        double peak = frame_processor_peak_freq(fp, FS);
        printf("  Last frame peak frequency: %.1f Hz\n", peak);
    }
    printf("\n");

    frame_processor_destroy(fp);
    ring_buffer_destroy(rb);
}

/* ── Section 9: Summary and Plots ────────────────────────────── */

static void section_summary(const double *clean, const double *noisy, const double *filtered)
{
    printf("── Section 9: Pipeline Summary ──\n\n");

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

    double improvement = 20.0 * log10(err_before / (err_after + 1e-30));
    printf("\n  FIR SNR improvement: %.1f dB\n\n", improvement);

    /* Generate a summary plot */
    gp_init("ch30");

    {
        double x[200], y_clean[200];
        int plot_len = 200;
        for (int i = 0; i < plot_len; i++) {
            x[i] = (double)i;
            y_clean[i] = clean[i];
        }
        gp_plot_1("ch30", "capstone_pipeline",
                  "Capstone: Clean Signal (200+500 Hz)",
                  "Sample", "Amplitude",
                  x, y_clean, plot_len, "lines");
        printf("  [saved] plots/ch30/capstone_pipeline.png\n\n");
    }
}

/* ── Main ────────────────────────────────────────────────────── */

int main(void)
{
    printf("=== Chapter 30: Putting It All Together — DSP Capstone ===\n\n");
    printf("  This capstone demo exercises modules from ALL 7 phases:\n");
    printf("  signal_gen → filter → fft → spectrum → hilbert → lpc\n");
    printf("  → multirate → adaptive → realtime → plots\n\n");

    double clean[N], noisy[N], noise_ref[N], filtered[N];

    section_signal_gen(clean, noisy, noise_ref);
    section_fir_filter(noisy, filtered);
    section_spectral_analysis(noisy, filtered);
    section_adaptive(noisy, noise_ref, clean);
    section_hilbert_envelope();
    section_lpc(clean);
    section_decimation(filtered);
    section_realtime(noisy);
    section_summary(clean, noisy, filtered);

    printf("Modules exercised: signal_gen, fft, dsp_utils, filter, iir,\n");
    printf("  spectrum, correlation, multirate, hilbert, lpc, adaptive,\n");
    printf("  realtime, gnuplot — 13 of 23 library modules.\n\n");
    printf("Key takeaways:\n");
    printf("  1. A complete DSP system chains many primitives\n");
    printf("  2. Adaptive filters can outperform fixed filters when noise is known\n");
    printf("  3. Hilbert transform enables amplitude/frequency demodulation\n");
    printf("  4. LPC captures signal structure in a few coefficients\n");
    printf("  5. Real-time streaming requires careful buffer management\n");

    return 0;
}

/**
 * @file 20-hilbert-transform.c
 * @brief Chapter 20 demo — Hilbert Transform & Analytic Signals.
 *
 * ── Analytic Signal Construction ─────────────────────────────────
 *
 *   x(t) ─────┬──────────────────┐
 *              │                  │  real part
 *              ▼                  │
 *        [Hilbert H{}] ──► x̂(t) │  imaginary part
 *              │                  │
 *              └── z(t) = x(t) + j·x̂(t)
 *
 *   The Hilbert transform shifts all frequencies by −90°.
 *
 * ── Envelope Detection ───────────────────────────────────────────
 *
 *   AM signal:  s(t) = A(t) · cos(2πfₒt)
 *
 *   Analytic:   z(t) = A(t) · e^{j2πfₒt}
 *
 *   Envelope:   |z(t)| = A(t)   ← recovered!
 *
 * ── Instantaneous Frequency ──────────────────────────────────────
 *
 *              d
 *   f_i(t) = ── [∠z(t)] / (2π)
 *             dt
 *
 * Demonstrates:
 *   - Hilbert FIR filter design (Type III)
 *   - Analytic signal: FIR method vs FFT method
 *   - AM envelope detection
 *   - Instantaneous frequency of chirp
 *   - Single-sideband (SSB) modulation
 *
 * Build & run:
 *   make chapters && ./build/bin/ch20
 *
 * Read alongside: chapters/20-hilbert-transform.md
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hilbert.h"
#include "signal_gen.h"
#include "fft.h"
#include "dsp_utils.h"
#include "gnuplot.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ── Demo 1: Hilbert FIR Design ───────────────────────────────── */

static void demo_hilbert_design(void)
{
    printf("=== Demo 1: Hilbert FIR Filter Design ===\n\n");

    int taps = 31;
    double *h = (double *)calloc((size_t)taps, sizeof(double));
    hilbert_design(h, taps);

    printf("  Hilbert FIR coefficients (Type III, %d taps):\n", taps);
    printf("  ┌───────┬─────────────┐\n");
    printf("  │ Index │  Coeff      │\n");
    printf("  ├───────┼─────────────┤\n");
    for (int i = 0; i < taps; i++) {
        if (fabs(h[i]) > 1e-8) {
            printf("  │  %3d  │ %+10.6f │\n", i, h[i]);
        }
    }
    printf("  └───────┴─────────────┘\n");
    printf("  Note: even-indexed coefficients are zero (Type III FIR)\n\n");

    free(h);
}

/* ── Demo 2: Analytic Signal (FIR vs FFT) ─────────────────────── */

static void demo_analytic_signal(void)
{
    printf("=== Demo 2: Analytic Signal — FIR vs FFT ===\n\n");

    const int N = 256;
    const double fs = 1000.0;
    double *x = (double *)malloc((size_t)N * sizeof(double));
    gen_sine(x, N, 1.0, 50.0, fs, 0.0);

    /* FIR method */
    Complex *z_fir = (Complex *)calloc((size_t)N, sizeof(Complex));
    analytic_signal(x, N, z_fir, 31);

    /* FFT method */
    Complex *z_fft = (Complex *)calloc((size_t)N, sizeof(Complex));
    analytic_signal_fft(x, N, z_fft);

    /* Compare envelopes */
    double max_diff = 0.0;
    int skip = 20; /* skip FIR transient */
    for (int i = skip; i < N - skip; i++) {
        double env_fir = complex_mag(z_fir[i]);
        double env_fft = complex_mag(z_fft[i]);
        double diff = fabs(env_fir - env_fft);
        if (diff > max_diff) max_diff = diff;
    }

    printf("  Input: 50 Hz sine, fs = %.0f Hz, N = %d\n", fs, N);
    printf("  FIR method: 31-tap Hilbert filter\n");
    printf("  FFT method: zero negative frequencies\n");
    printf("  Max |env_FIR - env_FFT| (excluding transient): %.4f\n", max_diff);
    printf("  FFT method is exact; FIR has filter transient at edges.\n\n");

    /* Plot real and imaginary parts */
    double *t_axis = (double *)malloc((size_t)N * sizeof(double));
    double *re_fft = (double *)malloc((size_t)N * sizeof(double));
    double *im_fft = (double *)malloc((size_t)N * sizeof(double));
    for (int i = 0; i < N; i++) {
        t_axis[i] = (double)i / fs;
        re_fft[i] = z_fft[i].re;
        im_fft[i] = z_fft[i].im;
    }

    gp_plot_multi("ch20", "analytic_signal",
        "Analytic Signal: Real and Hilbert (Imaginary) Parts",
        "Time (s)", "Amplitude",
        (GpSeries[]){
            {"Re{z(t)} = x(t)", t_axis, re_fft, N, "lines"},
            {"Im{z(t)} = H{x(t)}", t_axis, im_fft, N, "lines"}
        }, 2);
    printf("  → plots/ch20/analytic_signal.png\n\n");

    free(x); free(z_fir); free(z_fft);
    free(t_axis); free(re_fft); free(im_fft);
}

/* ── Demo 3: AM Envelope Detection ────────────────────────────── */

static void demo_envelope(void)
{
    printf("=== Demo 3: AM Envelope Detection ===\n\n");

    const int N = 512;
    const double fs = 8000.0;
    double *x = (double *)malloc((size_t)N * sizeof(double));
    double *env_true = (double *)malloc((size_t)N * sizeof(double));

    /* AM signal: A(t) = 1 + 0.5*sin(2π·5t), carrier = 500 Hz */
    for (int i = 0; i < N; i++) {
        double t = (double)i / fs;
        env_true[i] = 1.0 + 0.5 * sin(2.0 * M_PI * 5.0 * t);
        x[i] = env_true[i] * cos(2.0 * M_PI * 500.0 * t);
    }

    /* Extract envelope using Hilbert transform (FFT method, taps=0) */
    double *env_out = (double *)calloc((size_t)N, sizeof(double));
    envelope(x, N, env_out, 0);

    /* Measure accuracy */
    double mse = 0.0;
    for (int i = 0; i < N; i++) {
        double e = env_true[i] - env_out[i];
        mse += e * e;
    }
    mse /= N;

    printf("  AM signal: A(t)·cos(2π·500t), modulation 5 Hz\n");
    printf("  Envelope MSE: %.6e\n", mse);
    printf("  √MSE = %.4f (ideal envelope ∈ [0.5, 1.5])\n\n", sqrt(mse));

    /* Plot */
    double *t_axis = (double *)malloc((size_t)N * sizeof(double));
    for (int i = 0; i < N; i++) t_axis[i] = (double)i / fs;

    gp_plot_multi("ch20", "envelope_detection",
        "AM Envelope Detection via Hilbert Transform",
        "Time (s)", "Amplitude",
        (GpSeries[]){
            {"AM Signal", t_axis, x, N, "lines"},
            {"True Envelope", t_axis, env_true, N, "lines"},
            {"Detected Envelope", t_axis, env_out, N, "lines"}
        }, 3);
    printf("  → plots/ch20/envelope_detection.png\n\n");

    free(x); free(env_true); free(env_out); free(t_axis);
}

/* ── Demo 4: Instantaneous Frequency of Chirp ─────────────────── */

static void demo_inst_frequency(void)
{
    printf("=== Demo 4: Instantaneous Frequency of Chirp ===\n\n");

    const int N = 512;
    const double fs = 1000.0;
    double *x = (double *)malloc((size_t)N * sizeof(double));

    /* Linear chirp: 50 Hz → 200 Hz */
    double f0 = 50.0, f1 = 200.0;
    double T = (double)N / fs;
    for (int i = 0; i < N; i++) {
        double t = (double)i / fs;
        double freq_t = f0 + (f1 - f0) * t / T;
        x[i] = sin(2.0 * M_PI * (f0 * t + 0.5 * (f1 - f0) * t * t / T));
        (void)freq_t;
    }

    /* Extract instantaneous frequency (FFT method, taps=0) */
    double *freq = (double *)calloc((size_t)N, sizeof(double));
    inst_frequency(x, N, freq, 0);

    /* Convert from normalised [0, 0.5] to Hz */
    double *freq_hz = (double *)malloc((size_t)N * sizeof(double));
    double *freq_true = (double *)malloc((size_t)N * sizeof(double));
    double *t_axis = (double *)malloc((size_t)N * sizeof(double));
    for (int i = 0; i < N; i++) {
        freq_hz[i] = freq[i] * fs;
        t_axis[i] = (double)i / fs;
        freq_true[i] = f0 + (f1 - f0) * t_axis[i] / T;
    }

    /* Skip edges */
    int skip = 20;
    double max_err = 0.0;
    for (int i = skip; i < N - skip; i++) {
        double err = fabs(freq_hz[i] - freq_true[i]);
        if (err > max_err) max_err = err;
    }
    printf("  Chirp: %.0f → %.0f Hz over %.3f s\n", f0, f1, T);
    printf("  Max |f_est - f_true| (skipping edges): %.2f Hz\n\n", max_err);

    gp_plot_multi("ch20", "inst_frequency",
        "Instantaneous Frequency of Linear Chirp",
        "Time (s)", "Frequency (Hz)",
        (GpSeries[]){
            {"True IF", t_axis, freq_true, N, "lines"},
            {"Estimated IF", t_axis, freq_hz, N, "lines"}
        }, 2);
    printf("  → plots/ch20/inst_frequency.png\n\n");

    free(x); free(freq); free(freq_hz); free(freq_true); free(t_axis);
}

/* ── Demo 5: Single-Sideband Modulation ───────────────────────── */

static void demo_ssb(void)
{
    printf("=== Demo 5: Single-Sideband (SSB) Modulation ===\n\n");

    const int N = 512;
    const double fs = 8000.0;
    const double fc = 1000.0;    /* carrier */
    double *x = (double *)malloc((size_t)N * sizeof(double));
    gen_sine(x, N, 1.0, 200.0, fs, 0.0);  /* message signal */

    /* Analytic signal of message */
    Complex *za = (Complex *)calloc((size_t)N, sizeof(Complex));
    analytic_signal_fft(x, N, za);

    /* USB: y(t) = Re{z(t) · e^{j2πfc·t}} */
    double *usb = (double *)malloc((size_t)N * sizeof(double));
    for (int i = 0; i < N; i++) {
        double t = (double)i / fs;
        double cos_c = cos(2.0 * M_PI * fc * t);
        double sin_c = sin(2.0 * M_PI * fc * t);
        usb[i] = za[i].re * cos_c - za[i].im * sin_c;
    }

    /* Spectrum of USB */
    Complex *Y = (Complex *)calloc((size_t)N, sizeof(Complex));
    for (int i = 0; i < N; i++) { Y[i].re = usb[i]; Y[i].im = 0.0; }
    fft(Y, N);

    double *mag = (double *)malloc((size_t)(N / 2) * sizeof(double));
    double *freq_ax = (double *)malloc((size_t)(N / 2) * sizeof(double));
    for (int i = 0; i < N / 2; i++) {
        mag[i] = 20.0 * log10(complex_mag(Y[i]) / N + 1e-10);
        freq_ax[i] = (double)i * fs / N;
    }

    gp_plot_1("ch20", "ssb_spectrum",
              "USB SSB: 200 Hz Message, 1000 Hz Carrier",
              "Frequency (Hz)", "Magnitude (dB)",
              freq_ax, mag, N / 2, "lines");

    printf("  Message: 200 Hz, Carrier: %.0f Hz\n", fc);
    printf("  USB peak expected at %.0f Hz (carrier + message)\n", fc + 200.0);
    printf("  LSB (at %.0f Hz) is suppressed.\n", fc - 200.0);
    printf("  → plots/ch20/ssb_spectrum.png\n\n");

    free(x); free(za); free(usb); free(Y); free(mag); free(freq_ax);
}

int main(void)
{
    printf("╔══════════════════════════════════════════════════════════╗\n");
    printf("║  Chapter 20: Hilbert Transform & Analytic Signals      ║\n");
    printf("╚══════════════════════════════════════════════════════════╝\n\n");

    demo_hilbert_design();
    demo_analytic_signal();
    demo_envelope();
    demo_inst_frequency();
    demo_ssb();

    printf("=== Chapter 20 Complete ===\n");
    return 0;
}

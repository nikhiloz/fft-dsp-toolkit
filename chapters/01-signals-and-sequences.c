/**
 * Chapter 1 Demo: Discrete-Time Signals & Sequences
 *
 * Generates and displays the fundamental signal types:
 *   1. Unit impulse (Kronecker delta)
 *   2. Unit step
 *   3. Real exponential (decaying + growing)
 *   4. Sinusoid (cosine wave)
 *   5. Complex exponential (unit circle trace)
 *   6. Multi-tone composite
 *   7. Linear chirp
 *   8. White noise + Gaussian noise
 *
 * Build:  make chapters && ./build/bin/ch01s
 *
 * References:
 *   - Oppenheim & Willsky, Signals and Systems, Ch 1-2
 *   - Oppenheim & Schafer, Discrete-Time Signal Processing, Ch 2
 *   - Lyons, Understanding DSP, Ch 1
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include "signal_gen.h"
#include "dsp_utils.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Print first `count` samples of a signal, with an ASCII bar chart */
static void print_signal(const char *name, const double *x, int n, int count)
{
    if (count > n) count = n;
    printf("\n── %s (%d samples shown) ──\n", name, count);
    for (int i = 0; i < count; i++) {
        printf("  x[%2d] = %+8.4f  ", i, x[i]);
        /* Simple ASCII bar */
        int bar = (int)(x[i] * 20.0);
        if (bar >= 0) {
            printf("  |");
            for (int b = 0; b < bar && b < 40; b++) printf("█");
        } else {
            int spaces = 20 + bar; /* bar is negative */
            if (spaces < 0) spaces = 0;
            for (int b = 0; b < spaces; b++) printf(" ");
            for (int b = 0; b < -bar && b < 40; b++) printf("█");
            printf("|");
        }
        printf("\n");
    }
}

static void print_complex_signal(const char *name, const Complex *x, int n, int count)
{
    if (count > n) count = n;
    printf("\n── %s (%d samples shown) ──\n", name, count);
    for (int i = 0; i < count; i++) {
        printf("  x[%2d] = %+7.4f %+7.4fi   |x| = %.4f\n",
               i, x[i].re, x[i].im, complex_mag(x[i]));
    }
}

int main(void)
{
    printf("=== Chapter 1: Discrete-Time Signals & Sequences ===\n");

    /* ── Demo 1: Unit impulse ────────────────────────────────── */
    {
        double sig[16];
        gen_impulse(sig, 16, 0);
        print_signal("Unit impulse δ[n]", sig, 16, 8);
        printf("  → Single spike at n=0. FFT of this = flat spectrum (all freqs equal).\n");

        gen_impulse(sig, 16, 3);
        print_signal("Shifted impulse δ[n-3]", sig, 16, 8);
        printf("  → Same spike, delayed by 3 samples.\n");
    }

    /* ── Demo 2: Unit step ───────────────────────────────────── */
    {
        double sig[16];
        gen_step(sig, 16, 0);
        print_signal("Unit step u[n]", sig, 16, 8);
        printf("  → Turns on at n=0 and stays on.\n");
        printf("  → First difference: u[n] - u[n-1] = δ[n] (the impulse).\n");
    }

    /* ── Demo 3: Exponentials ────────────────────────────────── */
    {
        double sig[16];

        gen_exponential(sig, 16, 1.0, 0.85);
        print_signal("Decaying exponential (0.85)^n", sig, 16, 12);
        printf("  → |base| < 1: signal decays. Models damped systems.\n");

        gen_exponential(sig, 16, 0.01, 1.15);
        print_signal("Growing exponential (1.15)^n", sig, 16, 12);
        printf("  → |base| > 1: signal grows. Unstable — filter design avoids this!\n");

        gen_exponential(sig, 16, 1.0, -0.9);
        print_signal("Alternating exponential (-0.9)^n", sig, 16, 12);
        printf("  → Negative base: sign alternates each sample. At base=-1: Nyquist frequency.\n");
    }

    /* ── Demo 4: Sinusoid ────────────────────────────────────── */
    {
        double sig[32];
        double sr = 1000.0; /* sample rate */
        double freq = 100.0; /* 100 Hz */

        gen_cosine(sig, 32, 1.0, freq, sr, 0.0);
        print_signal("Cosine: 100 Hz, fs=1000 Hz", sig, 32, 20);
        printf("  → Period = fs/f = %.0f samples. One full cycle every 10 samples.\n",
               sr / freq);
        printf("  → Normalised frequency ω = 2π × %.0f/%.0f = %.4f rad/sample.\n",
               freq, sr, 2.0 * M_PI * freq / sr);
    }

    /* ── Demo 5: Complex exponential ─────────────────────────── */
    {
        Complex sig[16];
        gen_complex_exp(sig, 16, 1.0, 125.0, 1000.0, 0.0);
        print_complex_signal("Complex exponential: 125 Hz, fs=1000", sig, 16, 8);
        printf("  → All magnitudes = 1.0: traces a unit circle in the complex plane.\n");
        printf("  → Real part = cos(ωn), Imaginary part = sin(ωn).\n");
    }

    /* ── Demo 6: Multi-tone composite ────────────────────────── */
    {
        double sig[64];
        double freqs[] = {100.0, 250.0, 400.0};
        double amps[]  = {1.0,   0.5,   0.3};
        gen_multi_tone(sig, 64, freqs, amps, 3, 1000.0);
        print_signal("Multi-tone: 100+250+400 Hz", sig, 64, 20);
        printf("  → Superposition of 3 sinusoids. FFT reveals the 3 peaks.\n");

        /* Show RMS for reference */
        double r = rms(sig, 64);
        printf("  → RMS = %.4f (combine via N-tone formula: sqrt(sum(Ak²)/2)).\n", r);
    }

    /* ── Demo 7: Linear chirp ────────────────────────────────── */
    {
        double sig[64];
        gen_chirp(sig, 64, 1.0, 50.0, 450.0, 1000.0);
        print_signal("Chirp: 50→450 Hz sweep", sig, 64, 32);
        printf("  → Frequency increases with time. Used in radar & room measurement.\n");
    }

    /* ── Demo 8: Noise ───────────────────────────────────────── */
    {
        double white[32], gauss[32];
        gen_white_noise(white, 32, 0.5, 42);
        gen_gaussian_noise(gauss, 32, 0.0, 0.3, 42);

        print_signal("White noise (uniform, amp=0.5)", white, 32, 16);
        printf("  → Uniform distribution in [-0.5, +0.5]. Flat spectrum.\n");

        print_signal("Gaussian noise (μ=0, σ=0.3)", gauss, 32, 16);
        printf("  → Bell-curve distribution. Most natural noise model.\n");

        /* Show statistics */
        printf("\n  White noise:    RMS = %.4f\n", rms(white, 32));
        printf("  Gaussian noise: RMS = %.4f\n", rms(gauss, 32));
    }

    /* ── Demo 9: Signal operations ───────────────────────────── */
    {
        double tone[32], noise[32];
        gen_sine(tone, 32, 1.0, 100.0, 1000.0, 0.0);
        gen_white_noise(noise, 32, 0.3, 99);

        printf("\n── Signal operations ──\n");
        printf("  Clean signal RMS:   %.4f\n", rms(tone, 32));

        signal_add(tone, noise, 32);
        printf("  After adding noise: %.4f  (signal + noise)\n", rms(tone, 32));

        signal_scale(tone, 32, 0.5);
        printf("  After scaling ×0.5: %.4f\n", rms(tone, 32));
        printf("  → Linear operations: the building blocks of all DSP.\n");
    }

    printf("\n=== End of Chapter 1 Demo ===\n");
    return 0;
}

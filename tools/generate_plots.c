/**
 * @file generate_plots.c
 * @brief Generate all gnuplot visualisations for the DSP tutorial suite.
 *
 * Produces PNG plots in plots/chXX/ for every chapter.  These are
 * referenced from the chapter .md tutorial files.
 *
 * Build:  make release    (builds generate_plots alongside other targets)
 * Run:    make plots      (generates all PNGs in plots/)
 *   or:   ./build/bin/generate_plots
 *
 * Requires: gnuplot >= 5.0 (apt install gnuplot)
 *
 * Plot inventory (~30 plots across 14 chapters):
 *
 *   Ch01  signals       : impulse, exponentials, cosine, chirp, multitone
 *   Ch02  sampling      : aliasing, quantization, sinc reconstruction
 *   Ch03  complex       : twiddle factors on unit circle
 *   Ch04  LTI           : convolution smoothing, cross-correlation
 *   Ch05  z-transform   : lowpass response, resonance sharpness
 *   Ch06  freq response : FIR vs IIR, pole Q factor, group delay
 *   Ch07  DFT           : spectrum, zero-padding, standard signals
 *   Ch08  FFT           : two-tone magnitude spectrum
 *   Ch09  windows       : window shapes, spectral leakage comparison
 *   Ch10  FIR filters   : sinc kernel, noise reduction
 *   Ch11  IIR design    : Butterworth orders, Chebyshev ripple, filtering
 *   Ch12  structures    : DF1 vs DF2T, coefficient sensitivity
 *   Ch13  spectral      : rect vs Hann windowed spectrum
 *   Ch14  PSD/Welch     : periodogram vs Welch, resolution trade-off
 *   Ch15  correlation   : autocorr pitch, white noise autocorr
 *   Ch16  streaming     : OLA vs direct convolution
 *   Ch18  fixed-point   : Q15 SQNR, float vs Q15 FIR
 *   Ch19  advanced FFT  : Goertzel spectrum, sliding DFT
 *   Ch17  multirate     : decimation, interpolation, polyphase comparison
 *   Ch20  hilbert       : analytic signal envelope, instantaneous frequency
 *   Ch21  averaging     : coherent averaging SNR, EMA vs MA, median filter
 *   Ch22  advanced FIR  : Remez LP vs window, Remez bandpass
 *   Ch23  adaptive      : LMS learning curve, LMS vs NLMS convergence
 *   Ch24  linear pred   : AR spectrum (order 4/10/20), LPC round-trip
 *   Ch25  param spectral: MUSIC super-resolution, Capon vs FFT
 *   Ch26  cepstrum/MFCC : real cepstrum, liftering, Mel filterbank
 *   Ch27  2D DSP        : Gaussian blur, Sobel edges
 *   Ch30  capstone      : full pipeline time + frequency domain
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gnuplot.h"
#include "signal_gen.h"
#include "dsp_utils.h"
#include "fft.h"
#include "filter.h"
#include "iir.h"
#include "convolution.h"
#include "spectrum.h"
#include "correlation.h"
#include "fixed_point.h"
#include "advanced_fft.h"
#include "streaming.h"
#include "multirate.h"
#include "hilbert.h"
#include "averaging.h"
#include "remez.h"
#include "adaptive.h"
#include "lpc.h"
#include "spectral_est.h"
#include "cepstrum.h"
#include "dsp2d.h"
#include "realtime.h"
#include "optimization.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ================================================================== */
/*  Helper: evaluate H(z) at z = e^{jω} for a simple system          */
/* ================================================================== */

/*
 * Evaluate FIR polynomial B(z) = b[0] + b[1]*z^{-1} + ...
 * at z = e^{jω} using Horner's method.
 */
static Complex eval_poly(const double *b, int len, Complex z)
{
    double mag_sq = z.re * z.re + z.im * z.im;
    Complex z_inv;
    if (mag_sq < 1e-30) { z_inv.re = 0; z_inv.im = 0; }
    else { z_inv.re = z.re / mag_sq; z_inv.im = -z.im / mag_sq; }

    Complex r = { b[len - 1], 0.0 };
    for (int i = len - 2; i >= 0; i--) {
        r = complex_mul(r, z_inv);
        r.re += b[i];
    }
    return r;
}

/*
 * Evaluate IIR transfer function H(z) = B(z) / A(z)
 * where A(z) = 1 + a[0]*z^{-1} + a[1]*z^{-2} + ...
 */
static Complex eval_hz(const double *b, int blen,
                       const double *a, int alen, Complex z)
{
    Complex num = eval_poly(b, blen, z);
    if (alen <= 0) return num;

    double af[32];
    af[0] = 1.0;
    for (int i = 0; i < alen && i < 31; i++) af[i + 1] = a[i];
    Complex den = eval_poly(af, alen + 1, z);

    double d2 = den.re * den.re + den.im * den.im;
    if (d2 < 1e-30) { Complex inf = {1e10, 0}; return inf; }
    Complex h;
    h.re = (num.re * den.re + num.im * den.im) / d2;
    h.im = (num.im * den.re - num.re * den.im) / d2;
    return h;
}

static double quantize_sample(double x, int bits)
{
    double levels = (double)(1 << bits);
    double half = levels / 2.0;
    double q = floor(x * half + 0.5) / half;
    if (q > 1.0) q = 1.0;
    if (q < -1.0) q = -1.0;
    return q;
}

static double sinc_val(double x)
{
    if (fabs(x) < 1e-12) return 1.0;
    return sin(M_PI * x) / (M_PI * x);
}

/* ================================================================== */
/*  Chapter 1: Discrete-Time Signals & Sequences                      */
/* ================================================================== */

static void plot_ch01(void)
{
    printf("  Ch01: signals ...\n");
    gp_init("ch01");

    /* 1. Unit impulse */
    {
        double sig[32];
        gen_impulse(sig, 32, 0);
        gp_plot_1("ch01", "impulse", "Unit Impulse {/Symbol d}[n]",
                  "Sample n", "Amplitude", NULL, sig, 32, "impulses");
    }

    /* 2. Exponentials overlay: decaying, growing, alternating */
    {
        double decay[20], grow[20], alt[20];
        gen_exponential(decay, 20, 1.0, 0.85);
        gen_exponential(grow,  20, 0.01, 1.15);
        gen_exponential(alt,   20, 1.0, -0.9);
        GpSeries s[] = {
            { "Decaying (0.85)^n", NULL, decay, 20, "linespoints" },
            { "Growing (1.15)^n",  NULL, grow,  20, "linespoints" },
            { "Alternating (-0.9)^n", NULL, alt, 20, "linespoints" },
        };
        gp_plot_multi("ch01", "exponentials",
                      "Real Exponential Signals",
                      "Sample n", "x[n]", s, 3);
    }

    /* 3. Cosine wave */
    {
        double sig[40];
        gen_cosine(sig, 40, 1.0, 100.0, 1000.0, 0.0);
        gp_plot_1("ch01", "cosine", "Cosine: 100 Hz at f_s = 1000 Hz",
                  "Sample n", "x[n]", NULL, sig, 40, "linespoints");
    }

    /* 4. Chirp */
    {
        double sig[128];
        gen_chirp(sig, 128, 1.0, 50.0, 450.0, 1000.0);
        gp_plot_1("ch01", "chirp", "Linear Chirp: 50 {/Symbol \\256} 450 Hz",
                  "Sample n", "x[n]", NULL, sig, 128, "lines");
    }

    /* 5. Multi-tone */
    {
        double sig[128];
        double freqs[] = { 100.0, 250.0, 400.0 };
        double amps[]  = { 1.0,   0.5,   0.3  };
        gen_multi_tone(sig, 128, freqs, amps, 3, 1000.0);
        gp_plot_1("ch01", "multitone",
                  "Multi-tone: 100 + 250 + 400 Hz",
                  "Sample n", "x[n]", NULL, sig, 128, "lines");
    }
}

/* ================================================================== */
/*  Chapter 2: Sampling, Aliasing & Nyquist                           */
/* ================================================================== */

static void plot_ch02(void)
{
    printf("  Ch02: sampling ...\n");
    gp_init("ch02");

    /* 1. Aliasing: 300 Hz vs 700 Hz vs 1300 Hz at fs=1000 */
    {
        double s1[64], s2[64], s3[64];
        gen_sine(s1, 64, 1.0, 300.0, 1000.0, 0.0);
        gen_sine(s2, 64, 1.0, 700.0, 1000.0, 0.0);
        gen_sine(s3, 64, 1.0, 1300.0, 1000.0, 0.0);
        GpSeries s[] = {
            { "300 Hz (original)",   NULL, s1, 64, "linespoints" },
            { "700 Hz (alias)",      NULL, s2, 64, "linespoints" },
            { "1300 Hz (alias)",     NULL, s3, 64, "linespoints" },
        };
        gp_plot_multi("ch02", "aliasing",
                      "Aliasing: Three Frequencies, Same Samples (f_s = 1000 Hz)",
                      "Sample n", "x[n]", s, 3);
    }

    /* 2. Quantization: original vs 4-bit */
    {
        double orig[64], quant[64], noise[64];
        gen_sine(orig, 64, 0.9, 50.0, 1000.0, 0.0);
        for (int i = 0; i < 64; i++) {
            quant[i] = quantize_sample(orig[i], 4);
            noise[i] = orig[i] - quant[i];
        }
        GpSeries s[] = {
            { "Original",         NULL, orig,  64, "lines" },
            { "4-bit Quantized",  NULL, quant, 64, "linespoints" },
            { "Quant. Error",     NULL, noise, 64, "impulses" },
        };
        gp_plot_multi("ch02", "quantization",
                      "4-bit Quantization (16 levels)",
                      "Sample n", "Amplitude", s, 3);
    }

    /* 3. Sinc reconstruction */
    {
        int ns = 16, factor = 8;
        double fs = 200.0, freq = 50.0;
        double samples[16];
        gen_sine(samples, ns, 1.0, freq, fs, 0.0);

        int nr = ns * factor;
        double recon[128], true_sig[128], x_recon[128], x_samp[16];
        for (int k = 0; k < nr; k++) {
            double t = (double)k / factor;
            double val = 0.0;
            for (int n = 0; n < ns; n++) val += samples[n] * sinc_val(t - n);
            recon[k] = val;
            true_sig[k] = sin(2.0 * M_PI * freq * (double)k / (fs * factor));
            x_recon[k] = (double)k / factor;
        }
        for (int i = 0; i < ns; i++) x_samp[i] = (double)i;

        GpSeries s[] = {
            { "Discrete Samples", x_samp, samples,  ns, "points" },
            { "Sinc Reconstruction", x_recon, recon,  nr, "lines" },
            { "True Signal",     x_recon, true_sig, nr, "lines" },
        };
        gp_plot_multi("ch02", "reconstruction",
                      "Sinc Interpolation Reconstruction",
                      "Sample Index", "Amplitude", s, 3);
    }
}

/* ================================================================== */
/*  Chapter 3: Complex Numbers                                        */
/* ================================================================== */

static void plot_ch03(void)
{
    printf("  Ch03: complex numbers ...\n");
    gp_init("ch03");

    /* Twiddle factors W_8^k on unit circle */
    {
        double xr[8], xi[8];
        for (int k = 0; k < 8; k++) {
            double angle = -2.0 * M_PI * k / 8.0;
            xr[k] = cos(angle);
            xi[k] = sin(angle);
        }

        FILE *gp = gp_open("ch03", "twiddle_factors", 600, 600);
        if (gp) {
            fprintf(gp, "set title 'Twiddle Factors W_8^k on the Unit Circle'\n");
            fprintf(gp, "set xlabel 'Real'\n");
            fprintf(gp, "set ylabel 'Imaginary'\n");
            fprintf(gp, "set size ratio 1\n");
            fprintf(gp, "set xrange [-1.4:1.4]\n");
            fprintf(gp, "set yrange [-1.4:1.4]\n");
            /* Draw unit circle */
            fprintf(gp, "set parametric\n");
            fprintf(gp, "set trange [0:2*pi]\n");
            fprintf(gp, "plot cos(t), sin(t) w lines lw 1 lc rgb '#CCCCCC' "
                        "title 'Unit Circle', "
                        "'-' w points pt 7 ps 2 lc rgb '#2166AC' "
                        "title 'W_8^k'\n");
            gp_send_xy(gp, xr, xi, 8);
            gp_close(gp);
        }
    }
}

/* ================================================================== */
/*  Chapter 4: LTI Systems & Convolution                              */
/* ================================================================== */

static void plot_ch04(void)
{
    printf("  Ch04: LTI systems ...\n");
    gp_init("ch04");

    /* 1. Moving average convolution: input rect → smoothed output */
    {
        double x[] = {0, 0, 1, 1, 1, 1, 0, 0};
        double h[] = {1.0/3, 1.0/3, 1.0/3};
        double y[10];
        convolve(x, 8, h, 3, y);

        GpSeries s[] = {
            { "Input (rect pulse)",  NULL, x, 8,  "impulses" },
            { "3-pt Moving Average", NULL, y, 10, "linespoints" },
        };
        gp_plot_multi("ch04", "convolution",
                      "Convolution: Rectangular Pulse * Moving Average",
                      "Sample n", "Amplitude", s, 2);
    }

    /* 2. Cross-correlation delay estimation */
    {
        double x[32], y[32];
        memset(x, 0, sizeof(x));
        memset(y, 0, sizeof(y));
        x[4] = 1.0; x[5] = 0.8; x[6] = 0.5; x[7] = 0.2;
        y[14] = 1.0; y[15] = 0.8; y[16] = 0.5; y[17] = 0.2;

        double r[63];
        int rlen = cross_correlate(x, 32, y, 32, r);

        /* Extract lags -31..+31 mapped to r[0..62] */
        double lags[63], rval[63];
        for (int i = 0; i < rlen; i++) {
            lags[i] = (double)(i - 31);
            rval[i] = r[i];
        }
        GpSeries s[] = {
            { "Cross-Correlation", lags, rval, rlen, "lines" },
        };
        gp_plot_multi("ch04", "cross_correlation",
                      "Cross-Correlation: Delay = 10 Samples",
                      "Lag (samples)", "R_{xy}[l]", s, 1);
    }
}

/* ================================================================== */
/*  Chapter 5: Z-Transform                                            */
/* ================================================================== */

static void plot_ch05(void)
{
    printf("  Ch05: z-transform ...\n");
    gp_init("ch05");

    /* 1. Frequency response of 2-point average lowpass */
    {
        double b[] = {0.5, 0.5};
        int np = 200;
        double freq[200], mag[200];
        for (int i = 0; i < np; i++) {
            double omega = M_PI * (double)i / (double)(np - 1);
            freq[i] = omega / M_PI;
            Complex z = complex_from_polar(1.0, omega);
            Complex Hz = eval_poly(b, 2, z);
            double m = complex_mag(Hz);
            mag[i] = (m > 1e-10) ? 20.0 * log10(m) : -100.0;
        }
        gp_plot_spectrum("ch05", "lowpass_response",
                         "Frequency Response: (1 + z^{-1})/2  (2-Point Average)",
                         freq, mag, np);
    }

    /* 2. Resonance: pole radius {0.5, 0.8, 0.95, 0.99} at θ=π/4 */
    {
        double theta = M_PI / 4.0;
        double radii[] = {0.5, 0.8, 0.95, 0.99};
        char *labels[] = {"r = 0.50", "r = 0.80", "r = 0.95", "r = 0.99"};
        int np = 200;
        double freq[200];
        double mag0[200], mag1[200], mag2[200], mag3[200];
        double *mags[] = { mag0, mag1, mag2, mag3 };

        for (int i = 0; i < np; i++) {
            double omega = M_PI * (double)i / (double)(np - 1);
            freq[i] = omega / M_PI;
            for (int ri = 0; ri < 4; ri++) {
                double r = radii[ri];
                double a[] = { -2.0 * r * cos(theta), r * r };
                double b[] = { 1.0 };
                Complex z = complex_from_polar(1.0, omega);
                Complex Hz = eval_hz(b, 1, a, 2, z);
                double m = complex_mag(Hz);
                mags[ri][i] = (m > 1e-10) ? 20.0 * log10(m) : -100.0;
                if (mags[ri][i] < -40.0) mags[ri][i] = -40.0;
                if (mags[ri][i] > 40.0) mags[ri][i] = 40.0;
            }
        }
        GpSeries s[] = {
            { labels[0], freq, mag0, np, "lines" },
            { labels[1], freq, mag1, np, "lines" },
            { labels[2], freq, mag2, np, "lines" },
            { labels[3], freq, mag3, np, "lines" },
        };
        gp_plot_multi("ch05", "resonance",
                      "Resonance Sharpness vs Pole Radius ({/Symbol q} = {/Symbol p}/4)",
                      "Normalised Frequency ({/Symbol w}/{/Symbol p})",
                      "Magnitude (dB)", s, 4);
    }
}

/* ================================================================== */
/*  Chapter 6: Frequency Response                                     */
/* ================================================================== */

static void plot_ch06(void)
{
    printf("  Ch06: frequency response ...\n");
    gp_init("ch06");

    /* 1. FIR vs IIR magnitude + phase */
    {
        double b_fir[] = {0.2, 0.2, 0.2, 0.2, 0.2};
        double b_iir[] = {0.1};
        double a_iir[] = {-0.9};
        int np = 200;
        double freq[200], fir_mag[200], iir_mag[200];
        double fir_phase[200], iir_phase[200];

        for (int i = 0; i < np; i++) {
            double omega = M_PI * (double)i / (double)(np - 1);
            freq[i] = omega / M_PI;
            Complex z = complex_from_polar(1.0, omega);

            Complex Hfir = eval_poly(b_fir, 5, z);
            Complex Hiir = eval_hz(b_iir, 1, a_iir, 1, z);

            double mf = complex_mag(Hfir);
            double mi = complex_mag(Hiir);
            fir_mag[i] = (mf > 1e-10) ? 20.0 * log10(mf) : -60.0;
            iir_mag[i] = (mi > 1e-10) ? 20.0 * log10(mi) : -60.0;
            fir_phase[i] = complex_phase(Hfir) * 180.0 / M_PI;
            iir_phase[i] = complex_phase(Hiir) * 180.0 / M_PI;
        }

        /* Magnitude */
        GpSeries sm[] = {
            { "FIR 5-tap",    freq, fir_mag, np, "lines" },
            { "IIR 1st-order", freq, iir_mag, np, "lines" },
        };
        gp_plot_multi("ch06", "fir_vs_iir_magnitude",
                      "FIR vs IIR: Magnitude Response",
                      "Normalised Frequency ({/Symbol w}/{/Symbol p})",
                      "Magnitude (dB)", sm, 2);

        /* Phase */
        GpSeries sp[] = {
            { "FIR 5-tap",    freq, fir_phase, np, "lines" },
            { "IIR 1st-order", freq, iir_phase, np, "lines" },
        };
        gp_plot_multi("ch06", "fir_vs_iir_phase",
                      "FIR vs IIR: Phase Response",
                      "Normalised Frequency ({/Symbol w}/{/Symbol p})",
                      "Phase (degrees)", sp, 2);
    }

    /* 2. Pole radius Q factor */
    {
        double theta = M_PI / 4.0;
        double radii[] = {0.5, 0.8, 0.95, 0.99};
        char *labels[] = {"r = 0.50", "r = 0.80", "r = 0.95", "r = 0.99"};
        int np = 300;
        double freq[300], m0[300], m1[300], m2[300], m3[300];
        double *arrs[] = { m0, m1, m2, m3 };

        for (int i = 0; i < np; i++) {
            double omega = M_PI * (double)i / (double)(np - 1);
            freq[i] = omega / M_PI;
            for (int ri = 0; ri < 4; ri++) {
                double r = radii[ri];
                double a[] = { -2.0 * r * cos(theta), r * r };
                double b[] = { 1.0 };
                Complex z = complex_from_polar(1.0, omega);
                Complex Hz = eval_hz(b, 1, a, 2, z);
                double m = complex_mag(Hz);
                arrs[ri][i] = (m > 1e-10) ? 20.0 * log10(m) : -60.0;
                if (arrs[ri][i] < -40.0) arrs[ri][i] = -40.0;
                if (arrs[ri][i] > 50.0) arrs[ri][i] = 50.0;
            }
        }
        GpSeries s[] = {
            { labels[0], freq, m0, np, "lines" },
            { labels[1], freq, m1, np, "lines" },
            { labels[2], freq, m2, np, "lines" },
            { labels[3], freq, m3, np, "lines" },
        };
        gp_plot_multi("ch06", "pole_radius_q",
                      "Resonance Q Factor: Effect of Pole Radius",
                      "Normalised Frequency ({/Symbol w}/{/Symbol p})",
                      "Magnitude (dB)", s, 4);
    }

    /* 3. Group delay: FIR constant vs IIR varying */
    {
        double b_fir[] = {0.1, 0.2, 0.4, 0.2, 0.1};
        int np = 150;
        double freq[150], gd_fir[150], gd_iir[150];
        double dw = 0.001;

        for (int i = 0; i < np; i++) {
            double omega = M_PI * (double)i / (double)(np - 1);
            freq[i] = omega / M_PI;

            /* FIR group delay: d(phase)/d(omega) via central diff */
            Complex z1 = complex_from_polar(1.0, omega - dw);
            Complex z2 = complex_from_polar(1.0, omega + dw);
            Complex H1 = eval_poly(b_fir, 5, z1);
            Complex H2 = eval_poly(b_fir, 5, z2);
            double p1 = complex_phase(H1), p2 = complex_phase(H2);
            /* Unwrap */
            double dp = p2 - p1;
            while (dp > M_PI) dp -= 2.0 * M_PI;
            while (dp < -M_PI) dp += 2.0 * M_PI;
            gd_fir[i] = -dp / (2.0 * dw);

            /* IIR group delay */
            double a_iir[] = {-0.9};
            double b_iir[] = {0.1};
            H1 = eval_hz(b_iir, 1, a_iir, 1, z1);
            H2 = eval_hz(b_iir, 1, a_iir, 1, z2);
            p1 = complex_phase(H1); p2 = complex_phase(H2);
            dp = p2 - p1;
            while (dp > M_PI) dp -= 2.0 * M_PI;
            while (dp < -M_PI) dp += 2.0 * M_PI;
            gd_iir[i] = -dp / (2.0 * dw);
        }
        GpSeries s[] = {
            { "FIR (symmetric, constant)", freq, gd_fir, np, "lines" },
            { "IIR (1st-order, varying)",  freq, gd_iir, np, "lines" },
        };
        gp_plot_multi("ch06", "group_delay",
                      "Group Delay: FIR vs IIR",
                      "Normalised Frequency ({/Symbol w}/{/Symbol p})",
                      "Group Delay (samples)", s, 2);
    }
}

/* ================================================================== */
/*  Chapter 7: DFT Theory                                            */
/* ================================================================== */

static void plot_ch07(void)
{
    printf("  Ch07: DFT theory ...\n");
    gp_init("ch07");

    /* 1. DFT spectrum of 250 Hz sine at fs=2000, N=16 */
    {
        int N = 16;
        double sig[16];
        gen_sine(sig, N, 1.0, 250.0, 2000.0, 0.0);

        Complex X[16];
        fft_real(sig, X, N);

        double bins[9], mag[9]; /* 0..N/2 */
        for (int k = 0; k <= N / 2; k++) {
            bins[k] = (double)k;
            mag[k]  = complex_mag(X[k]) / (N / 2.0);
        }
        gp_plot_1("ch07", "dft_spectrum",
                  "16-point DFT of 250 Hz Sine (f_s = 2000 Hz)",
                  "Frequency Bin k", "|X[k]| (normalised)",
                  bins, mag, N / 2 + 1, "impulses");
    }

    /* 2. Zero-padding: N=8 vs N=32 */
    {
        double sig8[8];
        gen_sine(sig8, 8, 1.0, 1.0, 8.0, 0.0);  /* 1 Hz at fs=8 */

        /* N = 8 */
        Complex X8[8];
        fft_real(sig8, X8, 8);
        double b8[5], m8[5];
        for (int k = 0; k <= 4; k++) {
            b8[k] = (double)k / 8.0;
            m8[k] = complex_mag(X8[k]);
        }

        /* N = 32 (zero-padded) */
        double sig32[32];
        memset(sig32, 0, sizeof(sig32));
        memcpy(sig32, sig8, 8 * sizeof(double));
        Complex X32[32];
        fft_real(sig32, X32, 32);
        double b32[17], m32[17];
        for (int k = 0; k <= 16; k++) {
            b32[k] = (double)k / 32.0;
            m32[k] = complex_mag(X32[k]);
        }

        GpSeries s[] = {
            { "N = 8 (coarse)",     b8,  m8,  5,  "linespoints" },
            { "N = 32 (zero-padded)", b32, m32, 17, "linespoints" },
        };
        gp_plot_multi("ch07", "zero_padding",
                      "Zero-Padding: Spectral Interpolation",
                      "Normalised Frequency (f/f_s)",
                      "|X[k]|", s, 2);
    }

    /* 3. Standard signals DFT — 4-panel (impulse, DC, alternating, cosine) */
    {
        int N = 16;
        double imp[16], dc[16], alt[16], cos_sig[16];
        gen_impulse(imp, N, 0);
        for (int i = 0; i < N; i++) { dc[i] = 1.0; alt[i] = (i % 2 == 0) ? 1.0 : -1.0; }
        gen_cosine(cos_sig, N, 1.0, 2.0, 16.0, 0.0);

        Complex Xi[16], Xd[16], Xa[16], Xc[16];
        fft_real(imp, Xi, N);
        fft_real(dc, Xd, N);
        fft_real(alt, Xa, N);
        fft_real(cos_sig, Xc, N);

        FILE *gp = gp_open("ch07", "standard_signals_dft", 900, 700);
        if (gp) {
            fprintf(gp, "set multiplot layout 2,2 title "
                        "'DFT of Standard Signals (N = 16)'\n");
            const char *names[] = {"Impulse {/Symbol d}[n]", "DC (constant)",
                                   "Alternating (-1)^n", "Cosine (2 Hz)"};
            Complex *spectra[] = {Xi, Xd, Xa, Xc};
            for (int p = 0; p < 4; p++) {
                fprintf(gp, "set title '%s'\n", names[p]);
                fprintf(gp, "set xlabel 'Bin k'\n");
                fprintf(gp, "set ylabel '|X[k]|'\n");
                fprintf(gp, "plot '-' w impulses lw 2 notitle\n");
                for (int k = 0; k <= N / 2; k++)
                    fprintf(gp, "%d\t%.6f\n", k, complex_mag(spectra[p][k]));
                fprintf(gp, "e\n");
            }
            fprintf(gp, "unset multiplot\n");
            gp_close(gp);
        }
    }
}

/* ================================================================== */
/*  Chapter 8: FFT Fundamentals                                       */
/* ================================================================== */

static void plot_ch08(void)
{
    printf("  Ch08: FFT fundamentals ...\n");
    gp_init("ch08");

    /* Two-tone: 440 Hz + 1000 Hz at fs=4000, N=256 */
    {
        int N = 256;
        double sig[256];
        double f1[] = {440.0, 1000.0};
        double a1[] = {1.0, 0.7};
        gen_multi_tone(sig, N, f1, a1, 2, 4000.0);

        Complex X[256];
        fft_real(sig, X, N);

        int nb = N / 2 + 1;
        double freq[129], mag_db[129];
        for (int k = 0; k < nb; k++) {
            freq[k] = (double)k * 4000.0 / N;
            double m = complex_mag(X[k]) / (N / 2.0);
            mag_db[k] = (m > 1e-10) ? 20.0 * log10(m) : -100.0;
        }

        FILE *gp = gp_open("ch08", "fft_two_tones", 800, 500);
        if (gp) {
            fprintf(gp, "set title '256-Point FFT: 440 Hz + 1000 Hz'\n");
            fprintf(gp, "set xlabel 'Frequency (Hz)'\n");
            fprintf(gp, "set ylabel 'Magnitude (dB)'\n");
            fprintf(gp, "set xrange [0:2000]\n");
            fprintf(gp, "set yrange [-60:5]\n");
            fprintf(gp, "plot '-' w lines lw 2 notitle\n");
            gp_send_xy(gp, freq, mag_db, nb);
            gp_close(gp);
        }
    }
}

/* ================================================================== */
/*  Chapter 9: Window Functions                                       */
/* ================================================================== */

static void plot_ch09(void)
{
    printf("  Ch09: windows ...\n");
    gp_init("ch09");

    int N = 64;
    double rect[64], hann[64], hamm[64], black[64];

    /* Generate window shapes */
    for (int i = 0; i < N; i++) {
        rect[i]  = 1.0;
        hann[i]  = 0.5 * (1.0 - cos(2.0 * M_PI * i / (N - 1)));
        hamm[i]  = 0.54 - 0.46 * cos(2.0 * M_PI * i / (N - 1));
        black[i] = 0.42 - 0.5 * cos(2.0 * M_PI * i / (N - 1))
                        + 0.08 * cos(4.0 * M_PI * i / (N - 1));
    }

    /* 1. Window shapes */
    {
        GpSeries s[] = {
            { "Rectangular", NULL, rect,  N, "lines" },
            { "Hann",        NULL, hann,  N, "lines" },
            { "Hamming",     NULL, hamm,  N, "lines" },
            { "Blackman",    NULL, black, N, "lines" },
        };
        gp_plot_multi("ch09", "window_shapes",
                      "Window Functions (N = 64)",
                      "Sample n", "w[n]", s, 4);
    }

    /* 2. Spectral leakage comparison: 440 Hz off-bin at fs=4000, N=256 */
    {
        int SN = 256;
        double sig[256];
        gen_sine(sig, SN, 1.0, 440.0, 4000.0, 0.0);

        double w_rect[256], w_hann[256], w_hamm[256], w_black[256];
        for (int i = 0; i < SN; i++) {
            w_rect[i]  = 1.0;
            w_hann[i]  = 0.5 * (1.0 - cos(2.0 * M_PI * i / (SN - 1)));
            w_hamm[i]  = 0.54 - 0.46 * cos(2.0 * M_PI * i / (SN - 1));
            w_black[i] = 0.42 - 0.5 * cos(2.0 * M_PI * i / (SN - 1))
                              + 0.08 * cos(4.0 * M_PI * i / (SN - 1));
        }

        double *wins[] = { w_rect, w_hann, w_hamm, w_black };
        char *wnames[] = { "Rectangular", "Hann", "Hamming", "Blackman" };
        int nb = SN / 2 + 1;
        double freq[129], m0[129], m1[129], m2[129], m3[129];
        double *marrs[] = { m0, m1, m2, m3 };

        for (int wi = 0; wi < 4; wi++) {
            double windowed[256];
            for (int i = 0; i < SN; i++) windowed[i] = sig[i] * wins[wi][i];
            Complex X[256];
            fft_real(windowed, X, SN);
            for (int k = 0; k < nb; k++) {
                if (wi == 0) freq[k] = (double)k * 4000.0 / SN;
                double m = complex_mag(X[k]) / (SN / 2.0);
                marrs[wi][k] = (m > 1e-10) ? 20.0 * log10(m) : -100.0;
            }
        }

        GpSeries s[] = {
            { wnames[0], freq, m0, nb, "lines" },
            { wnames[1], freq, m1, nb, "lines" },
            { wnames[2], freq, m2, nb, "lines" },
            { wnames[3], freq, m3, nb, "lines" },
        };

        FILE *gp = gp_open("ch09", "spectral_leakage", 900, 500);
        if (gp) {
            fprintf(gp, "set title 'Spectral Leakage: 440 Hz with Different Windows'\n");
            fprintf(gp, "set xlabel 'Frequency (Hz)'\n");
            fprintf(gp, "set ylabel 'Magnitude (dB)'\n");
            fprintf(gp, "set xrange [0:1000]\n");
            fprintf(gp, "set yrange [-100:5]\n");
            fprintf(gp, "plot ");
            for (int i = 0; i < 4; i++) {
                if (i) fprintf(gp, ", ");
                fprintf(gp, "'-' w lines lw 2 title '%s'", s[i].label);
            }
            fprintf(gp, "\n");
            for (int i = 0; i < 4; i++)
                gp_send_xy(gp, s[i].x, s[i].y, s[i].n);
            gp_close(gp);
        }
    }
}

/* ================================================================== */
/*  Chapter 10: FIR Digital Filters                                   */
/* ================================================================== */

static void plot_ch10(void)
{
    printf("  Ch10: FIR filters ...\n");
    gp_init("ch10");

    /* 1. 31-tap lowpass kernel */
    {
        int ntaps = 31;
        double h[31];
        fir_lowpass(h, ntaps, 0.2);
        gp_plot_1("ch10", "lowpass_kernel",
                  "31-Tap Lowpass FIR Kernel (f_c = 0.2)",
                  "Tap n", "h[n]", NULL, h, ntaps, "impulses");
    }

    /* 2. Noise reduction: clean + noise → filter */
    {
        int N = 256;
        double clean[256], noisy[256], noise[256], filtered[256];
        gen_sine(clean, N, 1.0, 200.0, 4000.0, 0.0);
        gen_white_noise(noise, N, 0.5, 42);
        for (int i = 0; i < N; i++) noisy[i] = clean[i] + noise[i];

        int ntaps = 31;
        double h[31];
        fir_lowpass(h, ntaps, 0.2);
        fir_filter(noisy, filtered, N, h, ntaps);

        GpSeries s[] = {
            { "Clean Signal",  NULL, clean,    N, "lines" },
            { "Noisy Signal",  NULL, noisy,    N, "lines" },
            { "Filtered (FIR)", NULL, filtered, N, "lines" },
        };
        gp_plot_multi("ch10", "noise_reduction",
                      "FIR Lowpass Noise Reduction",
                      "Sample n", "Amplitude", s, 3);
    }
}

/* ================================================================== */
/*  Chapter 11: IIR Filter Design                                     */
/* ================================================================== */

static void plot_ch11(void)
{
    printf("  Ch11: IIR design ...\n");
    gp_init("ch11");

    /* 1. Butterworth LP: orders 2, 4, 6 */
    {
        int orders[] = {2, 4, 6};
        char *labels[] = {"Order 2", "Order 4", "Order 6"};
        int np = 300;
        double freq[300], m0[300], m1[300], m2[300];
        double *marrs[] = { m0, m1, m2 };

        for (int oi = 0; oi < 3; oi++) {
            SOSCascade sos;
            butterworth_lowpass(orders[oi], 0.2, &sos);
            sos_freq_response(&sos, marrs[oi], NULL, np);
            /* Convert to dB */
            for (int i = 0; i < np; i++) {
                if (i == 0) freq[i] = 0.0;
                else freq[i] = 0.5 * (double)i / (double)(np - 1);
                marrs[oi][i] = (marrs[oi][i] > 1e-10)
                    ? 20.0 * log10(marrs[oi][i]) : -100.0;
            }
        }
        GpSeries s[] = {
            { labels[0], freq, m0, np, "lines" },
            { labels[1], freq, m1, np, "lines" },
            { labels[2], freq, m2, np, "lines" },
        };
        FILE *gp = gp_open("ch11", "butterworth_orders", 800, 500);
        if (gp) {
            fprintf(gp, "set title 'Butterworth Lowpass: Order Comparison (f_c = 0.2)'\n");
            fprintf(gp, "set xlabel 'Normalised Frequency (f/f_s)'\n");
            fprintf(gp, "set ylabel 'Magnitude (dB)'\n");
            fprintf(gp, "set yrange [-80:5]\n");
            fprintf(gp, "set xrange [0:0.5]\n");
            /* Add -3 dB reference line */
            fprintf(gp, "set arrow from 0,-3 to 0.5,-3 nohead lt 0 lw 1\n");
            fprintf(gp, "set arrow from 0.2,-80 to 0.2,5 nohead lt 0 lw 1\n");
            fprintf(gp, "plot ");
            for (int i = 0; i < 3; i++) {
                if (i) fprintf(gp, ", ");
                fprintf(gp, "'-' w lines lw 2 title '%s'", s[i].label);
            }
            fprintf(gp, "\n");
            for (int i = 0; i < 3; i++)
                gp_send_xy(gp, s[i].x, s[i].y, s[i].n);
            gp_close(gp);
        }
    }

    /* 2. Chebyshev Type I: ripple comparison */
    {
        double ripples[] = {0.5, 1.0, 3.0};
        char *labels[] = {"0.5 dB ripple", "1.0 dB ripple", "3.0 dB ripple"};
        int np = 300;
        double freq[300], m0[300], m1[300], m2[300];
        double *marrs[] = { m0, m1, m2 };

        for (int ri = 0; ri < 3; ri++) {
            SOSCascade sos;
            chebyshev1_lowpass(4, ripples[ri], 0.2, &sos);
            sos_freq_response(&sos, marrs[ri], NULL, np);
            for (int i = 0; i < np; i++) {
                if (ri == 0) freq[i] = 0.5 * (double)i / (double)(np - 1);
                marrs[ri][i] = (marrs[ri][i] > 1e-10)
                    ? 20.0 * log10(marrs[ri][i]) : -100.0;
            }
        }
        GpSeries s[] = {
            { labels[0], freq, m0, np, "lines" },
            { labels[1], freq, m1, np, "lines" },
            { labels[2], freq, m2, np, "lines" },
        };
        FILE *gp = gp_open("ch11", "chebyshev_ripple", 800, 500);
        if (gp) {
            fprintf(gp, "set title 'Chebyshev Type I (Order 4, f_c = 0.2): "
                        "Ripple Comparison'\n");
            fprintf(gp, "set xlabel 'Normalised Frequency (f/f_s)'\n");
            fprintf(gp, "set ylabel 'Magnitude (dB)'\n");
            fprintf(gp, "set yrange [-80:5]\n");
            fprintf(gp, "set xrange [0:0.5]\n");
            fprintf(gp, "plot ");
            for (int i = 0; i < 3; i++) {
                if (i) fprintf(gp, ", ");
                fprintf(gp, "'-' w lines lw 2 title '%s'", s[i].label);
            }
            fprintf(gp, "\n");
            for (int i = 0; i < 3; i++)
                gp_send_xy(gp, s[i].x, s[i].y, s[i].n);
            gp_close(gp);
        }
    }

    /* 3. IIR noise filtering: time domain */
    {
        int N = 256;
        double clean[256], noisy[256], noise[256], filtered[256];
        gen_sine(clean, N, 1.0, 200.0, 4000.0, 0.0);
        gen_white_noise(noise, N, 0.4, 77);
        for (int i = 0; i < N; i++) noisy[i] = clean[i] + noise[i];

        SOSCascade sos;
        butterworth_lowpass(4, 0.15, &sos);
        sos_process_block(&sos, noisy, filtered, N);

        GpSeries s[] = {
            { "Clean",   NULL, clean,    N, "lines" },
            { "Noisy",   NULL, noisy,    N, "lines" },
            { "IIR Filtered", NULL, filtered, N, "lines" },
        };
        gp_plot_multi("ch11", "iir_filtering",
                      "IIR Butterworth Noise Reduction (Order 4)",
                      "Sample n", "Amplitude", s, 3);
    }
}

/* ================================================================== */
/*  Chapter 12: Filter Structures                                     */
/* ================================================================== */

static void plot_ch12(void)
{
    printf("  Ch12: filter structures ...\n");
    gp_init("ch12");

    /* 1. DF1 vs DF2T impulse response */
    {
        SOSCascade sos;
        butterworth_lowpass(4, 0.2, &sos);

        int N = 40;
        /* Use first section */
        Biquad *bq = &sos.sections[0];
        BiquadDF1State  st1 = {0};
        BiquadDF2TState st2 = {0};
        double y1[40], y2[40];

        for (int i = 0; i < N; i++) {
            double x = (i == 0) ? 1.0 : 0.0;
            y1[i] = biquad_process_df1(bq, &st1, x);
            y2[i] = biquad_process_df2t(bq, &st2, x);
        }

        GpSeries s[] = {
            { "Direct Form I",   NULL, y1, N, "linespoints" },
            { "Direct Form II Transposed", NULL, y2, N, "linespoints" },
        };
        gp_plot_multi("ch12", "df1_vs_df2t",
                      "Biquad Impulse Response: DF1 vs DF2T",
                      "Sample n", "y[n]", s, 2);
    }

    /* 2. Coefficient sensitivity */
    {
        SOSCascade orig, pert;
        butterworth_lowpass(4, 0.2, &orig);
        butterworth_lowpass(4, 0.2, &pert);

        /* Perturb a1 coefficients by 0.1% */
        for (int i = 0; i < pert.n_sections; i++) {
            pert.sections[i].a1 *= 1.001;
            pert.sections[i].a2 *= 1.001;
        }

        int np = 300;
        double mag_orig[300], mag_pert[300], freq[300];
        sos_freq_response(&orig, mag_orig, NULL, np);
        sos_freq_response(&pert, mag_pert, NULL, np);

        for (int i = 0; i < np; i++) {
            freq[i] = 0.5 * (double)i / (double)(np - 1);
            mag_orig[i] = (mag_orig[i] > 1e-10)
                ? 20.0 * log10(mag_orig[i]) : -100.0;
            mag_pert[i] = (mag_pert[i] > 1e-10)
                ? 20.0 * log10(mag_pert[i]) : -100.0;
        }

        GpSeries s[] = {
            { "Original Coefficients",    freq, mag_orig, np, "lines" },
            { "Perturbed (+0.1%)",        freq, mag_pert, np, "lines" },
        };
        gp_plot_multi("ch12", "coefficient_sensitivity",
                      "Coefficient Sensitivity: 0.1% Perturbation",
                      "Normalised Frequency (f/f_s)",
                      "Magnitude (dB)", s, 2);
    }
}

/* ================================================================== */
/*  Chapter 13: Spectral Analysis                                     */
/* ================================================================== */

static void plot_ch13(void)
{
    printf("  Ch13: spectral analysis ...\n");
    gp_init("ch13");

    /* Windowed vs unwindowed spectrum: 440 + 1000 + 2500 Hz */
    {
        int N = 256;
        double sig[256];
        double freqs[] = {440.0, 1000.0, 2500.0};
        double amps[]  = {1.0, 0.7, 0.4};
        gen_multi_tone(sig, N, freqs, amps, 3, 8000.0);

        /* Add noise */
        double noise[256];
        gen_white_noise(noise, N, 0.2, 42);
        for (int i = 0; i < N; i++) sig[i] += noise[i];

        int nb = N / 2 + 1;
        double freq_axis[129], mag_rect[129], mag_hann[129];

        /* Rectangular window */
        {
            Complex X[256];
            fft_real(sig, X, N);
            for (int k = 0; k < nb; k++) {
                freq_axis[k] = (double)k * 8000.0 / N;
                double m = complex_mag(X[k]) / (N / 2.0);
                mag_rect[k] = (m > 1e-10) ? 20.0 * log10(m) : -100.0;
            }
        }
        /* Hann window */
        {
            double windowed[256];
            for (int i = 0; i < N; i++) {
                double w = 0.5 * (1.0 - cos(2.0 * M_PI * i / (N - 1)));
                windowed[i] = sig[i] * w;
            }
            Complex X[256];
            fft_real(windowed, X, N);
            for (int k = 0; k < nb; k++) {
                double m = complex_mag(X[k]) / (N / 2.0);
                mag_hann[k] = (m > 1e-10) ? 20.0 * log10(m) : -100.0;
            }
        }

        GpSeries s[] = {
            { "Rectangular (no window)", freq_axis, mag_rect, nb, "lines" },
            { "Hann Window",             freq_axis, mag_hann, nb, "lines" },
        };
        FILE *gp = gp_open("ch13", "windowed_spectrum", 900, 500);
        if (gp) {
            fprintf(gp, "set title 'Spectral Analysis: Rectangular vs "
                        "Hann Window'\n");
            fprintf(gp, "set xlabel 'Frequency (Hz)'\n");
            fprintf(gp, "set ylabel 'Magnitude (dB)'\n");
            fprintf(gp, "set yrange [-80:5]\n");
            fprintf(gp, "set xrange [0:4000]\n");
            fprintf(gp, "plot ");
            for (int i = 0; i < 2; i++) {
                if (i) fprintf(gp, ", ");
                fprintf(gp, "'-' w lines lw 2 title '%s'", s[i].label);
            }
            fprintf(gp, "\n");
            for (int i = 0; i < 2; i++)
                gp_send_xy(gp, s[i].x, s[i].y, s[i].n);
            gp_close(gp);
        }
    }
}

/* ================================================================== */
/*  Chapter 30: Putting It All Together (Capstone)                    */
/* ================================================================== */

static void plot_ch30(void)
{
    printf("  Ch30: capstone pipeline ...\n");
    gp_init("ch30");

    int N = 512;
    double clean[512], noisy[512], noise[512], filtered[512];

    /* Signal: 200 + 500 Hz, noise at 2800 + 3500 Hz */
    double sig_f[] = {200.0, 500.0};
    double sig_a[] = {1.0, 0.6};
    gen_multi_tone(clean, N, sig_f, sig_a, 2, 8000.0);

    double nf[] = {2800.0, 3500.0};
    double na[] = {0.5, 0.3};
    gen_multi_tone(noise, N, nf, na, 2, 8000.0);
    for (int i = 0; i < N; i++) noisy[i] = clean[i] + noise[i];

    /* Filter: 31-tap lowpass at 800 Hz → cutoff = 800/4000 = 0.2 */
    int ntaps = 31;
    double h[31];
    fir_lowpass(h, ntaps, 0.2);
    fir_filter(noisy, filtered, N, h, ntaps);

    /* 1. Time domain: clean / noisy / filtered */
    {
        /* Show first 200 samples for clarity */
        GpSeries s[] = {
            { "Clean",    NULL, clean,    200, "lines" },
            { "Noisy",    NULL, noisy,    200, "lines" },
            { "Filtered", NULL, filtered, 200, "lines" },
        };
        gp_plot_multi("ch30", "pipeline_time",
                      "End-to-End Pipeline: Time Domain",
                      "Sample n", "Amplitude", s, 3);
    }

    /* 2. Spectrum: before and after filtering */
    {
        int nb = N / 2 + 1;
        double freq_axis[257], mag_before[257], mag_after[257];

        Complex Xb[512], Xa[512];
        fft_real(noisy, Xb, N);
        fft_real(filtered, Xa, N);

        for (int k = 0; k < nb; k++) {
            freq_axis[k] = (double)k * 8000.0 / N;
            double mb = complex_mag(Xb[k]) / (N / 2.0);
            double ma = complex_mag(Xa[k]) / (N / 2.0);
            mag_before[k] = (mb > 1e-10) ? 20.0 * log10(mb) : -100.0;
            mag_after[k]  = (ma > 1e-10) ? 20.0 * log10(ma) : -100.0;
        }

        GpSeries s[] = {
            { "Before Filtering", freq_axis, mag_before, nb, "lines" },
            { "After Filtering",  freq_axis, mag_after,  nb, "lines" },
        };
        FILE *gp = gp_open("ch30", "pipeline_spectrum", 900, 500);
        if (gp) {
            fprintf(gp, "set title 'Spectrum: Before vs After Lowpass "
                        "Filtering'\n");
            fprintf(gp, "set xlabel 'Frequency (Hz)'\n");
            fprintf(gp, "set ylabel 'Magnitude (dB)'\n");
            fprintf(gp, "set yrange [-80:5]\n");
            fprintf(gp, "set xrange [0:4000]\n");
            fprintf(gp, "plot ");
            for (int i = 0; i < 2; i++) {
                if (i) fprintf(gp, ", ");
                fprintf(gp, "'-' w lines lw 2 title '%s'", s[i].label);
            }
            fprintf(gp, "\n");
            for (int i = 0; i < 2; i++)
                gp_send_xy(gp, s[i].x, s[i].y, s[i].n);
            gp_close(gp);
        }
    }
}

/* ================================================================== */
/*  Ch14: PSD & Welch's Method plots                                  */
/*                                                                    */
/*  Generates:                                                        */
/*    - Periodogram vs Welch comparison                               */
/*    - Resolution trade-off (short/medium/long segments)             */
/* ================================================================== */

static void plot_ch14(void)
{
    printf("Ch14 PSD/Welch ...\n");

    const int    N  = 4096;
    const double fs = 8000.0;
    const double f1 = 500.0;

    double *x     = (double *)malloc((size_t)N * sizeof(double));
    double *noise = (double *)malloc((size_t)N * sizeof(double));

    gen_sine(x, N, 1.0, f1, fs, 0.0);
    gen_gaussian_noise(noise, N, 0.0, 2.0, 42);
    signal_add(x, noise, N);

    /* --- Periodogram --- */
    int nfft     = 4096;
    int n_bins_p = nfft / 2 + 1;
    double *psd_p   = (double *)calloc((size_t)n_bins_p, sizeof(double));
    double *db_p    = (double *)malloc((size_t)n_bins_p * sizeof(double));
    double *freq_p  = (double *)malloc((size_t)n_bins_p * sizeof(double));

    periodogram(x, N, psd_p, nfft);
    psd_to_db(psd_p, db_p, n_bins_p, -120.0);
    psd_freq_axis(freq_p, n_bins_p, fs);

    /* --- Welch --- */
    int seg_len  = 512;
    int nfft_w   = 512;
    int overlap  = 256;
    int n_bins_w = nfft_w / 2 + 1;
    double *psd_w  = (double *)calloc((size_t)n_bins_w, sizeof(double));
    double *db_w   = (double *)malloc((size_t)n_bins_w * sizeof(double));
    double *freq_w = (double *)malloc((size_t)n_bins_w * sizeof(double));

    welch_psd(x, N, psd_w, nfft_w, seg_len, overlap, hann_window);
    psd_to_db(psd_w, db_w, n_bins_w, -120.0);
    psd_freq_axis(freq_w, n_bins_w, fs);

    /* Plot both on one graph */
    GpSeries s[2];
    s[0].label = "Periodogram (high variance)";
    s[0].x     = freq_p;
    s[0].y     = db_p;
    s[0].n     = n_bins_p;
    s[0].style = "lines";

    s[1].label = "Welch (512-pt, 50%% overlap)";
    s[1].x     = freq_w;
    s[1].y     = db_w;
    s[1].n     = n_bins_w;
    s[1].style = "lines";

    gp_plot_multi("ch14", "periodogram_vs_welch",
                  "Periodogram vs Welch PSD (500 Hz + Noise)",
                  "Frequency (Hz)", "PSD (dB)",
                  s, 2);

    /* --- Resolution trade-off --- */
    int seg_lens[3]     = {128, 512, 2048};
    const char *labs[3] = {"128-pt", "512-pt", "2048-pt"};

    FILE *gp = gp_open("ch14", "resolution_tradeoff", 900, 500);
    if (gp) {
        fprintf(gp, "set title 'Welch PSD — Segment Length vs Resolution'\n");
        fprintf(gp, "set xlabel 'Frequency (Hz)'\n");
        fprintf(gp, "set ylabel 'PSD (dB)'\n");
        fprintf(gp, "set grid\n");
        fprintf(gp, "set xrange [0:2000]\n");
        fprintf(gp, "plot ");
        for (int si = 0; si < 3; si++) {
            if (si > 0) fprintf(gp, ", ");
            fprintf(gp, "'-' using 1:2 with lines lw 2 title '%s'", labs[si]);
        }
        fprintf(gp, "\n");

        /* Two close tones for resolution demo */
        double *x2    = (double *)malloc((size_t)N * sizeof(double));
        double *n2    = (double *)malloc((size_t)N * sizeof(double));
        double tones[2] = {900.0, 1100.0};
        double tamps[2] = {1.0, 1.0};
        gen_multi_tone(x2, N, tones, tamps, 2, fs);
        gen_gaussian_noise(n2, N, 0.0, 1.0, 99);
        signal_add(x2, n2, N);

        for (int si = 0; si < 3; si++) {
            int sl = seg_lens[si];
            int nf = next_power_of_2(sl);
            int nb = nf / 2 + 1;
            int ov = sl / 2;

            double *pw = (double *)calloc((size_t)nb, sizeof(double));
            double *dw = (double *)malloc((size_t)nb * sizeof(double));
            double *fw = (double *)malloc((size_t)nb * sizeof(double));

            welch_psd(x2, N, pw, nf, sl, ov, hann_window);
            psd_to_db(pw, dw, nb, -120.0);
            psd_freq_axis(fw, nb, fs);

            for (int k = 0; k < nb; k++)
                fprintf(gp, "%.2f %.4f\n", fw[k], dw[k]);
            fprintf(gp, "e\n");

            free(pw); free(dw); free(fw);
        }
        gp_close(gp);
        free(x2); free(n2);
    }

    free(x); free(noise);
    free(psd_p); free(db_p); free(freq_p);
    free(psd_w); free(db_w); free(freq_w);
}

/* ================================================================== */
/*  Ch15: Correlation plots                                           */
/*                                                                    */
/*  Generates:                                                        */
/*    - Autocorrelation of periodic signal (pitch detection)          */
/*    - White noise autocorrelation (delta)                           */
/* ================================================================== */

static void plot_ch15(void)
{
    printf("Ch15 Correlation ...\n");

    /* --- Autocorrelation for pitch: harmonic signal at 440 Hz --- */
    const int    N  = 2048;
    const double fs = 16000.0;
    const double f0 = 440.0;

    double *x     = (double *)malloc((size_t)N * sizeof(double));
    double *noise = (double *)malloc((size_t)N * sizeof(double));

    double freqs[3] = {f0, 2*f0, 3*f0};
    double amps[3]  = {1.0, 0.5, 0.25};
    gen_multi_tone(x, N, freqs, amps, 3, fs);
    gen_gaussian_noise(noise, N, 0.0, 0.3, 55);
    signal_add(x, noise, N);

    int r_len = 2 * N - 1;
    double *r = (double *)malloc((size_t)r_len * sizeof(double));
    autocorr_normalized(x, N, r);

    int centre   = N - 1;
    int max_lag  = (int)(fs / 50.0);
    int plot_len = max_lag + 50;
    if (plot_len > N) plot_len = N;

    double *lags = (double *)malloc((size_t)plot_len * sizeof(double));
    double *rpos = (double *)malloc((size_t)plot_len * sizeof(double));
    for (int i = 0; i < plot_len; i++) {
        lags[i] = (double)i;
        rpos[i] = r[centre + i];
    }

    gp_plot_1("ch15", "autocorr_pitch",
              "Autocorrelation — Pitch Detection (A4 = 440 Hz)",
              "Lag (samples)", "Normalised Autocorrelation",
              lags, rpos, plot_len, "lines");

    free(x); free(noise); free(r); free(lags); free(rpos);

    /* --- White noise autocorrelation — impulse at lag 0 --- */
    const int Nw = 4096;
    double *w = (double *)malloc((size_t)Nw * sizeof(double));
    gen_gaussian_noise(w, Nw, 0.0, 1.0, 123);

    int rw_len = 2 * Nw - 1;
    double *rw = (double *)malloc((size_t)rw_len * sizeof(double));
    autocorr_normalized(w, Nw, rw);

    int cw = Nw - 1;
    int ph = 100;
    int pw_len = 2 * ph + 1;
    double *wlags = (double *)malloc((size_t)pw_len * sizeof(double));
    double *wrp   = (double *)malloc((size_t)pw_len * sizeof(double));
    for (int i = 0; i < pw_len; i++) {
        wlags[i] = (double)(i - ph);
        wrp[i]   = rw[cw + i - ph];
    }

    gp_plot_1("ch15", "noise_autocorr",
              "Autocorrelation of White Noise — Impulse at Lag 0",
              "Lag (samples)", "Normalised Autocorrelation",
              wlags, wrp, pw_len, "impulses");

    free(w); free(rw); free(wlags); free(wrp);
}

/* ================================================================== */
/*  Chapter 16: Overlap-Add/Save Streaming                            */
/* ================================================================== */

static void plot_ch16(void)
{
    printf("  Ch16: OLA streaming ...\n");

    const int N = 1024, taps = 63, blk = 128;
    double h[63];
    fir_lowpass(h, taps, 0.2);

    double *x = (double *)malloc((size_t)N * sizeof(double));
    gen_sine(x, N, 1.0, 300.0, 8000.0, 0.0);

    /* Direct FIR */
    double *y_dir = (double *)calloc((size_t)N, sizeof(double));
    fir_filter(x, y_dir, N, h, taps);

    /* OLA */
    OlaState ola;
    ola_init(&ola, h, taps, blk);
    double *y_ola = (double *)calloc((size_t)N, sizeof(double));
    for (int b = 0; b < N / blk; b++)
        ola_process(&ola, x + b * blk, y_ola + b * blk);
    ola_free(&ola);

    /* Error */
    double *err = (double *)malloc((size_t)N * sizeof(double));
    double *idx = (double *)malloc((size_t)N * sizeof(double));
    for (int i = 0; i < N; i++) {
        err[i] = fabs(y_dir[i] - y_ola[i]);
        idx[i] = (double)i;
    }

    gp_plot_multi("ch16", "ola_vs_direct",
              "OLA vs Direct FIR Convolution",
              "Sample", "Amplitude",
              (GpSeries[]){
                  {"Direct FIR", idx, y_dir, N, "lines"},
                  {"OLA",        idx, y_ola, N, "lines"}
              }, 2);

    gp_plot_1("ch16", "ola_error",
              "OLA Reconstruction Error",
              "Sample", "|Error|",
              idx, err, N, "lines");

    free(x); free(y_dir); free(y_ola); free(err); free(idx);
}

/* ================================================================== */
/*  Chapter 18: Fixed-Point Arithmetic                                */
/* ================================================================== */

static void plot_ch18(void)
{
    printf("  Ch18: fixed-point SQNR ...\n");

    const int N = 512;
    double *x = (double *)malloc((size_t)N * sizeof(double));
    gen_sine(x, N, 0.9, 440.0, 8000.0, 0.0);

    /* Q15 quantisation error */
    q15_t *xq = (q15_t *)malloc((size_t)N * sizeof(q15_t));
    double *xr = (double *)malloc((size_t)N * sizeof(double));
    double *qerr = (double *)malloc((size_t)N * sizeof(double));
    double *idx  = (double *)malloc((size_t)N * sizeof(double));

    double_array_to_q15(x, xq, N);
    q15_array_to_double(xq, xr, N);

    for (int i = 0; i < N; i++) {
        qerr[i] = x[i] - xr[i];
        idx[i] = (double)i;
    }

    gp_plot_multi("ch18", "q15_quantisation",
              "Q15 Quantisation: Signal vs Recovered",
              "Sample", "Amplitude",
              (GpSeries[]){
                  {"Original (float)", idx, x,  N, "lines"},
                  {"Q15 recovered",    idx, xr, N, "lines"}
              }, 2);

    gp_plot_1("ch18", "q15_error",
              "Q15 Quantisation Error",
              "Sample", "Error",
              idx, qerr, N, "lines");

    free(x); free(xq); free(xr); free(qerr); free(idx);
}

/* ================================================================== */
/*  Chapter 19: Advanced FFT — Goertzel & Sliding DFT                 */
/* ================================================================== */

static void plot_ch19(void)
{
    printf("  Ch19: Goertzel spectrum ...\n");

    /* Goertzel power spectrum of a DTMF tone ('5': 770+1336 Hz) */
    const int N = 512;
    const double fs = 8000.0;
    double *tone = (double *)malloc((size_t)N * sizeof(double));
    double *tmp  = (double *)malloc((size_t)N * sizeof(double));
    gen_sine(tone, N, 0.5, 770.0,  fs, 0.0);
    gen_sine(tmp,  N, 0.5, 1336.0, fs, 0.0);
    signal_add(tone, tmp, N);

    /* Scan 0–4000 Hz in 10 Hz steps with Goertzel */
    int n_pts = 400;
    double *freq = (double *)malloc((size_t)n_pts * sizeof(double));
    double *power = (double *)malloc((size_t)n_pts * sizeof(double));
    for (int i = 0; i < n_pts; i++) {
        freq[i] = (double)i * 10.0;
        Complex G = goertzel_freq(tone, N, freq[i], fs);
        power[i] = 10.0 * log10(G.re * G.re + G.im * G.im + 1e-30);
    }

    gp_plot_1("ch19", "goertzel_dtmf",
              "Goertzel Spectrum - DTMF 5 (770+1336 Hz)",
              "Frequency (Hz)", "Power (dB)",
              freq, power, n_pts, "lines");

    free(tone); free(tmp); free(freq); free(power);
}

/* ================================================================== */
/*  Chapter 17: Multirate DSP                                         */
/* ================================================================== */

static void plot_ch17(void)
{
    printf("  Ch17: multirate decimation / interpolation ...\n");

    /* Original signal: 300 Hz sine at 8 kHz */
    const int N = 256;
    const double fs = 8000.0;
    double *x   = (double *)malloc((size_t)N * sizeof(double));
    double *idx = (double *)malloc((size_t)N * sizeof(double));
    gen_sine(x, N, 1.0, 300.0, fs, 0.0);
    for (int i = 0; i < N; i++) idx[i] = (double)i;

    /* Decimate by 4 */
    const int M = 4;
    int out_len = N / M;
    double *y_dec  = (double *)malloc((size_t)out_len * sizeof(double));
    double *idx_dec = (double *)malloc((size_t)out_len * sizeof(double));
    decimate(x, N, M, y_dec);
    for (int i = 0; i < out_len; i++) idx_dec[i] = (double)(i * M);

    gp_plot_multi("ch17", "decimation",
              "Decimation by 4: 300 Hz @ 8 kHz",
              "Sample", "Amplitude",
              (GpSeries[]){
                  {"Original",    idx,     x,     N,       "lines"},
                  {"Decimated×4", idx_dec, y_dec, out_len, "points"}
              }, 2);

    /* Interpolate by 4 */
    const int L = 4;
    int int_len = N * L;
    double *y_int  = (double *)malloc((size_t)int_len * sizeof(double));
    double *idx_int = (double *)malloc((size_t)int_len * sizeof(double));
    interpolate(x, N, L, y_int);
    for (int i = 0; i < int_len; i++) idx_int[i] = (double)i / (double)L;

    gp_plot_multi("ch17", "interpolation",
              "Interpolation by 4: 300 Hz @ 8 kHz",
              "Sample (original rate)", "Amplitude",
              (GpSeries[]){
                  {"Original",       idx,     x,     N,       "points"},
                  {"Interpolated×4", idx_int, y_int, int_len, "lines"}
              }, 2);

    /* Polyphase vs direct decimation — magnitude spectrum comparison */
    const int nfft = 512;
    double h_aa[33];
    fir_lowpass(h_aa, 33, 0.5 / M);

    double *y_poly = (double *)malloc((size_t)out_len * sizeof(double));
    polyphase_decimate(x, N, h_aa, 33, M, y_poly);

    /* FFT of polyphase output */
    Complex *H_poly = (Complex *)calloc((size_t)nfft, sizeof(Complex));
    for (int i = 0; i < out_len && i < nfft; i++) {
        H_poly[i].re = y_poly[i];
        H_poly[i].im = 0.0;
    }
    fft(H_poly, nfft);

    double *freq_ax = (double *)malloc((size_t)(nfft / 2) * sizeof(double));
    double *mag_poly = (double *)malloc((size_t)(nfft / 2) * sizeof(double));
    for (int i = 0; i < nfft / 2; i++) {
        freq_ax[i] = (double)i / (double)nfft;
        mag_poly[i] = 10.0 * log10(complex_mag(H_poly[i]) *
                                     complex_mag(H_poly[i]) + 1e-30);
    }

    gp_plot_1("ch17", "polyphase_spectrum",
              "Polyphase Decimated Spectrum",
              "Normalised Frequency", "Magnitude (dB)",
              freq_ax, mag_poly, nfft / 2, "lines");

    free(x); free(idx); free(y_dec); free(idx_dec);
    free(y_int); free(idx_int); free(y_poly);
    free(H_poly); free(freq_ax); free(mag_poly);
}

/* ================================================================== */
/*  Chapter 20: Hilbert Transform & Analytic Signal                   */
/* ================================================================== */

static void plot_ch20(void)
{
    printf("  Ch20: Hilbert analytic signal ...\n");

    /* AM-modulated signal: envelope = 1 + 0.5*cos(2π·5·t) */
    const int N = 512;
    const double fs = 1000.0;
    double *x   = (double *)malloc((size_t)N * sizeof(double));
    double *env_true  = (double *)malloc((size_t)N * sizeof(double));
    double *env_hilb  = (double *)malloc((size_t)N * sizeof(double));
    double *idx = (double *)malloc((size_t)N * sizeof(double));

    for (int i = 0; i < N; i++) {
        double t = (double)i / fs;
        double env_val = 1.0 + 0.5 * cos(2.0 * M_PI * 5.0 * t);
        x[i] = env_val * cos(2.0 * M_PI * 100.0 * t);
        env_true[i] = env_val;
        idx[i] = t * 1000.0;  /* ms */
    }

    envelope(x, N, env_hilb, 65);

    gp_plot_multi("ch20", "am_envelope",
              "Hilbert Envelope Detection — AM Signal (fc=100, fm=5 Hz)",
              "Time (ms)", "Amplitude",
              (GpSeries[]){
                  {"AM Signal",       idx, x,        N, "lines"},
                  {"True Envelope",   idx, env_true, N, "lines"},
                  {"Hilbert Envelope", idx, env_hilb, N, "lines"}
              }, 3);

    /* Instantaneous frequency of a chirp */
    double *chirp  = (double *)malloc((size_t)N * sizeof(double));
    double *ifreq  = (double *)malloc((size_t)N * sizeof(double));
    gen_chirp(chirp, N, 1.0, 50.0, 200.0, fs);
    inst_frequency(chirp, N, ifreq, 65);

    /* Scale to Hz */
    for (int i = 0; i < N; i++) {
        ifreq[i] *= fs;
        idx[i] = (double)i / fs * 1000.0;
    }

    gp_plot_multi("ch20", "inst_frequency",
              "Instantaneous Frequency — Linear Chirp (50→200 Hz)",
              "Time (ms)", "Frequency / Amplitude",
              (GpSeries[]){
                  {"Chirp Signal", idx, chirp, N, "lines"},
                  {"Inst. Freq (Hz)", idx, ifreq, N, "lines"}
              }, 2);

    free(x); free(env_true); free(env_hilb); free(idx);
    free(chirp); free(ifreq);
}

/* ================================================================== */
/*  Chapter 21: Signal Averaging & Noise Reduction                    */
/* ================================================================== */

static void plot_ch21(void)
{
    printf("  Ch21: signal averaging ...\n");

    const int N = 256;
    const int K = 32;

    /* Generate K noisy trials of a square pulse */
    double *clean = (double *)calloc((size_t)N, sizeof(double));
    for (int i = N / 4; i < 3 * N / 4; i++) clean[i] = 1.0;

    double **trials = (double **)malloc((size_t)K * sizeof(double *));
    for (int k = 0; k < K; k++) {
        trials[k] = (double *)malloc((size_t)N * sizeof(double));
        gen_white_noise(trials[k], N, 0.5, (unsigned)(42 + k));
        signal_add(trials[k], clean, N);
    }

    /* Coherent average */
    double *avg = (double *)malloc((size_t)N * sizeof(double));
    coherent_average((const double **)trials, K, N, avg);

    double *idx = (double *)malloc((size_t)N * sizeof(double));
    for (int i = 0; i < N; i++) idx[i] = (double)i;

    gp_plot_multi("ch21", "coherent_averaging",
              "Coherent Averaging (32 trials, σ_noise=0.5)",
              "Sample", "Amplitude",
              (GpSeries[]){
                  {"Noisy (trial 1)", idx, trials[0], N, "lines"},
                  {"Clean Signal",    idx, clean,      N, "lines"},
                  {"32-trial Average", idx, avg,       N, "lines"}
              }, 3);

    /* EMA vs MA on noisy signal */
    double *noisy = (double *)malloc((size_t)N * sizeof(double));
    gen_white_noise(noisy, N, 0.3, 99);
    signal_add(noisy, clean, N);

    double *y_ema = (double *)malloc((size_t)N * sizeof(double));
    double *y_ma  = (double *)malloc((size_t)N * sizeof(double));
    ema_filter(noisy, N, 0.1, y_ema);
    moving_average(noisy, N, 15, y_ma);

    gp_plot_multi("ch21", "ema_vs_ma",
              "EMA (α=0.1) vs Moving Average (M=15)",
              "Sample", "Amplitude",
              (GpSeries[]){
                  {"Noisy Input", idx, noisy, N, "lines"},
                  {"EMA α=0.1",   idx, y_ema, N, "lines"},
                  {"MA M=15",     idx, y_ma,  N, "lines"}
              }, 3);

    /* Median filter on impulse-corrupted signal */
    double *corrupted = (double *)malloc((size_t)N * sizeof(double));
    for (int i = 0; i < N; i++) corrupted[i] = clean[i];
    /* Add salt-and-pepper impulse noise */
    for (int i = 0; i < N; i += 7 + (i % 3))
        corrupted[i] += (i % 2 == 0) ? 2.0 : -2.0;

    double *y_med = (double *)malloc((size_t)N * sizeof(double));
    median_filter(corrupted, N, 5, y_med);

    gp_plot_multi("ch21", "median_filter",
              "Median Filter (M=5) — Impulse Noise Removal",
              "Sample", "Amplitude",
              (GpSeries[]){
                  {"Corrupted",      idx, corrupted, N, "lines"},
                  {"Median Filtered", idx, y_med,    N, "lines"},
                  {"Clean",          idx, clean,     N, "lines"}
              }, 3);

    for (int k = 0; k < K; k++) free(trials[k]);
    free(trials); free(clean); free(avg); free(idx);
    free(noisy); free(y_ema); free(y_ma);
    free(corrupted); free(y_med);
}

/* ================================================================== */
/*  Chapter 22: Advanced FIR Design (Remez / IRLS)                    */
/* ================================================================== */

static void plot_ch22(void)
{
    printf("  Ch22: Remez FIR design ...\n");

    const int taps = 51;
    const int nfft = 1024;

    /* Remez lowpass: passband 0–0.2, stopband 0.3–0.5 */
    double h_remez[51];
    remez_lowpass(h_remez, taps, 0.2, 0.3, 1.0, 1.0);

    /* Window-method lowpass for comparison */
    double h_win[51];
    fir_lowpass(h_win, taps, 0.25);

    /* Frequency responses */
    Complex *H_r = (Complex *)calloc((size_t)nfft, sizeof(Complex));
    Complex *H_w = (Complex *)calloc((size_t)nfft, sizeof(Complex));
    for (int i = 0; i < taps; i++) {
        H_r[i].re = h_remez[i]; H_r[i].im = 0.0;
        H_w[i].re = h_win[i];   H_w[i].im = 0.0;
    }
    fft(H_r, nfft);
    fft(H_w, nfft);

    int half = nfft / 2;
    double *freq  = (double *)malloc((size_t)half * sizeof(double));
    double *mag_r = (double *)malloc((size_t)half * sizeof(double));
    double *mag_w = (double *)malloc((size_t)half * sizeof(double));

    for (int i = 0; i < half; i++) {
        freq[i] = (double)i / (double)nfft;
        double mr = complex_mag(H_r[i]);
        double mw = complex_mag(H_w[i]);
        mag_r[i] = 20.0 * log10(mr + 1e-30);
        mag_w[i] = 20.0 * log10(mw + 1e-30);
    }

    gp_plot_multi("ch22", "remez_vs_window",
              "Remez vs Window-Method Lowpass (51 taps, fc≈0.25)",
              "Normalised Frequency (0–0.5)", "Magnitude (dB)",
              (GpSeries[]){
                  {"Remez (pb=0.2, sb=0.3)", freq, mag_r, half, "lines"},
                  {"Window (Hann, fc=0.25)",  freq, mag_w, half, "lines"}
              }, 2);

    /* Remez bandpass filter */
    double h_bp[51];
    remez_bandpass(h_bp, taps, 0.1, 0.2, 0.35, 0.45);

    Complex *H_bp = (Complex *)calloc((size_t)nfft, sizeof(Complex));
    for (int i = 0; i < taps; i++) {
        H_bp[i].re = h_bp[i]; H_bp[i].im = 0.0;
    }
    fft(H_bp, nfft);

    double *mag_bp = (double *)malloc((size_t)half * sizeof(double));
    for (int i = 0; i < half; i++) {
        double m = complex_mag(H_bp[i]);
        mag_bp[i] = 20.0 * log10(m + 1e-30);
    }

    gp_plot_1("ch22", "remez_bandpass",
              "Remez Bandpass FIR (51 taps, pass=0.2–0.35)",
              "Normalised Frequency (0–0.5)", "Magnitude (dB)",
              freq, mag_bp, half, "lines");

    free(H_r); free(H_w); free(H_bp);
    free(freq); free(mag_r); free(mag_w); free(mag_bp);
}

/* ================================================================== */
/*  Chapter 23: Adaptive Filters                                      */
/* ================================================================== */

static void plot_ch23(void)
{
    printf("  Ch23: adaptive filters ...\n");
    gp_init("ch23");

    /* LMS system identification — learning curve */
    {
        int N = 1000, L = 8;
        double plant[8] = {0.1, -0.3, 0.5, -0.2, 0.4, -0.1, 0.3, -0.15};
        double *x = (double *)malloc((size_t)N * sizeof(double));
        double *d = (double *)malloc((size_t)N * sizeof(double));
        double *y = (double *)malloc((size_t)N * sizeof(double));
        double *e = (double *)malloc((size_t)N * sizeof(double));
        double w[8];
        gen_white_noise(x, N, 1.0, 42);
        /* Generate desired = x filtered through plant */
        for (int i = 0; i < N; i++) {
            d[i] = 0;
            for (int j = 0; j < L; j++)
                if (i - j >= 0) d[i] += plant[j] * x[i - j];
        }
        lms_filter(x, d, N, L, 0.01, y, e, w);

        /* Plot squared error (smoothed) */
        double *e2 = (double *)malloc((size_t)N * sizeof(double));
        for (int i = 0; i < N; i++) e2[i] = e[i] * e[i];
        /* Moving average smoothing (window=20) */
        double *e2s = (double *)calloc((size_t)N, sizeof(double));
        int win = 20;
        for (int i = 0; i < N; i++) {
            int cnt = 0;
            for (int j = i - win + 1; j <= i; j++)
                if (j >= 0) { e2s[i] += e2[j]; cnt++; }
            e2s[i] /= cnt;
        }
        gp_plot_1("ch23", "lms_learning",
                  "LMS Learning Curve (System Identification)",
                  "Iteration", "MSE (smoothed)", NULL, e2s, N, "lines");
        free(x); free(d); free(y); free(e); free(e2); free(e2s);
    }

    /* LMS vs NLMS convergence comparison */
    {
        int N = 500, L = 4;
        double plant[4] = {0.5, -0.3, 0.2, -0.1};
        double *x = (double *)malloc((size_t)N * sizeof(double));
        double *d = (double *)malloc((size_t)N * sizeof(double));
        double *y_l = (double *)malloc((size_t)N * sizeof(double));
        double *y_n = (double *)malloc((size_t)N * sizeof(double));
        double *e_l = (double *)malloc((size_t)N * sizeof(double));
        double *e_n = (double *)malloc((size_t)N * sizeof(double));
        double w_l[4], w_n[4];
        gen_white_noise(x, N, 1.0, 99);
        for (int i = 0; i < N; i++) {
            d[i] = 0;
            for (int j = 0; j < L; j++)
                if (i - j >= 0) d[i] += plant[j] * x[i - j];
        }
        lms_filter(x, d, N, L, 0.02, y_l, e_l, w_l);
        nlms_filter(x, d, N, L, 0.5, 1e-8, y_n, e_n, w_n);

        double *el2 = (double *)malloc((size_t)N * sizeof(double));
        double *en2 = (double *)malloc((size_t)N * sizeof(double));
        for (int i = 0; i < N; i++) {
            el2[i] = e_l[i] * e_l[i];
            en2[i] = e_n[i] * e_n[i];
        }
        GpSeries s[] = {
            { "LMS",  NULL, el2, N, "lines" },
            { "NLMS", NULL, en2, N, "lines" },
        };
        gp_plot_multi("ch23", "lms_vs_nlms",
                      "LMS vs NLMS Convergence",
                      "Iteration", "Squared Error", s, 2);
        free(x); free(d); free(y_l); free(y_n);
        free(e_l); free(e_n); free(el2); free(en2);
    }
}

/* ================================================================== */
/*  Chapter 24: Linear Prediction                                     */
/* ================================================================== */

static void plot_ch24(void)
{
    printf("  Ch24: linear prediction ...\n");
    gp_init("ch24");

    /* AR spectrum at different model orders */
    {
        int N = 512, nfft = 512;
        double *x = (double *)malloc((size_t)N * sizeof(double));
        /* Two-tone signal */
        for (int i = 0; i < N; i++)
            x[i] = sin(2.0 * M_PI * 0.1 * i) + 0.7 * sin(2.0 * M_PI * 0.3 * i);

        int half = nfft / 2;
        double *freq = (double *)malloc((size_t)half * sizeof(double));
        for (int i = 0; i < half; i++)
            freq[i] = (double)i / (double)nfft;

        double *s4 = (double *)malloc((size_t)nfft * sizeof(double));
        double *s10 = (double *)malloc((size_t)nfft * sizeof(double));
        double *s20 = (double *)malloc((size_t)nfft * sizeof(double));
        {
            double a[4], E;
            lpc_coefficients(x, N, 4, a, &E);
            lpc_spectrum(a, 4, E, s4, nfft);
        }
        {
            double a[10], E;
            lpc_coefficients(x, N, 10, a, &E);
            lpc_spectrum(a, 10, E, s10, nfft);
        }
        {
            double a[20], E;
            lpc_coefficients(x, N, 20, a, &E);
            lpc_spectrum(a, 20, E, s20, nfft);
        }
        /* Convert to dB */
        double *db4 = (double *)malloc((size_t)half * sizeof(double));
        double *db10 = (double *)malloc((size_t)half * sizeof(double));
        double *db20 = (double *)malloc((size_t)half * sizeof(double));
        for (int i = 0; i < half; i++) {
            db4[i]  = 10.0 * log10(s4[i] + 1e-30);
            db10[i] = 10.0 * log10(s10[i] + 1e-30);
            db20[i] = 10.0 * log10(s20[i] + 1e-30);
        }
        GpSeries s[] = {
            { "Order 4",  freq, db4,  half, "lines" },
            { "Order 10", freq, db10, half, "lines" },
            { "Order 20", freq, db20, half, "lines" },
        };
        gp_plot_multi("ch24", "ar_spectrum",
                      "AR Spectral Envelope (Two-Tone Signal)",
                      "Normalised Frequency", "Power (dB)", s, 3);
        free(x); free(freq); free(s4); free(s10); free(s20);
        free(db4); free(db10); free(db20);
    }

    /* Levinson-Durbin: original vs reconstructed */
    {
        int N = 128, order = 10;
        double *x = (double *)malloc((size_t)N * sizeof(double));
        for (int i = 0; i < N; i++)
            x[i] = sin(2.0 * M_PI * 0.15 * i) + 0.3 * cos(2.0 * M_PI * 0.35 * i);

        double a[10], E;
        lpc_coefficients(x, N, order, a, &E);

        double *res = (double *)malloc((size_t)N * sizeof(double));
        double *rec = (double *)malloc((size_t)N * sizeof(double));
        lpc_residual(x, N, a, order, res);
        lpc_synthesise(res, N, a, order, rec);

        GpSeries s[] = {
            { "Original",      NULL, x,   N, "lines" },
            { "Reconstructed", NULL, rec, N, "lines" },
        };
        gp_plot_multi("ch24", "lpc_roundtrip",
                      "LPC Analysis/Synthesis Round-Trip (order=10)",
                      "Sample n", "x[n]", s, 2);
        free(x); free(res); free(rec);
    }
}

/* ================================================================== */
/*  Chapter 25: Parametric Spectral Estimation                        */
/* ================================================================== */

static void plot_ch25(void)
{
    printf("  Ch25: parametric spectral ...\n");
    gp_init("ch25");

    /* MUSIC super-resolution */
    {
        double f1 = 0.12, f2 = 0.14;  /* close frequencies */
        int N = 128;
        double *x = (double *)malloc((size_t)N * sizeof(double));
        unsigned int seed = 10;
        for (int i = 0; i < N; i++)
            x[i] = sin(2.0 * M_PI * f1 * i) + sin(2.0 * M_PI * f2 * i)
                  + 0.1 * ((double)(seed = seed * 1103515245 + 12345) / 2147483648.0 - 0.5);

        int nfft = 1024, half = nfft / 2;
        double *spec = (double *)malloc((size_t)half * sizeof(double));
        double *freq = (double *)malloc((size_t)half * sizeof(double));
        music_spectrum(x, N, 20, 2, spec, nfft);
        for (int i = 0; i < half; i++) {
            freq[i] = (double)i / (double)nfft;
            spec[i] = 10.0 * log10(spec[i] + 1e-30);
        }
        gp_plot_1("ch25", "music_spectrum",
                  "MUSIC Pseudospectrum (f1=0.12, f2=0.14)",
                  "Normalised Frequency", "MUSIC (dB)",
                  freq, spec, half, "lines");
        free(x); free(spec); free(freq);
    }

    /* Capon vs FFT */
    {
        double f1 = 0.15;
        int N = 256;
        double *x = (double *)malloc((size_t)N * sizeof(double));
        unsigned int seed2 = 77;
        for (int i = 0; i < N; i++)
            x[i] = sin(2.0 * M_PI * f1 * i) + sin(2.0 * M_PI * 0.35 * i)
                  + 0.1 * ((double)(seed2 = seed2 * 1103515245 + 12345) / 2147483648.0 - 0.5);

        int nfft = 512, half = nfft / 2;
        double *capon = (double *)malloc((size_t)half * sizeof(double));
        double *freq = (double *)malloc((size_t)half * sizeof(double));
        capon_spectrum(x, N, 16, capon, nfft);
        for (int i = 0; i < half; i++) {
            freq[i] = (double)i / (double)nfft;
            capon[i] = 10.0 * log10(capon[i] + 1e-30);
        }

        /* Also compute FFT-based periodogram for comparison */
        int nfft2 = 512;
        Complex *X = (Complex *)calloc((size_t)nfft2, sizeof(Complex));
        for (int i = 0; i < N && i < nfft2; i++) X[i].re = x[i];
        fft(X, nfft2);
        double *fft_db = (double *)malloc((size_t)half * sizeof(double));
        for (int i = 0; i < half; i++)
            fft_db[i] = 10.0 * log10(complex_mag(X[i]) * complex_mag(X[i]) / nfft2 + 1e-30);

        GpSeries s[] = {
            { "Capon (MVDR)", freq, capon,  half, "lines" },
            { "FFT Periodogram", freq, fft_db, half, "lines" },
        };
        gp_plot_multi("ch25", "capon_vs_fft",
                      "Capon vs FFT Spectrum",
                      "Normalised Frequency", "Power (dB)", s, 2);
        free(x); free(capon); free(freq); free(X); free(fft_db);
    }
}

/* ================================================================== */
/*  Chapter 26: Cepstrum & MFCC                                       */
/* ================================================================== */

static void plot_ch26(void)
{
    printf("  Ch26: cepstrum & MFCC ...\n");
    gp_init("ch26");

    /* Real cepstrum and liftered envelope */
    {
        int N = 512;
        double *x = (double *)malloc((size_t)N * sizeof(double));
        /* Simulate speech-like: periodic pulse + formant-like resonance */
        for (int i = 0; i < N; i++) {
            double pulse = (i % 50 == 0) ? 1.0 : 0.0;
            x[i] = pulse;
        }
        /* Simple formant filter: IIR low-order */
        for (int i = 1; i < N; i++)
            x[i] += 0.9 * x[i - 1];

        double *cep = (double *)malloc((size_t)N * sizeof(double));
        cepstrum_real(x, N, cep, N);

        gp_plot_1("ch26", "cepstrum_real",
                  "Real Cepstrum (Pulse Train + Formant Filter)",
                  "Quefrency (samples)", "Cepstral Value",
                  NULL, cep, N / 2, "lines");  /* show first half */

        /* Liftered envelope */
        double *liftered = (double *)malloc((size_t)N * sizeof(double));
        cepstrum_lifter(cep, N, 30, liftered);

        GpSeries s[] = {
            { "Full Cepstrum", NULL, cep,      N / 2, "lines" },
            { "Liftered (L=30)", NULL, liftered, N / 2, "lines" },
        };
        gp_plot_multi("ch26", "cepstrum_lifter",
                      "Cepstrum vs Liftered (Spectral Envelope)",
                      "Quefrency (samples)", "Value", s, 2);

        free(x); free(cep); free(liftered);
    }

    /* Mel filterbank shape */
    {
        double fs = 8000.0;
        int n_filters = 26;
        int nfft = 512, half = nfft / 2;
        double *freq = (double *)malloc((size_t)half * sizeof(double));
        for (int i = 0; i < half; i++)
            freq[i] = (double)i * fs / (double)nfft;

        FILE *gp = gp_open("ch26", "mel_filterbank", 800, 400);
        if (gp) {
            fprintf(gp, "set title 'Mel Filterbank (%d filters, fs=%.0f Hz)'\n",
                    n_filters, fs);
            fprintf(gp, "set xlabel 'Frequency (Hz)'\nset ylabel 'Weight'\n");
            fprintf(gp, "set grid\n");

            double mel_lo = hz_to_mel(0.0);
            double mel_hi = hz_to_mel(fs / 2.0);
            int n_pts = n_filters + 2;
            fprintf(gp, "plot ");
            for (int f = 0; f < n_filters; f++) {
                double c_lo = mel_to_hz(mel_lo + (mel_hi - mel_lo) * f / (n_pts - 1));
                double c_mid = mel_to_hz(mel_lo + (mel_hi - mel_lo) * (f + 1) / (n_pts - 1));
                double c_hi = mel_to_hz(mel_lo + (mel_hi - mel_lo) * (f + 2) / (n_pts - 1));
                if (f > 0) fprintf(gp, ", ");
                fprintf(gp, "'-' with lines notitle");
                (void)c_lo; (void)c_mid; (void)c_hi;
            }
            fprintf(gp, "\n");
            /* Send each triangular filter as 3 points */
            for (int f = 0; f < n_filters; f++) {
                double c_lo = mel_to_hz(mel_lo + (mel_hi - mel_lo) * f / (n_pts - 1));
                double c_mid = mel_to_hz(mel_lo + (mel_hi - mel_lo) * (f + 1) / (n_pts - 1));
                double c_hi = mel_to_hz(mel_lo + (mel_hi - mel_lo) * (f + 2) / (n_pts - 1));
                fprintf(gp, "%.2f 0\n%.2f 1\n%.2f 0\ne\n", c_lo, c_mid, c_hi);
            }
            gp_close(gp);
        }
        free(freq);
    }
}

/* ================================================================== */
/*  Chapter 27: 2-D DSP                                               */
/* ================================================================== */

static void plot_ch27(void)
{
    printf("  Ch27: 2D DSP ...\n");
    gp_init("ch27");

    /* 2D test image: gradient + circle */
    int rows = 64, cols = 64;
    double *img = (double *)calloc((size_t)(rows * cols), sizeof(double));
    for (int r = 0; r < rows; r++)
        for (int c = 0; c < cols; c++) {
            double cx = c - 32.0, cy = r - 32.0;
            double dist = sqrt(cx * cx + cy * cy);
            img[r * cols + c] = (dist < 15.0) ? 1.0 : (double)c / cols;
        }

    /* Gaussian blur */
    {
        double kernel[25];
        kernel_gaussian(kernel, 5, 1.5);
        double *out = (double *)malloc((size_t)(rows * cols) * sizeof(double));
        conv2d(img, rows, cols, kernel, 5, 5, out);

        FILE *gp = gp_open("ch27", "gaussian_blur", 600, 500);
        if (gp) {
            fprintf(gp, "set title 'Gaussian Blur (5x5, σ=1.5)'\n");
            fprintf(gp, "set pm3d map\nset palette grey\n");
            fprintf(gp, "set xlabel 'Column'\nset ylabel 'Row'\n");
            fprintf(gp, "splot '-' matrix with image\n");
            for (int r = 0; r < rows; r++) {
                for (int c = 0; c < cols; c++)
                    fprintf(gp, "%.4f ", out[r * cols + c]);
                fprintf(gp, "\n");
            }
            fprintf(gp, "e\ne\n");
            gp_close(gp);
        }
        free(out);
    }

    /* Sobel edges */
    {
        double *mag = (double *)malloc((size_t)(rows * cols) * sizeof(double));
        sobel_magnitude(img, rows, cols, mag);

        FILE *gp = gp_open("ch27", "sobel_edges", 600, 500);
        if (gp) {
            fprintf(gp, "set title 'Sobel Edge Detection'\n");
            fprintf(gp, "set pm3d map\nset palette grey\n");
            fprintf(gp, "set xlabel 'Column'\nset ylabel 'Row'\n");
            fprintf(gp, "splot '-' matrix with image\n");
            for (int r = 0; r < rows; r++) {
                for (int c = 0; c < cols; c++)
                    fprintf(gp, "%.4f ", mag[r * cols + c]);
                fprintf(gp, "\n");
            }
            fprintf(gp, "e\ne\n");
            gp_close(gp);
        }
        free(mag);
    }

    free(img);
}

/* ── Ch28: Real-Time Streaming ──────────────────────────────────── */

static void plot_ch28(void)
{
    printf("  Ch28: Real-Time Streaming ...\n");
    gp_init("ch28");

    /* Ring buffer fill level over time */
    {
        RingBuffer *rb = ring_buffer_create(256);
        int n_pts = 200;
        double x[200], y[200];
        unsigned int seed = 42;
        for (int i = 0; i < n_pts; i++) {
            /* Write 1-4 samples, read 1-3 */
            seed = seed * 1103515245 + 12345;
            int nw = 1 + (int)((seed >> 16) % 4);
            double tmp[4] = {0};
            ring_buffer_write(rb, tmp, nw);
            seed = seed * 1103515245 + 12345;
            int nr = 1 + (int)((seed >> 16) % 3);
            if (nr > ring_buffer_available(rb)) nr = ring_buffer_available(rb);
            double discard[4];
            ring_buffer_read(rb, discard, nr);
            x[i] = (double)i;
            y[i] = (double)ring_buffer_available(rb);
        }
        gp_plot_1("ch28", "ring_buffer_fill",
                  "Ring Buffer Fill Level Over Time",
                  "Iteration", "Samples in Buffer",
                  x, y, n_pts, "lines");
        ring_buffer_destroy(rb);
    }

    /* Streaming FFT: peak frequency per frame */
    {
        int sr = 8000, total = 4096;
        double *sig = (double *)malloc((size_t)total * sizeof(double));
        gen_chirp(sig, total, 1.0, 200.0, 2000.0, (double)sr);

        int frame_sz = 256, hop = 128;
        FrameProcessor *fp = frame_processor_create(frame_sz, hop);
        int max_frames = (total - frame_sz) / hop + 1;
        double *fx = (double *)malloc((size_t)max_frames * sizeof(double));
        double *fy = (double *)malloc((size_t)max_frames * sizeof(double));
        int fi = 0;
        for (int i = 0; i < total; i++) {
            if (frame_processor_feed(fp, &sig[i], 1)) {
                fx[fi] = (double)fi;
                fy[fi] = frame_processor_peak_freq(fp, (double)sr);
                fi++;
            }
        }
        gp_plot_1("ch28", "streaming_fft_peak",
                  "Streaming FFT: Peak Frequency per Frame",
                  "Frame", "Frequency (Hz)",
                  fx, fy, fi, "linespoints");
        frame_processor_destroy(fp);
        free(sig); free(fx); free(fy);
    }
}

/* ── Ch29: Optimisation ─────────────────────────────────────────── */

static void plot_ch29(void)
{
    printf("  Ch29: Optimisation ...\n");
    gp_init("ch29");

    /* Benchmark: radix-2 vs radix-4 throughput across sizes */
    {
        int sizes[] = { 64, 256, 1024, 4096 };
        int nsizes = 4, runs = 50;
        double xs[4], y_r2[4], y_r4[4];
        for (int i = 0; i < nsizes; i++) {
            BenchResult r2 = bench_fft_radix2(sizes[i], runs);
            BenchResult r4 = bench_fft_radix4(sizes[i], runs);
            xs[i]   = log2((double)sizes[i]);
            y_r2[i] = r2.mflops;
            y_r4[i] = r4.mflops;
        }
        FILE *gp = gp_open("ch29", "radix_comparison", 640, 480);
        if (gp) {
            fprintf(gp, "set title 'FFT Throughput: Radix-2 vs Radix-4'\n");
            fprintf(gp, "set xlabel 'log2(N)'\nset ylabel 'MFLOP/s'\n");
            fprintf(gp, "set grid\nset key top left\n");
            fprintf(gp, "plot '-' with linespoints title 'Radix-2', "
                        "'-' with linespoints title 'Radix-4'\n");
            for (int i = 0; i < nsizes; i++)
                fprintf(gp, "%.1f %.2f\n", xs[i], y_r2[i]);
            fprintf(gp, "e\n");
            for (int i = 0; i < nsizes; i++)
                fprintf(gp, "%.1f %.2f\n", xs[i], y_r4[i]);
            fprintf(gp, "e\n");
            gp_close(gp);
        }
    }

    /* Twiddle vs direct: timing ratio */
    {
        int n = 1024, runs = 100;
        BenchResult direct = bench_fft_radix2(n, runs);
        double x_vals[2] = { 0.0, 1.0 };
        double y_vals[2] = { direct.avg_us, 0.0 };

        TwiddleTable *tt = twiddle_create(n);
        Complex *buf = (Complex *)malloc((size_t)n * sizeof(Complex));
        double sum_us = 0.0;
        for (int r = 0; r < runs; r++) {
            for (int i = 0; i < n; i++) { buf[i].re = 0.0; buf[i].im = 0.0; }
            double t0 = timer_usec();
            fft_with_twiddles(buf, n, tt);
            double t1 = timer_usec();
            sum_us += t1 - t0;
        }
        y_vals[1] = sum_us / runs;
        twiddle_destroy(tt);
        free(buf);

        gp_plot_1("ch29", "twiddle_speedup",
                  "FFT Timing: Direct sin/cos vs Twiddle Table (N=1024)",
                  "Method (0=Direct, 1=Twiddle)", "Average Time (us)",
                  x_vals, y_vals, 2, "boxes");
    }
}

/* ================================================================== */
/*  Main: generate all plots                                          */
/* ================================================================== */

int main(void)
{
    printf("╔══════════════════════════════════════════════════════════╗\n");
    printf("║  DSP Tutorial Suite — Generating All Gnuplot Plots     ║\n");
    printf("╚══════════════════════════════════════════════════════════╝\n\n");

    plot_ch01();
    plot_ch02();
    plot_ch03();
    plot_ch04();
    plot_ch05();
    plot_ch06();
    plot_ch07();
    plot_ch08();
    plot_ch09();
    plot_ch10();
    plot_ch11();
    plot_ch12();
    plot_ch13();
    plot_ch14();
    plot_ch15();
    plot_ch16();
    plot_ch18();
    plot_ch19();
    plot_ch17();
    plot_ch20();
    plot_ch21();
    plot_ch22();
    plot_ch23();
    plot_ch24();
    plot_ch25();
    plot_ch26();
    plot_ch27();
    plot_ch28();
    plot_ch29();
    plot_ch30();

    printf("\n  Done! All plots saved to plots/\n");
    printf("  View with: eog plots/ch01/impulse.png\n");
    printf("  or:        xdg-open plots/ch11/butterworth_orders.png\n\n");

    return 0;
}

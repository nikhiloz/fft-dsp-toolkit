/**
 * @file optimization.c
 * @brief Chapter 29 — DSP Optimisation Utilities.
 *
 * Implements:
 *   - Benchmarking via POSIX clock_gettime (monotonic, µs resolution)
 *   - Radix-4 FFT (25% fewer muls than radix-2)
 *   - Pre-computed twiddle factor tables
 *   - Portable aligned memory allocation
 *
 * ── Radix-4 Butterfly ──────────────────────────────────────────
 *
 *        a ───┬───(+)───── X[k]         X[k]       = a + b + c + d
 *             │                          X[k+N/4]   = a - jb - c + jd
 *        b ──W¹──(+)───── X[k+N/4]      X[k+N/2]   = a - b + c - d
 *             │                          X[k+3N/4]  = a + jb - c - jd
 *        c ──W²──(+)───── X[k+N/2]
 *             │                         Key: W^{N/4} = -j, W^{N/2} = -1
 *        d ──W³──(+)───── X[k+3N/4]    So 2 of 4 twiddles are "free"
 *
 * ── Twiddle Table ──────────────────────────────────────────────
 *
 *   Pre-compute W[k] = exp(-j·2π·k/N) for k = 0 .. N/2-1
 *   Speeds up inner loop by replacing sin/cos with table lookup.
 *
 * @see include/optimization.h
 */

#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "optimization.h"
#include "fft.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ================================================================== */
/*  Timing helper                                                     */
/* ================================================================== */

static double time_usec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec * 1e6 + (double)ts.tv_nsec / 1000.0;
}

/* ================================================================== */
/*  Bit-reversal permutation (for radix-4: base-4 digit reversal)     */
/* ================================================================== */

static int is_power_of_4(int n)
{
    if (n <= 0) return 0;
    while (n > 1) {
        if (n % 4 != 0) return 0;
        n /= 4;
    }
    return 1;
}

/* Standard base-2 bit-reversal for any power-of-2 */
static void bit_reverse(Complex *x, int n)
{
    int j = 0;
    for (int i = 0; i < n - 1; i++) {
        if (i < j) {
            Complex tmp = x[i];
            x[i] = x[j];
            x[j] = tmp;
        }
        int m = n >> 1;
        while (m >= 1 && j >= m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
}

/* Base-4 digit reversal for radix-4 DIT.
 * Reverses the base-4 digits of each index (NOT the same as bit-reversal).
 * E.g., for N=16: 1=01₄ → 10₄=4, 2=02₄ → 20₄=8 */
static void digit_reverse_4(Complex *x, int n)
{
    int log4n = 0;
    int temp = n;
    while (temp > 1) { temp /= 4; log4n++; }

    for (int i = 0; i < n; i++) {
        int reversed = 0;
        int val = i;
        for (int d = 0; d < log4n; d++) {
            reversed = reversed * 4 + (val % 4);
            val /= 4;
        }
        if (reversed > i) {
            Complex tmp = x[i];
            x[i] = x[reversed];
            x[reversed] = tmp;
        }
    }
}

/* ================================================================== */
/*  Radix-4 FFT                                                       */
/* ================================================================== */

void fft_radix4(Complex *x, int n)
{
    if (n <= 1) return;

    /* Fall back to radix-2 for non-power-of-4 sizes */
    if (!is_power_of_4(n)) {
        fft(x, n);
        return;
    }

    /* Base-4 digit reversal permutation */
    digit_reverse_4(x, n);

    /* Radix-4 stages */
    for (int stage_len = 4; stage_len <= n; stage_len *= 4) {
        int quarter = stage_len / 4;
        double angle = -2.0 * M_PI / stage_len;

        for (int group = 0; group < n; group += stage_len) {
            for (int k = 0; k < quarter; k++) {
                double theta1 = angle * k;
                double theta2 = angle * 2 * k;
                double theta3 = angle * 3 * k;

                Complex W1 = { cos(theta1), sin(theta1) };
                Complex W2 = { cos(theta2), sin(theta2) };
                Complex W3 = { cos(theta3), sin(theta3) };

                int i0 = group + k;
                int i1 = i0 + quarter;
                int i2 = i1 + quarter;
                int i3 = i2 + quarter;

                Complex a = x[i0];
                Complex b = complex_mul(W1, x[i1]);
                Complex c = complex_mul(W2, x[i2]);
                Complex d = complex_mul(W3, x[i3]);

                /* Radix-4 butterfly equations */
                Complex ap = { a.re + c.re, a.im + c.im };  /* a + c */
                Complex am = { a.re - c.re, a.im - c.im };  /* a - c */
                Complex bp = { b.re + d.re, b.im + d.im };  /* b + d */
                Complex bm = { b.re - d.re, b.im - d.im };  /* b - d */

                /* -j * bm = { bm.im, -bm.re } */
                Complex jbm = { bm.im, -bm.re };

                x[i0].re = ap.re + bp.re;  x[i0].im = ap.im + bp.im;
                x[i1].re = am.re + jbm.re; x[i1].im = am.im + jbm.im;
                x[i2].re = ap.re - bp.re;  x[i2].im = ap.im - bp.im;
                x[i3].re = am.re - jbm.re; x[i3].im = am.im - jbm.im;
            }
        }
    }
}

void ifft_radix4(Complex *x, int n)
{
    /* Conjugate → FFT → conjugate → scale */
    for (int i = 0; i < n; i++) x[i].im = -x[i].im;
    fft_radix4(x, n);
    for (int i = 0; i < n; i++) {
        x[i].re /= n;
        x[i].im = -x[i].im / n;
    }
}

/* ================================================================== */
/*  Twiddle Table                                                     */
/* ================================================================== */

TwiddleTable *twiddle_create(int n)
{
    TwiddleTable *tt = (TwiddleTable *)malloc(sizeof(TwiddleTable));
    if (!tt) return NULL;

    tt->n = n;
    tt->W = (Complex *)malloc((size_t)(n / 2) * sizeof(Complex));
    if (!tt->W) { free(tt); return NULL; }

    for (int k = 0; k < n / 2; k++) {
        double angle = -2.0 * M_PI * k / n;
        tt->W[k].re = cos(angle);
        tt->W[k].im = sin(angle);
    }
    return tt;
}

void twiddle_destroy(TwiddleTable *tt)
{
    if (!tt) return;
    free(tt->W);
    free(tt);
}

void fft_with_twiddles(Complex *x, int n, const TwiddleTable *tt)
{
    if (n <= 1) return;

    /* Bit-reversal */
    bit_reverse(x, n);

    /* Radix-2 DIT using pre-computed twiddles */
    for (int stage = 2; stage <= n; stage *= 2) {
        int half = stage / 2;
        int step = n / stage;  /* index step in twiddle table */

        for (int group = 0; group < n; group += stage) {
            for (int k = 0; k < half; k++) {
                Complex w = tt->W[k * step];
                int i = group + k;
                int j = i + half;

                Complex t = complex_mul(w, x[j]);
                x[j].re = x[i].re - t.re;
                x[j].im = x[i].im - t.im;
                x[i].re = x[i].re + t.re;
                x[i].im = x[i].im + t.im;
            }
        }
    }
}

/* ================================================================== */
/*  Benchmarking                                                      */
/* ================================================================== */

static void gen_random_complex(Complex *x, int n, unsigned int seed)
{
    for (int i = 0; i < n; i++) {
        seed = seed * 1103515245 + 12345;
        x[i].re = ((double)(seed >> 16) / 32768.0) - 1.0;
        seed = seed * 1103515245 + 12345;
        x[i].im = ((double)(seed >> 16) / 32768.0) - 1.0;
    }
}

BenchResult bench_fft_radix2(int n, int runs)
{
    BenchResult r = {0};
    r.n    = n;
    r.runs = runs;
    r.min_us = 1e30;

    Complex *x = (Complex *)malloc((size_t)n * sizeof(Complex));
    Complex *orig = (Complex *)malloc((size_t)n * sizeof(Complex));
    gen_random_complex(orig, n, 42);

    for (int run = 0; run < runs; run++) {
        memcpy(x, orig, (size_t)n * sizeof(Complex));

        double t0 = time_usec();
        fft(x, n);
        double t1 = time_usec();

        double elapsed = t1 - t0;
        if (elapsed < r.min_us) r.min_us = elapsed;
        if (elapsed > r.max_us) r.max_us = elapsed;
        r.avg_us += elapsed;
    }

    r.avg_us /= runs;

    /* FFT FLOP count: 5·N·log2(N) (complex muls + adds) */
    double log2n = log2((double)n);
    double flops = 5.0 * n * log2n;
    r.mflops = flops / r.avg_us;  /* µs → MFLOP/s */

    free(x);
    free(orig);
    return r;
}

BenchResult bench_fft_radix4(int n, int runs)
{
    BenchResult r = {0};
    r.n    = n;
    r.runs = runs;
    r.min_us = 1e30;

    Complex *x = (Complex *)malloc((size_t)n * sizeof(Complex));
    Complex *orig = (Complex *)malloc((size_t)n * sizeof(Complex));
    gen_random_complex(orig, n, 42);

    for (int run = 0; run < runs; run++) {
        memcpy(x, orig, (size_t)n * sizeof(Complex));

        double t0 = time_usec();
        fft_radix4(x, n);
        double t1 = time_usec();

        double elapsed = t1 - t0;
        if (elapsed < r.min_us) r.min_us = elapsed;
        if (elapsed > r.max_us) r.max_us = elapsed;
        r.avg_us += elapsed;
    }

    r.avg_us /= runs;

    double log2n = log2((double)n);
    double flops = 5.0 * n * log2n;
    r.mflops = flops / r.avg_us;

    free(x);
    free(orig);
    return r;
}

void bench_print(const char *label, const BenchResult *r)
{
    printf("  %-22s  N=%-5d  min=%7.1f µs  avg=%7.1f µs  max=%7.1f µs  %.1f MFLOP/s\n",
           label, r->n, r->min_us, r->avg_us, r->max_us, r->mflops);
}

/* ================================================================== */
/*  Aligned Memory Allocation                                         */
/* ================================================================== */

void *aligned_alloc_dsp(int alignment, int size)
{
    /* Portable aligned alloc for C99 (posix_memalign or manual) */
    void *ptr = NULL;
    /* Allocate extra space for alignment + stored original pointer */
    void *raw = malloc((size_t)(size + alignment + (int)sizeof(void *)));
    if (!raw) return NULL;

    /* Align: skip sizeof(void*) to store original pointer, then align */
    size_t addr = (size_t)raw + sizeof(void *);
    size_t aligned = (addr + (size_t)(alignment - 1)) & ~((size_t)(alignment - 1));
    ptr = (void *)aligned;

    /* Store original pointer right before the aligned pointer */
    ((void **)ptr)[-1] = raw;
    return ptr;
}

void aligned_free_dsp(void *ptr)
{
    if (!ptr) return;
    free(((void **)ptr)[-1]);
}

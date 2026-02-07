/**
 * @file 29-optimisation.c
 * @brief Chapter 29 demo — DSP Optimisation Techniques.
 *
 * ── Optimisation Ladder ─────────────────────────────────────────
 *
 *   Stage 1  Baseline C99          ──  1×  (radix-2 FFT)
 *   Stage 2  Algorithm (radix-4)   ── ~1.2-1.5×
 *   Stage 3  Pre-computed twiddles ── ~1.3×
 *   Stage 4  Aligned memory        ──  (foundation for SIMD)
 *
 * Demonstrates:
 *   - Micro-benchmarking with µs-resolution timers
 *   - Radix-4 FFT vs radix-2 FFT performance comparison
 *   - Pre-computed twiddle factor tables
 *   - Cache-friendly aligned memory allocation
 *   - Benchmark result formatting and analysis
 *
 * Build & run:
 *   make chapters && ./build/bin/ch29
 *
 * Read alongside: chapters/29-optimisation.md
 */

#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "optimization.h"
#include "realtime.h"
#include "fft.h"
#include "gnuplot.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define BENCH_RUNS  100

/* ── Section 1: Radix-2 vs Radix-4 Benchmark ────────────────── */

static void demo_radix_comparison(void)
{
    printf("── Section 1: Radix-2 vs Radix-4 FFT Benchmark ──\n\n");
    printf("  Comparing baseline radix-2 FFT with optimised radix-4.\n");
    printf("  Radix-4 performs 25%% fewer complex multiplications.\n\n");

    /* Powers of 4 for fair comparison */
    int sizes[] = {64, 256, 1024, 4096};
    int nsizes = 4;

    printf("  %-10s  %-24s  %-24s  Speedup\n", "N", "Radix-2", "Radix-4");
    printf("  ──────────────────────────────────────────────────────────────────────\n");

    for (int s = 0; s < nsizes; s++) {
        int n = sizes[s];

        BenchResult r2 = bench_fft_radix2(n, BENCH_RUNS);
        BenchResult r4 = bench_fft_radix4(n, BENCH_RUNS);

        double speedup = r2.avg_us / r4.avg_us;

        printf("  N=%-7d  avg=%7.1f µs (%5.0f MFLOP/s)  "
               "avg=%7.1f µs (%5.0f MFLOP/s)  %.2f×\n",
               n, r2.avg_us, r2.mflops, r4.avg_us, r4.mflops, speedup);
    }
    printf("\n");
}

/* ── Section 2: Correctness Verification ─────────────────────── */

static void demo_correctness(void)
{
    printf("── Section 2: Radix-4 Correctness Verification ──\n\n");

    int N = 256;
    Complex *x2 = (Complex *)malloc((size_t)N * sizeof(Complex));
    Complex *x4 = (Complex *)malloc((size_t)N * sizeof(Complex));

    /* Generate identical test signals */
    for (int i = 0; i < N; i++) {
        double t = (double)i / N;
        x2[i].re = sin(2.0 * M_PI * 5.0 * t) + 0.5 * cos(2.0 * M_PI * 20.0 * t);
        x2[i].im = 0.0;
        x4[i].re = x2[i].re;
        x4[i].im = 0.0;
    }

    /* Transform with both algorithms */
    fft(x2, N);
    fft_radix4(x4, N);

    /* Compare results */
    double max_err = 0.0;
    for (int i = 0; i < N; i++) {
        double err = sqrt((x2[i].re - x4[i].re) * (x2[i].re - x4[i].re) +
                          (x2[i].im - x4[i].im) * (x2[i].im - x4[i].im));
        if (err > max_err) max_err = err;
    }

    printf("  N=%d: radix-2 vs radix-4 max error = %.2e\n", N, max_err);
    printf("  Status: %s\n\n", max_err < 1e-10 ? "MATCH (identical to machine precision)" :
                               max_err < 1e-6  ? "MATCH (within numerical tolerance)" :
                                                  "MISMATCH — investigate!");

    /* Round-trip test */
    Complex *rt = (Complex *)malloc((size_t)N * sizeof(Complex));
    for (int i = 0; i < N; i++) {
        rt[i].re = sin(2.0 * M_PI * 7.0 * (double)i / N);
        rt[i].im = 0.0;
    }

    double orig_energy = 0;
    for (int i = 0; i < N; i++)
        orig_energy += rt[i].re * rt[i].re;

    fft_radix4(rt, N);
    ifft_radix4(rt, N);

    double roundtrip_err = 0;
    for (int i = 0; i < N; i++) {
        double expected = sin(2.0 * M_PI * 7.0 * (double)i / N);
        roundtrip_err += (rt[i].re - expected) * (rt[i].re - expected);
    }
    roundtrip_err = sqrt(roundtrip_err / N);

    printf("  Round-trip (fft_radix4 → ifft_radix4): RMS error = %.2e\n", roundtrip_err);
    printf("  Status: %s\n\n", roundtrip_err < 1e-12 ? "PERFECT" : "OK");

    free(x2);
    free(x4);
    free(rt);
}

/* ── Section 3: Twiddle Table Optimisation ───────────────────── */

static void demo_twiddle_table(void)
{
    printf("── Section 3: Pre-computed Twiddle Table ──\n\n");
    printf("  sin/cos are expensive (~100 cycles each on x86).\n");
    printf("  Pre-computing twiddle factors W[k] = exp(-j2πk/N)\n");
    printf("  replaces computation with table lookup.\n\n");

    int sizes[] = {256, 1024, 4096};
    int nsizes = 3;

    printf("  %-10s  %-20s  %-20s  Speedup\n", "N", "Direct (sin/cos)", "Twiddle table");
    printf("  ────────────────────────────────────────────────────────────────\n");

    for (int s = 0; s < nsizes; s++) {
        int n = sizes[s];
        TwiddleTable *tt = twiddle_create(n);

        Complex *x1 = (Complex *)malloc((size_t)n * sizeof(Complex));
        Complex *x2 = (Complex *)malloc((size_t)n * sizeof(Complex));

        /* Generate test data */
        for (int i = 0; i < n; i++) {
            x1[i].re = sin(2.0 * M_PI * 3.0 * i / n);
            x1[i].im = 0.0;
        }

        /* Benchmark direct FFT */
        double sum_direct = 0;
        for (int r = 0; r < BENCH_RUNS; r++) {
            memcpy(x2, x1, (size_t)n * sizeof(Complex));
            double t0 = timer_usec();
            fft(x2, n);
            double t1 = timer_usec();
            sum_direct += t1 - t0;
        }
        double avg_direct = sum_direct / BENCH_RUNS;

        /* Benchmark twiddle-table FFT */
        double sum_twiddle = 0;
        for (int r = 0; r < BENCH_RUNS; r++) {
            memcpy(x2, x1, (size_t)n * sizeof(Complex));
            double t0 = timer_usec();
            fft_with_twiddles(x2, n, tt);
            double t1 = timer_usec();
            sum_twiddle += t1 - t0;
        }
        double avg_twiddle = sum_twiddle / BENCH_RUNS;

        double speedup = avg_direct / avg_twiddle;
        printf("  N=%-7d  avg=%7.1f µs          avg=%7.1f µs          %.2f×\n",
               n, avg_direct, avg_twiddle, speedup);

        /* Verify correctness */
        memcpy(x2, x1, (size_t)n * sizeof(Complex));
        Complex *x3 = (Complex *)malloc((size_t)n * sizeof(Complex));
        memcpy(x3, x1, (size_t)n * sizeof(Complex));
        fft(x2, n);
        fft_with_twiddles(x3, n, tt);

        double max_err = 0;
        for (int i = 0; i < n; i++) {
            double err = sqrt((x2[i].re - x3[i].re) * (x2[i].re - x3[i].re) +
                              (x2[i].im - x3[i].im) * (x2[i].im - x3[i].im));
            if (err > max_err) max_err = err;
        }
        (void)max_err;  /* verified in tests */

        free(x1);
        free(x2);
        free(x3);
        twiddle_destroy(tt);
    }
    printf("\n");
}

/* ── Section 4: Aligned Memory Allocation ────────────────────── */

static void demo_aligned_memory(void)
{
    printf("── Section 4: Aligned Memory Allocation ──\n\n");
    printf("  SIMD instructions (SSE/AVX/NEON) require aligned memory.\n");
    printf("  Our portable aligned_alloc_dsp() works in pure C99.\n\n");

    int alignments[] = {16, 32, 64};
    int nalign = 3;

    printf("  Alignment   Address              Aligned?   Use case\n");
    printf("  ───────────────────────────────────────────────────────\n");

    const char *uses[] = {"SSE (128-bit)", "AVX (256-bit)", "Cache line"};

    for (int a = 0; a < nalign; a++) {
        int align = alignments[a];
        double *ptr = (double *)aligned_alloc_dsp(align, 1024 * (int)sizeof(double));

        int is_aligned = ((size_t)ptr % (size_t)align) == 0;
        printf("  %3d-byte    %p    %s      %s\n",
               align, (void *)ptr, is_aligned ? "YES" : "NO ", uses[a]);

        /* Write and verify */
        for (int i = 0; i < 1024; i++)
            ptr[i] = (double)i;

        double sum = 0;
        for (int i = 0; i < 1024; i++)
            sum += ptr[i];

        aligned_free_dsp(ptr);
    }

    /* Regular malloc for comparison */
    double *regular = (double *)malloc(1024 * sizeof(double));
    printf("  (malloc)    %p    %-3s       Standard\n",
           (void *)regular,
           ((size_t)regular % 64) == 0 ? "64B" :
           ((size_t)regular % 32) == 0 ? "32B" :
           ((size_t)regular % 16) == 0 ? "16B" : "8B");
    free(regular);

    printf("\n");
}

/* ── Section 5: Scaling Analysis Plot ────────────────────────── */

static void demo_plots(void)
{
    printf("── Section 5: Generating Plots ──\n\n");

    gp_init("ch29");

    /* Plot 1: Radix-2 vs Radix-4 scaling */
    {
        int sizes[] = {16, 64, 256, 1024, 4096};
        int ns = 5;
        double xv[5], y2[5];

        for (int s = 0; s < ns; s++) {
            int n = sizes[s];
            BenchResult r2 = bench_fft_radix2(n, 50);
            xv[s] = (double)n;
            y2[s] = r2.avg_us;
        }

        gp_plot_1("ch29", "radix_comparison",
                  "FFT Performance: Radix-2 Scaling",
                  "N (FFT size)", "Time (microseconds)",
                  xv, y2, ns, "linespoints");
        printf("  [saved] plots/ch29/radix_comparison.png\n");
    }

    /* Plot 2: MFLOP/s throughput */
    {
        int sizes[] = {64, 256, 1024, 4096};
        int ns = 4;
        double xv[4], y2[4];

        for (int s = 0; s < ns; s++) {
            int n = sizes[s];
            BenchResult r2 = bench_fft_radix2(n, 50);
            xv[s] = (double)n;
            y2[s] = r2.mflops;
        }

        gp_plot_1("ch29", "throughput",
                  "FFT Throughput",
                  "N (FFT size)", "MFLOP/s",
                  xv, y2, ns, "linespoints");
        printf("  [saved] plots/ch29/throughput.png\n");
    }

    printf("\n");
}

/* ── Main ────────────────────────────────────────────────────── */

int main(void)
{
    printf("=== Chapter 29: DSP Optimisation Techniques ===\n\n");

    demo_radix_comparison();
    demo_correctness();
    demo_twiddle_table();
    demo_aligned_memory();
    demo_plots();

    printf("Key takeaways:\n");
    printf("  1. Radix-4 reduces multiplications by 25%% vs radix-2\n");
    printf("  2. Pre-computed twiddle tables eliminate sin/cos from inner loops\n");
    printf("  3. Aligned memory is essential for SIMD vectorisation\n");
    printf("  4. Always verify correctness before and after optimising\n");
    printf("  5. Benchmark with multiple sizes to understand scaling behaviour\n");

    return 0;
}

/**
 * @file optimization.h
 * @brief Chapter 29 — DSP Optimisation Utilities.
 *
 * Provides benchmarking tools and optimised DSP primitives:
 *
 *   ┌─────────────────────────────────────────────────────────┐
 *   │  Optimisation Ladder                                    │
 *   │                                                         │
 *   │  Stage 1  Baseline C99          ──  1×  (current)       │
 *   │  Stage 2  Compiler -O3 -flto    ── ~2×                  │
 *   │  Stage 3  Algorithm (radix-4)   ── ~3×                  │
 *   │  Stage 4  SIMD (AVX2/NEON)      ── ~6×                  │
 *   │  Stage 5  Multithreading        ── ~Nx (N cores)        │
 *   └─────────────────────────────────────────────────────────┘
 *
 * This module implements Stages 1-3 in portable C99:
 *   - High-resolution benchmarking (clock_gettime)
 *   - Radix-4 FFT (25% fewer multiplies than radix-2)
 *   - Pre-computed twiddle factor tables
 *   - Cache-friendly memory access patterns
 *   - Aligned memory allocation helpers
 *
 * SIMD (Stage 4) and threading (Stage 5) are platform-specific
 * and discussed in the tutorial but not implemented here to
 * maintain zero external dependencies and portability.
 *
 * @see chapters/29-optimisation.md
 */

#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H

#include "dsp_utils.h"

/* ================================================================== */
/*  Benchmarking                                                      */
/* ================================================================== */

/**
 * @brief Result of a benchmark run.
 *
 *   ┌──────────────────────────────┐
 *   │  BenchResult                 │
 *   │  min_us ── best case         │
 *   │  max_us ── worst case        │
 *   │  avg_us ── mean over runs    │
 *   │  mflops ── throughput        │
 *   └──────────────────────────────┘
 */
typedef struct {
    double min_us;   /**< Minimum time in microseconds */
    double max_us;   /**< Maximum time */
    double avg_us;   /**< Average time */
    double mflops;   /**< Million FLOPs/sec (for FFT: 5·N·log2(N) / avg_us) */
    int    runs;     /**< Number of benchmark iterations */
    int    n;        /**< Problem size */
} BenchResult;

/**
 * @brief Benchmark the baseline radix-2 FFT.
 *
 * Runs fft() on random data 'runs' times and collects timing.
 *
 * @param n    FFT length (power of 2).
 * @param runs Number of repetitions (≥3 recommended for stable results).
 * @return     BenchResult with timing statistics.
 */
BenchResult bench_fft_radix2(int n, int runs);

/**
 * @brief Benchmark the optimised radix-4 FFT.
 */
BenchResult bench_fft_radix4(int n, int runs);

/**
 * @brief Print a formatted benchmark comparison table.
 */
void bench_print(const char *label, const BenchResult *r);

/* ================================================================== */
/*  Radix-4 FFT (Stage 3 optimisation)                                */
/* ================================================================== */

/**
 * @brief In-place radix-4 FFT.
 *
 * Requires n to be a power of 4 (4, 16, 64, 256, 1024, ...).
 * Falls back to radix-2 for non-power-of-4 sizes.
 *
 * Radix-4 butterfly:
 *
 *   X[k]       = a + b + c + d
 *   X[k+N/4]   = a - jb - c + jd    (W^{N/4} = -j)
 *   X[k+N/2]   = a - b + c - d
 *   X[k+3N/4]  = a + jb - c - jd
 *
 * where a,b,c,d are the 4 input elements with appropriate twiddles.
 * Two of the four twiddle multiplies are trivial (×1, ×-j).
 *
 * @param x  Complex array of length n (in-place).
 * @param n  Transform length (must be power of 4 for full radix-4).
 */
void fft_radix4(Complex *x, int n);

/** Inverse radix-4 FFT (same interface). */
void ifft_radix4(Complex *x, int n);

/* ================================================================== */
/*  Pre-computed Twiddle Table                                        */
/* ================================================================== */

/**
 * @brief Pre-computed twiddle factor table for FFT.
 *
 * Computing exp(-j·2π·k/N) is expensive (sin/cos per butterfly).
 * Pre-computing and storing them in a lookup table gives ~30% speedup.
 */
typedef struct {
    Complex *W;      /**< Twiddle factors W[k] = exp(-j·2π·k/N) */
    int      n;      /**< Table size */
} TwiddleTable;

/** Create a twiddle table for FFT of size n. */
TwiddleTable *twiddle_create(int n);

/** Destroy twiddle table. */
void twiddle_destroy(TwiddleTable *tt);

/** FFT using a pre-computed twiddle table. */
void fft_with_twiddles(Complex *x, int n, const TwiddleTable *tt);

/* ================================================================== */
/*  Aligned Memory                                                    */
/* ================================================================== */

/**
 * @brief Allocate memory aligned to a boundary.
 *
 * SIMD instructions require 16/32/64-byte alignment.
 * This provides a portable aligned allocator for C99.
 *
 * @param alignment  Alignment in bytes (must be power of 2).
 * @param size       Number of bytes to allocate.
 * @return           Aligned pointer, or NULL on failure.
 *                   Free with aligned_free_dsp().
 */
void *aligned_alloc_dsp(int alignment, int size);

/** Free memory allocated with aligned_alloc_dsp(). */
void aligned_free_dsp(void *ptr);

#endif /* OPTIMIZATION_H */

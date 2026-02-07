# Chapter 7 — Optimisation

Once the algorithms are correct, we make them fast. This chapter covers
the optimisation strategy for DSP code on modern CPUs.

> **Status:** This chapter describes the planned approach. The SIMD
> kernels and multithreading are planned for Phase 4 of the project.

---

## 7.1 The Five-Stage Approach

> **📊 Optimisation Roadmap** — [View full-size diagram →](../reference/diagrams/optimization_roadmap.png)

| Stage | Technique | Expected Speedup |
|-------|-----------|-----------------|
| 1. Baseline | Pure C99 scalar code | 1× |
| 2. Compiler | `-O3`, LTO, PGO | 1.5–2× |
| 3. Algorithm | Radix-4, cache-friendly layout | 2–3× |
| 4. SIMD | AVX2 (x86), NEON (ARM) | 4–8× |
| 5. Multithreading | OpenMP parallel stages | N× (N cores) |

**Rule:** Always optimise in this order. Algorithmic improvements
compound with SIMD, and SIMD compounds with threading.

## 7.2 Stage 1: The Baseline

Our current `fft()` in [`src/fft.c`](../src/fft.c) is stage 1:
correct, readable C99 with no tricks. This is the foundation that
all optimisations build upon.

Key hotspot: the butterfly inner loop (lines 101–115).
A 1024-point FFT runs this loop $512 \times 10 = 5120$ times,
each iteration doing one `complex_mul` (4 multiplies, 2 adds).

## 7.3 Stage 2: Compiler Optimisations

Already in our Makefile:

```makefile
RELEASE_FLAGS = -O3 -DNDEBUG
```

Additional flags to explore:
- **`-march=native`** — Use all CPU instructions available
- **`-ffast-math`** — Allow FP reordering (breaks strict IEEE 754)
- **`-flto`** — Link-Time Optimisation (inline across .c files)
- **`-fprofile-generate` / `-fprofile-use`** — Profile-Guided Optimisation

## 7.4 Stage 3: Algorithmic Improvements

### Radix-4 FFT

Instead of 2-element butterflies, use 4-element butterflies.
Reduces the number of complex multiplies by 25%:

- Radix-2: $\frac{N}{2} \log_2 N$ multiplies
- Radix-4: $\frac{3N}{8} \log_2 N$ multiplies

### Split-Radix

Combines radix-2 for odd stages and radix-4 for even stages.
Achieves the lowest known multiply count for power-of-2 FFTs.

### Cache-Friendly Memory Layout

Modern CPUs have L1 caches of 32–64 KB. A 1024-point complex array
is $1024 \times 16 = 16$ KB — fits in L1. But for larger FFTs:

- Process data in blocks that fit in L1
- Use **stride-1 access patterns** (butterfly partners should be
  adjacent in memory)
- Consider **out-of-order** (Stockham) FFT for better locality

## 7.5 Stage 4: SIMD Vectorisation

SIMD (Single Instruction, Multiple Data) processes 4–8 values
simultaneously.

### Complex Multiply in AVX2

Standard complex multiply: 4 multiplies + 2 adds per element.
With AVX2, we process 4 complex numbers at once:

```c
/* Pseudocode — AVX2 intrinsics */
__m256d ar = _mm256_load_pd(/* 4 real parts of a */);
__m256d ai = _mm256_load_pd(/* 4 imag parts of a */);
__m256d br = _mm256_load_pd(/* 4 real parts of b */);
__m256d bi = _mm256_load_pd(/* 4 imag parts of b */);

/* (a.re * b.re - a.im * b.im) */
__m256d cr = _mm256_sub_pd(_mm256_mul_pd(ar, br),
                           _mm256_mul_pd(ai, bi));
/* (a.re * b.im + a.im * b.re) */
__m256d ci = _mm256_add_pd(_mm256_mul_pd(ar, bi),
                           _mm256_mul_pd(ai, br));
```

**Key requirement:** Data must be aligned to 32-byte boundaries.
Use `aligned_alloc(32, size)` instead of `malloc`.

### ARM NEON

For Raspberry Pi and embedded targets:

```c
float32x4_t ar = vld1q_f32(/* 4 real parts */);
/* Similar pattern with vmulq_f32, vsubq_f32, vaddq_f32 */
```

NEON works on 4 floats (128-bit) vs AVX2's 4 doubles (256-bit).

## 7.6 Stage 5: Multithreading

For very large FFTs or batch processing:

```c
#pragma omp parallel for
for (int group = 0; group < n; group += stage_size) {
    /* each group is independent — perfect for parallelism */
}
```

The outer two loops of `fft()` have independent groups within each
stage, so they parallelise cleanly.

**Caution:** For small FFTs (<4096 points), threading overhead exceeds
the benefit. Only parallelise when the data is large enough.

## 7.7 Benchmarking Targets

> **📊 Performance Benchmarks** — [View full-size diagram →](../reference/diagrams/benchmarks.png)

| Implementation | 1024-pt Latency | Goal |
|----------------|---------------:|------|
| Baseline C99 | ~12 ms | Stage 1 (current) |
| + Compiler opts | ~6 ms | Stage 2 |
| + Radix-4 | ~3 ms | Stage 3 |
| + AVX2 SIMD | ~0.8 ms | Stage 4 |
| + OpenMP (4 cores) | ~0.3 ms | Stage 5 |
| FFTW3 (reference) | ~3.5 ms | Comparison |

## 7.8 Exercises

1. **Benchmark:** Time our current `fft()` for N=256, 1024, 4096,
   16384. Plot N vs. time. Does it follow $O(N \log N)$?

2. **Compiler test:** Build with `-O0`, `-O2`, `-O3`,
   `-O3 -march=native -ffast-math`. Compare FFT timing for each.

3. **Code exercise:** Implement the Radix-4 butterfly. The key
   equation uses $W^0, W^{N/4}, W^{N/2}, W^{3N/4}$ — note that
   $W^{N/2} = -1$ and $W^{N/4} = -i$, so two of the four multiplies
   are free.

4. **Thinking question:** Why does `-ffast-math` help FFT performance?
   What precision do we sacrifice?

---

**Previous:** [Chapter 06 — Real-Time Streaming](06-real-time-streaming.md)
| **Next:** [Chapter 08 — Putting It Together →](08-putting-it-together.md)

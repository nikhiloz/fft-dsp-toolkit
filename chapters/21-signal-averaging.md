# Chapter 21 — Signal Averaging & Noise Reduction

## Overview

Signal averaging exploits the fact that noise is random while the signal of
interest is deterministic (or repeatable). By combining multiple observations,
the signal adds coherently while noise partially cancels, improving the
**signal-to-noise ratio (SNR)**.

## Key Concepts

### Coherent (Synchronous) Averaging

Given K time-aligned trials of the same signal corrupted by independent noise:

    y[n] = (1/K) · Σ_{k=1}^{K} (s[n] + w_k[n])
         = s[n] + (1/K) · Σ w_k[n]

Since noise is zero-mean and uncorrelated:
- Signal power unchanged: P_s.
- Noise power reduced: P_w / K.
- **SNR improvement**: 10·log₁₀(K) dB.

Requirements:
- Signal must be periodic or repeatable (e.g., evoked potentials, radar pulses).
- Trials must be time-aligned (synchronised trigger).

### Exponential Moving Average (EMA)

    y[n] = α·x[n] + (1−α)·y[n−1]

- α near 0: heavy smoothing, large group delay.
- α near 1: tracks input closely, minimal smoothing.
- Cutoff frequency: fc ≈ α·fs / (2π) for small α.
- Memory: infinite (IIR), but effective window ≈ 2/α samples.

### Moving Average

    y[n] = (1/M) · Σ_{k=0}^{M-1} x[n-k]

- FIR filter with equal weights.
- Frequency response: sinc-like, first null at fs/M.
- Good for random noise; poor for periodic interference.

### Median Filter

    y[n] = median{x[n-M/2], ..., x[n], ..., x[n+M/2]}

- Nonlinear filter (not representable as convolution).
- Excellent for impulse (salt-and-pepper) noise.
- Preserves edges (step transitions) that moving average blurs.
- Computational cost: O(M log M) per sample.

## SNR Concepts

    SNR = 10·log₁₀(P_signal / P_noise) dB

| Method             | SNR Gain           | Best For                |
|-------------------|--------------------|-------------------------|
| Coherent avg (K)  | 10·log₁₀(K) dB   | Repeatable signals      |
| Moving average (M)| ~10·log₁₀(M) dB  | White noise             |
| EMA (α)           | ~10·log₁₀(2/α) dB| Online/streaming data   |
| Median (M)        | Robust vs impulse  | Outliers, spike removal |

## API Reference

```c
#include "averaging.h"

void coherent_average(const double **trials, int K, int n, double *out);
void ema_filter(const double *x, int n, double alpha, double *y);
void moving_average(const double *x, int n, int M, double *y);
void median_filter(const double *x, int n, int M, double *y);
void compute_snr_improvement(const double *x_noisy, const double *x_clean,
                             const double *x_avg, int n,
                             double *snr_before, double *snr_after);
```

## Running the Demo

```bash
make chapters
./build/bin/ch21
```

Plots are written to `plots/ch21/`.

## Cross-References

- [Chapter 10: FIR Filters](10-fir-filter-design.md) — moving average as FIR
- [Chapter 11: IIR Filters](11-iir-filter-design.md) — EMA as first-order IIR
- [Chapter 15: Correlation](15-correlation-xcorr.md) — SNR estimation

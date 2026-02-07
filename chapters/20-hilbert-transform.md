# Chapter 20 — Hilbert Transform & Analytic Signals

## Overview

The **Hilbert transform** shifts all frequency components of a signal by −90°.
When combined with the original signal, it forms the **analytic signal**, which
has powerful applications in envelope detection, instantaneous frequency
estimation, and single-sideband modulation.

## Key Concepts

### The Hilbert Transform

For a real signal x(t), the Hilbert transform x̂(t) is:

    x̂(t) = H{x(t)} = (1/πt) * x(t)

In the frequency domain:
- Positive frequencies: multiplied by −j (phase shift −90°)
- Negative frequencies: multiplied by +j (phase shift +90°)
- DC and Nyquist: unchanged

### Analytic Signal

    z(t) = x(t) + j·x̂(t)

Properties:
- **One-sided spectrum**: only positive frequencies survive.
- **Envelope**: |z(t)| = A(t) (instantaneous amplitude).
- **Phase**: ∠z(t) = φ(t) (instantaneous phase).
- **Instantaneous frequency**: f_i(t) = (1/2π) · dφ/dt.

### FIR Hilbert Filter (Type III)

A causal approximation of the ideal Hilbert transform:

    h[n] = { 2/(πn)  for n odd,   window applied
              0        for n even               }

Type III FIR: antisymmetric coefficients, odd length.

### FFT Method

1. Compute X(k) = FFT{x(n)}.
2. Zero out negative frequencies: X(k) = 0 for k = N/2+1 ... N−1.
3. Double positive frequencies: X(k) ×= 2 for k = 1 ... N/2−1.
4. Z(k) = IFFT{X(k)} → analytic signal.

This method is exact (no transient), but requires block processing.

### Applications

| Application             | Uses                          |
|------------------------|-------------------------------|
| AM Demodulation        | Envelope = |z(t)|             |
| FM Demodulation        | f_i = d∠z(t)/dt / (2π)       |
| SSB Modulation         | USB = Re{z(t)·e^{j2πfct}}    |
| Radar pulse analysis   | Instantaneous amplitude/phase |
| Speech processing      | Pitch tracking, formants      |

## API Reference

```c
#include "hilbert.h"

void    hilbert_design(double *h, int taps);
void    analytic_signal(const double *x, int n, Complex *z, int taps);
void    analytic_signal_fft(const double *x, int n, Complex *z);
void    envelope(const double *x, int n, double *env, int taps);
void    inst_frequency(const double *x, int n, double *freq, int taps);
```

When `taps = 0`, functions use the FFT method automatically.

## Running the Demo

```bash
make chapters
./build/bin/ch20
```

Plots are written to `plots/ch20/`.

## Cross-References

- [Chapter 10: FIR Filters](10-fir-filter-design.md) — Type III FIR structure
- [Chapter 08: FFT](08-fft.md) — FFT-based analytic signal
- [Chapter 03: Complex Numbers](03-complex-numbers.md) — complex magnitude and phase

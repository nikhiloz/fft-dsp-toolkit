# Chapter 22 — Advanced FIR Design (Parks-McClellan / Remez)

## Overview

The **Parks-McClellan algorithm** (also called the **Remez exchange**) designs
FIR filters that are **optimal in the Chebyshev (minimax) sense**: the maximum
error between the desired and actual frequency response is minimised. This
produces **equiripple** filters — the passband and stopband ripples are equal
in magnitude and alternate in sign.

## Key Concepts

### Why Equiripple?

Window methods (Hamming, Kaiser, etc.) produce filters where the maximum error
is in the transition band. The Remez algorithm instead:

1. Distributes the error **equally** across the passband and stopband.
2. Achieves the **narrowest transition band** for a given filter order.
3. Provides explicit control over passband/stopband ripple via weights.

### The Chebyshev Alternation Theorem

An optimal length-N filter has at least N/2 + 2 points where the weighted
error alternates between +δ and −δ (equiripple extremals). This property
guarantees the design is minimax optimal.

### Remez Exchange Algorithm

1. **Initialise**: Choose (N/2 + 2) trial extremal frequencies.
2. **Interpolate**: Find the Lagrange polynomial through these points.
3. **Compute delta**: Calculate the equiripple deviation δ.
4. **Evaluate error**: Compute error on a dense frequency grid.
5. **Update extremals**: Replace trial frequencies with actual error peaks.
6. **Repeat** until δ converges.

### Design Parameters

| Parameter        | Effect                                   |
|-----------------|------------------------------------------|
| Filter order N  | More taps → narrower transition           |
| Passband weight | Higher → tighter passband ripple          |
| Stopband weight | Higher → deeper stopband rejection        |
| Transition Δf   | Narrower → needs more taps               |

### Comparison: Window vs Remez

| Property          | Window Method      | Remez (Equiripple)   |
|------------------|--------------------|----------------------|
| Transition width | Wider              | Narrower (optimal)   |
| Ripple pattern   | Decays from edge   | Equal throughout     |
| Design control   | Window type + order| Weights + band edges |
| Computational    | Fast (direct)      | Iterative            |
| Optimality       | Not guaranteed     | Minimax optimal      |

## API Reference

```c
#include "remez.h"

typedef struct {
    double low;      /* lower band edge (normalised 0–0.5) */
    double high;     /* upper band edge */
    double desired;  /* desired gain (1.0 = passband, 0.0 = stopband) */
    double weight;   /* error weight for this band */
} RemezBand;

int remez_fir(double *h, int taps, const RemezBand *bands,
              int n_bands, int max_iter);

int remez_lowpass(double *h, int taps,
                  double fpass, double fstop,
                  double wpass, double wstop);

int remez_bandpass(double *h, int taps,
                   double fstop1, double fpass1,
                   double fpass2, double fstop2);
```

Returns 0 on success, −1 if convergence fails.

## Running the Demo

```bash
make chapters
./build/bin/ch22
```

Plots are written to `plots/ch22/`.

## Cross-References

- [Chapter 10: FIR Filters](10-fir-filter-design.md) — Window method comparison
- [Chapter 09: Window Functions](09-window-functions.md) — Window shapes
- [Chapter 17: Multirate DSP](17-multirate-dsp.md) — Remez filters for anti-alias

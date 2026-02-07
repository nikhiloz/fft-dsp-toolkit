# Chapter 23 — Adaptive Filters: LMS, NLMS, RLS

## Overview

**Adaptive filters** are filters whose coefficients are adjusted automatically
to minimise an error signal. Unlike fixed FIR/IIR designs, adaptive filters
**learn** the optimal coefficients by processing data sample-by-sample,
making them ideal for environments where signal statistics change over time.

## Key Concepts

### The Adaptive Filtering Problem

```
x[n] → ┌────────┐ → y[n] ─→ ⊕ → e[n]
        │  w[k]  │           ↑
        └────────┘         d[n] (desired signal)
             ↑
        e[n] updates w[k]
```

The adaptive filter produces output y[n] = Σ w[k]·x[n−k]. The error
e[n] = d[n] − y[n] drives the weight update. The goal: minimise E{e²[n]}.

### LMS (Least Mean Squares)

The simplest and most widely-used adaptive algorithm:

**Update rule**: w[k] ← w[k] + μ · e[n] · x[n−k]

| Property | Value |
|----------|-------|
| Complexity | O(L) per sample |
| Step size μ | 0 < μ < 2/λ_max |
| Convergence | Slow, depends on eigenvalue spread |
| Stability | Guaranteed for small μ |

The step size μ controls the speed/accuracy trade-off:
- Large μ → fast convergence, higher misadjustment
- Small μ → slow convergence, lower steady-state error

### NLMS (Normalised LMS)

Normalises the step size by the input power to handle varying signal levels:

**Update rule**: w[k] ← w[k] + μ/(‖x‖² + ε) · e[n] · x[n−k]

The regularisation ε prevents division by zero. NLMS converges faster than
LMS when the input signal has varying power, and the step size is less
sensitive to the input statistics.

### RLS (Recursive Least Squares)

Minimises the **total weighted least-squares error** rather than the
instantaneous squared error:

$$J[n] = \sum_{i=0}^{n} \lambda^{n-i} e^2[i]$$

The forgetting factor λ (typically 0.95–1.0) exponentially down-weights
older data, allowing the filter to track non-stationary environments.

**Key equations**:
- Gain vector: k[n] = P[n−1]·x[n] / (λ + x[n]ᵀ·P[n−1]·x[n])
- Error: e[n] = d[n] − wᵀ[n−1]·x[n] (a priori)
- Weight update: w[n] = w[n−1] + k[n]·e[n]
- Matrix update: P[n] = (P[n−1] − k[n]·x[n]ᵀ·P[n−1]) / λ

| Property | LMS | NLMS | RLS |
|----------|-----|------|-----|
| Complexity | O(L) | O(L) | O(L²) |
| Convergence | Slow | Medium | Fast |
| Tracking | Poor | Good | Excellent |
| Memory | L | L | L² |

## Applications

### System Identification

Estimate an unknown system's impulse response by driving it with a known
input and comparing outputs:
- **x[n]**: input to both unknown system and adaptive filter
- **d[n]**: output of unknown system
- After convergence: w[k] ≈ h[k] (unknown system coefficients)

### Noise Cancellation

Remove correlated noise from a signal:
- **Primary input**: signal + noise passed through unknown path
- **Reference input**: noise reference (correlated with the noise)
- **Error output**: clean signal estimate

### Echo Cancellation

In telecommunications, remove the echo of the far-end speaker:
- The adaptive filter models the echo path (speaker → microphone)
- Subtracts the estimated echo from the microphone signal

## Implementation Notes

- Initialise weights to zero; P₀ = δ·I for RLS (δ = 100 typical)
- Circular buffer for delay line avoids data shifting
- For NLMS, ε = 10⁻⁶ is typical; prevents instability with silence
- RLS converges in ~2L samples (L = filter length)

## Demo

Run the Chapter 23 demo:
```bash
make chapters && ./build/bin/ch23
```

### Generated Plots

![LMS Learning Curve](../plots/ch23/lms_learning_curve.png)

![NLMS Noise Cancellation](../plots/ch23/nlms_noise_cancel.png)

![Convergence Comparison](../plots/ch23/convergence_compare.png)

## Further Reading

- Haykin, *Adaptive Filter Theory* (5th ed.), Chapters 5 (LMS), 9 (RLS)
- Widrow & Stearns, *Adaptive Signal Processing*
- Sayed, *Fundamentals of Adaptive Filtering*

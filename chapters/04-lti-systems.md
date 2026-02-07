# Chapter 4: LTI Systems & Discrete Convolution

> **Prerequisites:** [Ch 01 — Signals & Sequences](01-signals-and-sequences.md), [Ch 02 — Sampling & Aliasing](02-sampling-and-aliasing.md)
> **Library module:** [`convolution.h`](../include/convolution.h) — convolution, correlation, energy/power
> **Runnable demo:** [`04-lti-systems.c`](04-lti-systems.c) — `make chapters && ./build/bin/ch04s`

---

## Why This Chapter Matters

Nearly all useful DSP operations — filtering, echo, reverb, equalization,
noise removal — are **LTI** (Linear Time-Invariant) systems. If you understand
two properties (linearity + time-invariance) and one operation (convolution),
you understand the core mechanism of every filter in this tutorial.

> **Textbook references:**
> - Oppenheim & Willsky, *Signals and Systems*, Ch 1–2
> - Oppenheim & Schafer, *Discrete-Time Signal Processing*, Ch 2
> - Proakis & Manolakis, *Digital Signal Processing*, Ch 2
> - Smith, *Scientist & Engineer's Guide to DSP*, Ch 5–7
> - Lyons, *Understanding DSP*, Ch 5

---

## 4.1 What Is a System?

A **discrete-time system** is anything that takes an input sequence `x[n]`
and produces an output sequence `y[n]`:

$$
y[n] = T\{x[n]\}
$$

Examples: a moving average filter, an echo effect, a gain block.

---

## 4.2 Linearity

A system `T` is **linear** if it obeys both:

1. **Additivity:** $T\{x_1[n] + x_2[n]\} = T\{x_1[n]\} + T\{x_2[n]\}$
2. **Homogeneity (scaling):** $T\{a \cdot x[n]\} = a \cdot T\{x[n]\}$

Combined: $T\{a \cdot x_1[n] + b \cdot x_2[n]\} = a \cdot T\{x_1[n]\} + b \cdot T\{x_2[n]\}$

This is the **superposition principle**. It means the output to a sum of
inputs is the sum of the individual outputs. Enormously powerful — we can
analyze complex signals component by component.

**Non-linear example:** $y[n] = x[n]^2$ — scaling by `a` gives $a^2 x[n]^2$,
not $a \cdot x[n]^2$.

---

## 4.3 Time-Invariance

A system is **time-invariant** if a delay in the input produces exactly the
same delay in the output:

$$
\text{If } y[n] = T\{x[n]\}, \text{ then } T\{x[n-n_0]\} = y[n-n_0].
$$

The system's behaviour doesn't change over time. A filter that works today
will work the same way tomorrow.

**Non-TI example:** $y[n] = n \cdot x[n]$ — the output depends on the
absolute time index.

---

## 4.4 The Impulse Response

The **impulse response** `h[n]` of an LTI system is its output when the
input is a unit impulse `δ[n]`:

$$
h[n] = T\{\delta[n]\}
$$

This single sequence completely characterises the system. Every possible
output can be derived from `h[n]`:

$$
y[n] = \sum_{k=-\infty}^{\infty} x[k] \cdot h[n-k] = (x * h)[n]
$$

This is the **convolution sum**. It's the most important equation in LTI
system theory.

---

## 4.5 Convolution: The Core Operation

### The Convolution Sum

$$
y[n] = \sum_{k=0}^{M-1} h[k] \cdot x[n-k]
$$

(For FIR systems with impulse response of length `M`.)

### How to Think About It

For each output sample `y[n]`:
1. Flip the impulse response `h` to get `h[-k]`
2. Shift it by `n` to get `h[n-k]`
3. Multiply point-by-point with `x[k]`
4. Sum all the products

### In C

```c
int convolve(const double *x, int x_len,
             const double *h, int h_len,
             double *y);
```

The output has length `x_len + h_len - 1`. This "growth" makes sense —
the filter's transient response adds samples at the beginning and end.

### Properties of Convolution

| Property | Meaning |
|----------|---------|
| Commutative | $x * h = h * x$ |
| Associative | $(x * h_1) * h_2 = x * (h_1 * h_2)$ — cascading filters |
| Distributive | $x * (h_1 + h_2) = x * h_1 + x * h_2$ — parallel filters |
| Identity | $x * \delta = x$ — convolving with impulse gives original |

---

## 4.6 BIBO Stability

A system is **Bounded-Input Bounded-Output (BIBO) stable** if every bounded
input produces a bounded output. For LTI systems, this happens if and only if:

$$
\sum_{n=-\infty}^{\infty} |h[n]| < \infty
$$

(The impulse response is **absolutely summable**.)

All FIR filters are BIBO stable (finite sums are always bounded). IIR filters
must be checked — their impulse response is infinitely long, and stability
depends on pole locations (covered in Ch 05: Z-Transform).

---

## 4.7 Causality

A system is **causal** if the output at time `n` depends only on current and
past inputs:

$$
h[n] = 0, \quad n < 0
$$

All real-time systems must be causal — you can't use future samples that
haven't arrived yet. For offline processing (batch), non-causal systems are
fine (e.g., zero-phase filtering).

---

## 4.8 FIR vs. IIR

| Property | FIR (Finite Impulse Response) | IIR (Infinite Impulse Response) |
|----------|------|------|
| h[n] length | Finite (M samples) | Infinite (decaying) |
| Stability | Always BIBO stable | Must verify (pole locations) |
| Linear phase | Possible (symmetric h) | Generally not |
| Efficiency | Needs more taps | Fewer coefficients for same spec |

This chapter implements everything as FIR (direct convolution). IIR design
and difference equations come in [Ch 11: IIR Filter Design](11-iir-design.md).

---

## 4.9 Cross-Correlation

Correlation measures similarity between two signals at various time shifts:

$$
r_{xy}[\ell] = \sum_{n} x[n] \cdot y[n + \ell]
$$

It's related to convolution: $r_{xy}[\ell] = x[-\ell] * y[\ell]$ (convolution
with one signal time-reversed).

**Auto-correlation** ($r_{xx}$) measures how a signal resembles a delayed
version of itself. It's central to periodicity detection, pitch estimation,
and spectral analysis.

---

## 4.10 Signal Energy and Power

| Quantity | Formula | Meaning |
|----------|---------|---------|
| Energy | $E = \sum |x[n]|^2$ | Total energy (finite-length signals) |
| Power | $P = \frac{1}{N} \sum |x[n]|^2$ | Average power per sample |

These connect to Parseval's theorem (Ch 07: DFT) — energy in time equals
energy in frequency.

---

## Exercises

1. **Convolution by hand:** Convolve `x = {1, 2, 3}` with `h = {1, 1, 1}`.
   Verify with the demo. (Answer: `{1, 3, 6, 5, 3}`.)

2. **Identity property:** Convolve any signal with `δ[n] = {1}`. Verify the
   output is identical to the input.

3. **Commutativity test:** Convolve x with h, then h with x. Are the results
   identical? (Use the demo or write code to confirm.)

---

## What's Next?

- [Ch 05: Z-Transform](05-z-transform.md) — the algebraic framework that turns convolution into multiplication.
- [Ch 06: Frequency Response, Poles & Zeros](06-frequency-response.md) — what convolution looks like in the frequency domain.

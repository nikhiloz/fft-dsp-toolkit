# Chapter 5: The Z-Transform

> **Prerequisites:** [Ch 01 — Signals & Sequences](01-signals-and-sequences.md), [Ch 03 — Complex Numbers](01-complex-numbers.md), [Ch 04 — LTI Systems](04-lti-systems.md)
> **Library module:** [`dsp_utils.h`](../include/dsp_utils.h) (complex arithmetic for evaluation)
> **Runnable demo:** [`05-z-transform.c`](05-z-transform.c) — `make chapters && ./build/bin/ch05s`

---

## Why This Chapter Matters

The Z-transform does for discrete-time signals what the Laplace transform
does for continuous-time: it converts **convolution into multiplication** and
turns difference equations into polynomials.

With the Z-transform, you can:
- Analyse filter stability by inspecting **poles**
- Cascade filters by multiplying transfer functions
- Design filters by placing poles and zeros strategically
- Understand why IIR filters can be unstable while FIR filters never are

> **Textbook references:**
> - Oppenheim & Willsky, *Signals and Systems*, Ch 10
> - Oppenheim & Schafer, *Discrete-Time Signal Processing*, Ch 3
> - Proakis & Manolakis, *Digital Signal Processing*, Ch 3
> - Lyons, *Understanding DSP*, Ch 6–7

---

## 5.1 Definition

The **Z-transform** of a discrete-time signal `x[n]` is:

$$
X(z) = \sum_{n=-\infty}^{\infty} x[n] \, z^{-n}
$$

where $z$ is a complex variable. This is a power series in $z^{-1}$.

For causal signals (x[n] = 0 for n < 0):

$$
X(z) = \sum_{n=0}^{\infty} x[n] \, z^{-n} = x[0] + x[1] z^{-1} + x[2] z^{-2} + \cdots
$$

---

## 5.2 Common Z-Transform Pairs

| Signal x[n] | X(z) | ROC |
|-------------|------|-----|
| δ[n] | 1 | All z |
| u[n] (unit step) | $\frac{1}{1 - z^{-1}}$ | \|z\| > 1 |
| $a^n u[n]$ | $\frac{1}{1 - a z^{-1}}$ | \|z\| > \|a\| |
| $n a^n u[n]$ | $\frac{a z^{-1}}{(1 - a z^{-1})^2}$ | \|z\| > \|a\| |
| FIR: {b₀, b₁, ..., bₘ} | $b_0 + b_1 z^{-1} + \cdots + b_M z^{-M}$ | All z ≠ 0 |

**Key insight:** An FIR filter's Z-transform is just a polynomial in $z^{-1}$.
Its coefficients are the filter taps.

---

## 5.3 Region of Convergence (ROC)

The Z-transform sum only converges for certain values of `z`. The set of `z`
values where it converges is the **ROC**.

Rules:
- The ROC is always an annular region: $r_1 < |z| < r_2$
- For **causal** signals: ROC is the exterior of a circle
- For **anti-causal** signals: ROC is the interior of a circle
- For **finite-length** signals: ROC is the entire z-plane (except possibly z=0 or z=∞)

**Why it matters:** Two different signals can have the same algebraic X(z)
but different ROCs — the ROC disambiguates them.

---

## 5.4 Properties of the Z-Transform

| Property | Time Domain | Z-Domain |
|----------|-------------|----------|
| Linearity | $a\,x_1[n] + b\,x_2[n]$ | $a\,X_1(z) + b\,X_2(z)$ |
| Time delay | $x[n - k]$ | $z^{-k} X(z)$ |
| Time advance | $x[n + k]$ | $z^{k} X(z)$ |
| Convolution | $(x * h)[n]$ | $X(z) \cdot H(z)$ |
| Multiplication by n | $n \, x[n]$ | $-z \frac{dX(z)}{dz}$ |
| Initial value | $x[0]$ | $\lim_{z \to \infty} X(z)$ |

The **convolution property** is the most important: it means filtering in
the time domain = multiplication in the Z-domain.

---

## 5.5 The Transfer Function H(z)

For an LTI system with impulse response `h[n]`:

$$
H(z) = \frac{Y(z)}{X(z)} = \sum_{n=0}^{\infty} h[n] z^{-n}
$$

For a system described by a **difference equation**:

$$
y[n] = \sum_{k=0}^{M} b_k \, x[n-k] - \sum_{k=1}^{N} a_k \, y[n-k]
$$

the transfer function is:

$$
H(z) = \frac{B(z)}{A(z)} = \frac{b_0 + b_1 z^{-1} + \cdots + b_M z^{-M}}{1 + a_1 z^{-1} + \cdots + a_N z^{-N}}
$$

- **FIR filters:** A(z) = 1 (no denominator), so $H(z) = B(z)$
- **IIR filters:** A(z) ≠ 1 (feedback terms), rational function

---

## 5.6 Poles and Zeros

Factor H(z) into:

$$
H(z) = \frac{b_0}{a_0} \cdot \frac{(1 - q_1 z^{-1})(1 - q_2 z^{-1}) \cdots}{(1 - p_1 z^{-1})(1 - p_2 z^{-1}) \cdots}
$$

- **Zeros** ($q_i$): values of z where H(z) = 0 — output is blocked
- **Poles** ($p_i$): values of z where H(z) → ∞ — system resonates

### Stability from Poles

> A causal LTI system is **BIBO stable** if and only if all poles lie
> **inside** the unit circle: $|p_i| < 1$ for all $i$.

This is the Z-domain equivalent of the absolute summability condition
from Ch 04. Poles inside the unit circle → decaying impulse response → stable.

---

## 5.7 Frequency Response from H(z)

The **frequency response** is H(z) evaluated on the unit circle:

$$
H(e^{j\omega}) = H(z)\big|_{z = e^{j\omega}} = |H(e^{j\omega})| \, e^{j\angle H(e^{j\omega})}
$$

where $\omega = 2\pi f / f_s$ is the normalised frequency.

In C, we evaluate this by substituting $z = e^{j\omega}$ into the polynomial:

```c
Complex z = complex_from_polar(1.0, omega);
Complex Hz = evaluate_transfer_function(b, M, a, N, z);
double magnitude = complex_mag(Hz);
double phase = complex_phase(Hz);
```

This gives us the magnitude and phase response at every frequency — the
core of filter analysis.

---

## 5.8 Inverse Z-Transform

Going from X(z) back to x[n] is the **inverse Z-transform**. Methods:

1. **Partial fraction expansion** — decompose into simple terms, look up each
2. **Power series expansion** — long division of B(z)/A(z) gives coefficients
3. **Residue theorem** — contour integration (theoretical, not for C code)

For our C implementations, we typically don't compute inverse Z-transforms
in closed form. Instead, we compute the impulse response by running the
difference equation with δ[n] as input.

---

## 5.9 Connecting It All Together

The Z-transform is the unifying framework:

```
Time domain:          y[n] = x[n] * h[n]     (convolution — Ch 04)
Z-domain:             Y(z) = X(z) · H(z)      (multiplication)
Frequency response:   H(e^jω)                 (evaluate on unit circle)
Stability:            |poles| < 1              (check pole locations)
Filter design:        place poles & zeros      (Ch 06, Ch 11)
```

Everything connects through H(z).

---

## Exercises

1. **Compute by hand:** Find the Z-transform of ``h[n] = {1, -2, 1}``.
   (Answer: $H(z) = 1 - 2z^{-1} + z^{-2}$. It's an FIR — just a polynomial.)

2. **Find the poles:** For $H(z) = \frac{1}{1 - 0.9z^{-1}}$, where is the
   pole? Is the system stable? (Answer: pole at z = 0.9, inside unit circle → stable.)

3. **Frequency response:** Run the demo and observe how moving a pole closer
   to the unit circle sharpens the resonance peak.

---

## What's Next?

- [Ch 06: Frequency Response, Poles & Zeros](06-frequency-response.md) — visual analysis of magnitude/phase plots.
- [Ch 07: The DFT](07-dft-theory.md) — the Z-transform evaluated at N equally-spaced points on the unit circle.

# Chapter 12: Filter Structures

> **Prerequisites:** [Ch 06 — Frequency Response](06-frequency-response.md), [Ch 11 — IIR Filter Design](11-iir-filter-design.md)
> **Library module:** [`iir.h`](../include/iir.h) (`biquad_process_df1`, `biquad_process_df2t`, `sos_process_block`)
> **Runnable demo:** [`12-filter-structures.c`](12-filter-structures.c) — `make release && ./build/bin/ch12`

---

## Why This Chapter Matters

The *same* transfer function $H(z)$ can be implemented in many different ways.
The choice of **structure** affects:

- **Numerical precision** (especially in fixed-point)
- **Memory usage** (number of state variables)
- **Coefficient sensitivity** (how much the response changes when coefficients are rounded)
- **Computational cost** (multiplies and adds per sample)

For a 2nd-order section (biquad), the differences are small in double-precision
floating-point. But for high-order filters or fixed-point embedded systems,
the structure matters enormously.

> **Textbook references:**
> - Oppenheim & Schafer, *Discrete-Time Signal Processing*, Ch 6
> - Proakis & Manolakis, *Digital Signal Processing*, Ch 9
> - Lyons, *Understanding DSP*, Ch 6–7

---

## 12.1 The Biquad (Second-Order Section)

Every IIR filter can be factored into a cascade of **2nd-order sections**
(biquads). The general biquad transfer function is:

$$
H(z) = \frac{b_0 + b_1 z^{-1} + b_2 z^{-2}}{1 + a_1 z^{-1} + a_2 z^{-2}}
$$

This has 5 coefficients ($b_0, b_1, b_2, a_1, a_2$) since $a_0 = 1$ by
convention (normalised form).

---

## 12.2 Direct Form I

The most straightforward implementation — directly evaluates the difference equation:

$$
y[n] = b_0 x[n] + b_1 x[n-1] + b_2 x[n-2] - a_1 y[n-1] - a_2 y[n-2]
$$

**Signal flow graph:**

```
x[n] ──►[×b0]──┬──►(+)──────────────────────► y[n]
            z⁻¹ │       ▲               ▲
x[n-1]──►[×b1]─┘       │               │
            z⁻¹ │   [×(-a1)]◄──y[n-1]  │
x[n-2]──►[×b2]─┘       │               │
                    [×(-a2)]◄──y[n-2]───┘
```

| Property | Value |
|----------|-------|
| State variables | 4 ($x_1, x_2, y_1, y_2$) |
| Multiplies/sample | 5 |
| Adds/sample | 4 |
| Memory | More than DF2 |
| Overflow risk | Lower (no intermediate accumulation) |

**Advantage:** In fixed-point, DF1 cannot overflow internally if the output
is within range, because intermediate sums don't exceed the output range.

---

## 12.3 Direct Form II Transposed (DF2T)

A mathematically equivalent structure that uses only **2 state variables**:

$$
\begin{aligned}
y[n] &= b_0 x[n] + s_1 \\
s_1 &= b_1 x[n] - a_1 y[n] + s_2 \\
s_2 &= b_2 x[n] - a_2 y[n]
\end{aligned}
$$

**Signal flow graph:**

```
x[n]──►[×b0]──►(+)──────────────► y[n]
   │            ▲             │
   │           s1          [×(-a1)]
   │            ▲             │
   ├──►[×b1]──►(+)◄──────────┘
   │            ▲             │
   │           s2          [×(-a2)]
   │            ▲             │
   └──►[×b2]──►(+)◄──────────┘
```

| Property | Value |
|----------|-------|
| State variables | 2 ($s_1, s_2$) |
| Multiplies/sample | 5 |
| Adds/sample | 4 |
| Memory | Minimal |
| Numerical behaviour | Best for floating-point |

**Advantage:** DF2T is the recommended form for floating-point implementations.
The state variables store *differences* rather than accumulated sums, which
reduces rounding error accumulation for narrow-band (high-Q) filters.

---

## 12.4 Cascaded Biquads (Second-Order Sections)

For an Nth-order IIR filter, factor into $\lceil N/2 \rceil$ biquad sections:

$$
H(z) = G \cdot \prod_{k=0}^{K-1} \frac{b_{0k} + b_{1k} z^{-1} + b_{2k} z^{-2}}{1 + a_{1k} z^{-1} + a_{2k} z^{-2}}
$$

```
x[n] ──► [SOS 0] ──► [SOS 1] ──► ··· ──► [SOS K−1] ──► ×G ──► y[n]
```

### Why Cascade?

A single 6th-order section has coefficients like $a_6 = r^6$ where $r$ is
the pole radius. If $r = 0.99$, then $a_6 = 0.941$, but the individual
coefficients interact in complex ways. Rounding any one coefficient can shift
all poles simultaneously.

In cascade form, each biquad handles at most **one conjugate pole pair**.
Rounding a coefficient in section $k$ only affects that section's poles,
not the others.

| Structure | Sensitivity | Memory | Recommendation |
|-----------|------------|--------|----------------|
| Single high-order | Very high | Minimal | ✗ Never use for order > 2 |
| Cascaded SOS | Low | K × 2 states | ✓ Always use |
| Parallel SOS | Low | K × 2 states | ✓ Special cases |

---

## 12.5 Section Ordering in a Cascade

The biquad sections can be ordered in different ways. The choice affects
intermediate signal dynamic range:

1. **Poles closest to unit circle first** — tends to amplify, then attenuate
2. **Poles farthest from unit circle first** — attenuates early, lower dynamic range
3. **Interleave high-Q and low-Q sections** — balanced dynamic range (often best)

For floating-point with sufficient headroom, ordering rarely matters. For
fixed-point, careful ordering prevents overflow.

---

## 12.6 Comparison Summary

```
                    ┌─────────────────────────────────────────┐
                    │         WHICH STRUCTURE TO USE?          │
                    ├─────────────────────────────────────────┤
                    │                                         │
                    │  Floating-point?                        │
                    │    ├── Yes → DF2T + SOS cascade         │
                    │    └── No  → Fixed-point?               │
                    │              ├── Yes → DF1 + SOS cascade│
                    │              └── Other → DF1 + SOS      │
                    │                                         │
                    │  Single biquad (2nd order)?              │
                    │    └── Any form is fine                  │
                    │                                         │
                    │  Order > 2?                             │
                    │    └── ALWAYS use SOS cascade            │
                    │                                         │
                    └─────────────────────────────────────────┘
```

---

## 12.7 Exercises

1. **DF1 vs DF2T memory:** For a 10th-order filter implemented as 5 cascaded
   biquads, how many state variables does DF1 require total? DF2T?

2. **Polynomial expansion:** Take a 4th-order Butterworth designed with
   `butterworth_lowpass()`. Expand the SOS cascade to a single numerator and
   denominator polynomial. Compare the impulse responses of both forms for
   100 samples. How large is the maximum difference?

3. **Coefficient quantisation:** Quantise all biquad coefficients to 16-bit
   fixed-point (divide by $2^{-15}$, round, multiply back). Process a sine
   wave through both the original and quantised versions. How much does the
   output differ?

---

## 12.8 Summary

| Structure | States | Best For | Avoid When |
|-----------|--------|----------|------------|
| Direct Form I | 4 per biquad | Fixed-point | — |
| DF2 Transposed | 2 per biquad | Floating-point | — |
| Single high-order | $2N$ | Never (order > 2) | Always |
| Cascaded SOS | $2K$ (DF2T) | All practical IIR | — |

**Key rule:** *Never implement an IIR filter of order > 2 as a single section.
Always use cascaded second-order sections.*

**Next:** [Ch 13 — Spectral Analysis](13-spectral-analysis.md) — apply filters
and window functions to real-world spectral estimation.

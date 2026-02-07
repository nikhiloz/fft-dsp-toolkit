# Chapter 18 — Fixed-Point Arithmetic & Quantisation

| Navigation | |
|---|---|
| **Previous:** [Chapter 15 — Correlation](15-correlation.md) | **Next:** [Chapter 19 — Advanced FFT](19-advanced-fft.md) |
| **Code:** [18-fixed-point.c](18-fixed-point.c) | **API:** [fixed_point.h](../include/fixed_point.h) |

---

## Why Fixed-Point?

Most desktop CPUs have hardware floating-point units (FPUs), but many embedded
processors — ARM Cortex-M0, TI C2000, older DSP chips — either lack an FPU or
have limited FP throughput.  Fixed-point arithmetic performs all DSP operations
using **integer instructions** with an _implied binary point_.

```
  Floating-point (IEEE 754):        Fixed-point (Q1.15):
  ┌────┬──────────┬───────────┐     ┌───┬──────────────────┐
  │sign│ exponent │ mantissa  │     │ S │  15 frac bits     │
  └────┴──────────┴───────────┘     └───┴──────────────────┘
   32 bits, ~7 sig digits            16 bits, ~4.5 sig digits
   Huge dynamic range                Fixed range [-1, +1)
   Needs FPU for speed               Integer ALU only
```

**Key trade-off**: fixed-point is faster and more deterministic on integer-only
hardware, but has limited dynamic range and introduces quantisation noise.

---

## Q-Format Notation

The **Q-format** specifies how many bits represent the fractional part:

| Format | Bits | Sign | Integer | Fraction | Range | Resolution |
|--------|------|------|---------|----------|-------|------------|
| Q1.15  | 16   | 1    | 0       | 15       | [-1, +1) | 3.05×10⁻⁵ |
| Q1.31  | 32   | 1    | 0       | 31       | [-1, +1) | 4.66×10⁻¹⁰ |
| Q8.8   | 16   | 0    | 8       | 8        | [0, 256) | 3.91×10⁻³ |

For DSP, **Q1.15** (abbreviated Q15) is the most common: it represents numbers
in the range [-1.0, +1.0), which is the natural range for normalised audio
samples and filter coefficients.

### Conversion

```
double → Q15:    q15_value = (int16_t)(x × 32768)
Q15 → double:    double_value = q15_value / 32768.0
```

---

## Saturating Arithmetic

In two's complement integer arithmetic, **overflow wraps around**:
`32767 + 1 = -32768`.  This creates harsh distortion.  **Saturation** clamps
the result to the maximum representable value instead:

```
  Without saturation (wrap):       With saturation (clamp):

  +1.0 ┬                           +1.0 ┬━━━━━━━━━━
       │\                                │         ╱
       │ \  ← wrap to -1.0!             │        ╱
  -1.0 ┼──\──────                  -1.0 ┼───────╱
         input →                          input →
```

Our `q15_add()`, `q15_sub()`, `q15_mul()` all use saturation.

---

## Quantisation Noise

Converting a continuous (or high-resolution) signal to fixed-point introduces
**quantisation error**:

```
  x(t) ───► Quantise ───► x̂[n]
                 │
                 ▼
            e[n] = x[n] - x̂[n]    (noise floor)
```

For uniformly distributed quantisation error with **B** fractional bits:

$$\text{SQNR} \approx 6.02 \cdot B + 1.76 \;\text{dB}$$

| Format | Bits (B) | Theoretical SQNR |
|--------|----------|-------------------|
| Q15    | 15       | ~92 dB            |
| Q31    | 31       | ~188 dB           |
| 8-bit  | 7        | ~44 dB            |

---

## Fixed-Point FIR Filter

A standard FIR filter in floating-point:

```c
y[n] = Σ h[k] × x[n-k]    for k = 0..taps-1
```

In Q15 fixed-point, the multiply produces a Q30 result (15+15 frac bits).
We accumulate in a 32-bit register, then shift right by 15 to get Q15 output:

```
  Accumulator (32-bit / Q30):
  ┌─────────────────────────────────┐
  │ SSign │      30 fractional bits │
  └─────────────────────────────────┘
      acc += (int32_t)h[k] * (int32_t)x[n-k]
                                         ▼
  Output = saturate_q15( acc >> 15 )     
  ┌───┬──────────────────┐
  │ S │  15 frac bits     │   Q15
  └───┴──────────────────┘
```

This approach maintains precision across multiply-accumulate operations
and is standard practice on DSP processors like the TI C5000/C6000 series.

---

## Demo Walkthrough

### Demo 1: Conversion Round-Trip
Shows double → Q15 → double for various values.  Small values like 0.001
suffer visible quantisation error (~3×10⁻⁵ step size).

### Demo 2: Saturating Arithmetic
Demonstrates that 0.75 + 0.5 = 0.999969 (saturated), not 1.25.

### Demo 3: SQNR Measurement
Quantises a 440 Hz sine to both Q15 and Q31, measures actual SQNR.
Compare with the theoretical 6.02B + 1.76 formula.

### Demo 4: Float vs Q15 FIR
Applies the same 31-tap lowpass filter using both `fir_filter()` (double)
and `fir_filter_q15()` (Q15), then computes the SQNR between outputs.

![Q15 Quantisation Error](../plots/ch18/quantisation_error.png)
![Q15 Signal vs Recovered](../plots/ch18/q15_quantisation.png)
![Q15 Error](../plots/ch18/q15_error.png)
![FIR Float vs Q15](../plots/ch18/fir_float_vs_q15.png)

### Demo 5: Saturation Visual
Scales a sine wave at different gains (0.5×, 0.9×, 1.5×, 2.0×) through
Q15 conversion.  At gains > 1.0, clipping is clearly visible.

![Saturation at Different Gains](../plots/ch18/saturation.png)

---

## Key Takeaways

1. **Q15** represents [-1, +1) in 16 bits — perfect for audio and filter coefficients
2. **Saturation** prevents wrap-around artefacts at the cost of harmonic distortion
3. **SQNR ≈ 6.02B + 1.76 dB** — each extra bit gives ~6 dB of dynamic range
4. **32-bit accumulator** is essential for FIR filters in Q15 to avoid precision loss
5. Always **benchmark float vs fixed** — on modern CPUs, float may be faster

---

| Navigation | |
|---|---|
| **Previous:** [Chapter 15 — Correlation](15-correlation.md) | **Next:** [Chapter 19 — Advanced FFT](19-advanced-fft.md) |

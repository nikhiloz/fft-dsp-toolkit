# Chapter 11: IIR Filter Design

> **Prerequisites:** [Ch 05 — Z-Transform](05-z-transform.md), [Ch 06 — Frequency Response](06-frequency-response.md), [Ch 10 — FIR Filters](10-digital-filters.md)
> **Library module:** [`iir.h`](../include/iir.h) (`butterworth_lowpass`, `chebyshev1_lowpass`, etc.)
> **Runnable demo:** [`11-iir-filter-design.c`](11-iir-filter-design.c) — `make release && ./build/bin/ch11`

---

## Why This Chapter Matters

FIR filters (Ch 10) are always stable and have linear phase, but they need
**many taps** for a sharp transition. IIR filters achieve the same selectivity
with far fewer coefficients — often 4 vs 50+ multiplies per sample.

The trade-off: IIR filters have non-linear phase and can be unstable if
designed incorrectly. This chapter teaches you to design them properly.

> **Textbook references:**
> - Oppenheim & Schafer, *Discrete-Time Signal Processing*, Ch 8
> - Proakis & Manolakis, *Digital Signal Processing*, Ch 8
> - Lyons, *Understanding DSP*, Ch 7
> - Smith, *Scientist & Engineer's Guide to DSP*, Ch 19–21

---

## 11.1 The Design Method: Analog Prototype → Bilinear Transform

IIR filter design follows a well-established recipe:

1. **Choose an analog prototype** (Butterworth, Chebyshev, Elliptic)
2. **Pre-warp** the digital cutoff frequency to compensate for bilinear transform distortion
3. **Apply the bilinear transform** to convert analog poles/zeros to digital
4. **Factor into biquad sections** for numerical stability

```
┌─────────────┐     ┌──────────┐     ┌─────────────┐     ┌──────────┐
│   Analog     │────►│ Pre-warp │────►│  Bilinear   │────►│ Biquad   │
│  Prototype   │     │  cutoff  │     │  Transform  │     │ Cascade  │
│  (s-domain)  │     │          │     │  s → z      │     │ (SOS)    │
└─────────────┘     └──────────┘     └─────────────┘     └──────────┘
```

---

## 11.2 The Bilinear Transform

The bilinear transform maps the s-plane to the z-plane:

$$
s = \frac{2}{T} \cdot \frac{z - 1}{z + 1}
\quad \Leftrightarrow \quad
z = \frac{1 + sT/2}{1 - sT/2}
$$

With $T = 1$ (normalised):

$$
s = 2 \cdot \frac{z - 1}{z + 1}
$$

**Key properties:**
- Maps the $j\Omega$ axis → unit circle (frequency axis preserved)
- Maps left half-plane → inside unit circle (stability preserved)
- Maps right half-plane → outside unit circle

**Frequency warping:**

$$
\Omega_a = 2 \tan\left(\frac{\omega_d}{2}\right)
$$

Near DC, the mapping is approximately linear ($\Omega_a \approx \omega_d$).
Near Nyquist, the entire positive analog axis is compressed onto $[0, \pi]$,
causing severe non-linear distortion.

**Pre-warping** compensates by designing the analog filter at $\Omega_a$
instead of $\omega_d$, so the −3 dB point falls exactly at the desired
digital frequency.

---

## 11.3 Butterworth Filters

The Butterworth filter has the **maximally flat** magnitude response:

$$
|H_a(j\Omega)|^2 = \frac{1}{1 + (\Omega/\Omega_c)^{2N}}
$$

Properties:
- **No passband ripple** — monotonically decreasing
- **−3 dB at cutoff** — by definition
- **−20N dB/decade** rolloff in the stopband

### Analog Prototype Poles

The Nth-order Butterworth poles lie on a circle of radius $\Omega_c$:

$$
s_k = \Omega_c \cdot e^{j\pi(2k + N + 1)/(2N)}, \quad k = 0, 1, \ldots, N-1
$$

Only the left-half-plane poles (with $\text{Re}(s) < 0$) are used.

### Design Example (4th-order, cutoff = 0.2)

```
Step 1: Pre-warp  →  Ωa = 2·tan(π·0.2) = 1.453
Step 2: Poles at θ₀ = 5π/8, θ₁ = 7π/8 on circle of radius 1.453
Step 3: Bilinear transform each pole pair → z-domain biquad
Step 4: Cascade 2 biquad sections
```

---

## 11.4 Chebyshev Type I Filters

Chebyshev filters allow **passband ripple** in exchange for a steeper
transition band:

$$
|H_a(j\Omega)|^2 = \frac{1}{1 + \varepsilon^2 T_N^2(\Omega/\Omega_c)}
$$

where $T_N$ is the Nth Chebyshev polynomial and
$\varepsilon = \sqrt{10^{R_{\text{dB}}/10} - 1}$.

### Poles on an Ellipse

Unlike Butterworth (circle), Chebyshev poles lie on an **ellipse**:

$$
\text{Re}(s_k) = -\Omega_c \sinh(a) \sin(\theta_k), \quad
\text{Im}(s_k) = \Omega_c \cosh(a) \cos(\theta_k)
$$

where $a = \frac{1}{N} \text{arcsinh}(1/\varepsilon)$.

| Parameter | Butterworth | Chebyshev I |
|-----------|-------------|-------------|
| Passband | Maximally flat | Equiripple (±R dB) |
| Transition | Gradual | Steeper (same order) |
| Stopband | Monotonic | Monotonic |
| Poles | Circle | Ellipse |

---

## 11.5 Choosing the Right Filter

| Criterion | Choose | Why |
|-----------|--------|-----|
| No passband distortion | Butterworth | Maximally flat |
| Sharpest rolloff for given order | Elliptic | Best transition/order ratio |
| Sharp rolloff, moderate ripple OK | Chebyshev I | Good compromise |
| Phase matters more than magnitude | FIR | Linear phase guaranteed |

---

## 11.6 Implementation in C

Our library implements the full design chain:

```c
SOSCascade sos;

/* Butterworth lowpass: order 4, cutoff 0.2 (normalised) */
butterworth_lowpass(4, 0.2, &sos);

/* Process a signal block */
sos_process_block(&sos, input, output, n_samples);
```

The output `sos` contains biquad coefficients ready for real-time processing.
See [Ch 12 — Filter Structures](12-filter-structures.md) for the different
ways to implement each biquad section.

---

## 11.7 Exercises

1. **Order selection:** For a Butterworth lowpass with cutoff 0.1 and −40 dB
   at 0.2, what minimum order is needed? Verify with the library.

2. **Chebyshev trade-off:** Design a 4th-order Chebyshev I at cutoff 0.15
   with 0.5 dB, 1 dB, and 3 dB ripple. Compare the stopband attenuation at
   0.3.

3. **Highpass design:** Design a 6th-order Butterworth highpass at cutoff 0.3.
   Filter a signal containing 0.1·fs and 0.4·fs components. Verify the low
   component is removed.

4. **Pre-warping effect:** Design a 4th-order Butterworth at cutoff 0.45
   (near Nyquist) with and without pre-warping. How far is the −3 dB point
   from the target without pre-warping?

---

## 11.8 Summary

| Concept | Key Point |
|---------|-----------|
| Bilinear transform | Maps analog → digital, preserves stability |
| Pre-warping | $\Omega_a = 2\tan(\omega_d/2)$ compensates non-linear mapping |
| Butterworth | Maximally flat, −20N dB/decade rolloff |
| Chebyshev I | Equiripple passband, steeper rolloff for same order |
| SOS cascade | Factor into biquads for numerical robustness |

**Next:** [Ch 12 — Filter Structures](12-filter-structures.md) — how to
implement biquads efficiently (Direct Form I/II, Transposed, Cascaded).

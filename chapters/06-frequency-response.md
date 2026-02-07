# Chapter 6: Frequency Response, Poles & Zeros

> **Prerequisites:** [Ch 03 — Complex Numbers](03-complex-numbers.md), [Ch 05 — Z-Transform](05-z-transform.md)
> **Library module:** [`iir.h`](../include/iir.h) (`freq_response`, `group_delay_at`)
> **Runnable demo:** [`06-frequency-response.c`](06-frequency-response.c) — `make release && ./build/bin/ch06`

---

## Why This Chapter Matters

Every filter — FIR or IIR — has a **frequency response** that tells you exactly
what it does to each frequency in your signal. Understanding this is the
single most important skill in practical DSP:

- You can **predict** whether a filter passes or blocks a frequency
- You can **design** filters by manipulating poles and zeros
- You can **verify** implementations by comparing measured vs theoretical response
- You can **debug** audio glitches and spectral artefacts

> **Textbook references:**
> - Oppenheim & Schafer, *Discrete-Time Signal Processing*, Ch 2–3, 5
> - Proakis & Manolakis, *Digital Signal Processing*, Ch 3–4
> - Lyons, *Understanding DSP*, Ch 5–7

---

## 6.1 Definition: Frequency Response

For an LTI system with transfer function $H(z)$, the **frequency response**
is obtained by evaluating on the unit circle:

$$
H(e^{j\omega}) = H(z)\bigg|_{z = e^{j\omega}} = \frac{B(e^{j\omega})}{A(e^{j\omega})}
$$

where:

$$
B(e^{j\omega}) = \sum_{k=0}^{M} b_k \, e^{-jk\omega}, \quad
A(e^{j\omega}) = \sum_{k=0}^{N} a_k \, e^{-jk\omega}
$$

The frequency response is complex-valued. We decompose it:

- **Magnitude response:** $|H(e^{j\omega})|$ — gain at each frequency
- **Phase response:** $\angle H(e^{j\omega})$ — phase shift at each frequency
- **Group delay:** $\tau(\omega) = -\frac{d\phi(\omega)}{d\omega}$ — delay in samples

---

## 6.2 Magnitude Response

The magnitude $|H(e^{j\omega})|$ tells you how much the filter amplifies or
attenuates each frequency:

| Value | Meaning |
|-------|---------|
| $|H| = 1$ (0 dB) | Frequency passed unchanged |
| $|H| > 1$ (> 0 dB) | Frequency amplified |
| $|H| < 1$ (< 0 dB) | Frequency attenuated |
| $|H| = 0$ ($-\infty$ dB) | Frequency completely blocked |

In decibels: $|H|_{\text{dB}} = 20 \log_{10} |H(e^{j\omega})|$

**Example:** The 2-point average $H(z) = \frac{1 + z^{-1}}{2}$ has:

$$
|H(e^{j\omega})| = \left|\cos\left(\frac{\omega}{2}\right)\right|
$$

This equals 1 at DC ($\omega = 0$) and 0 at Nyquist ($\omega = \pi$) — a lowpass response.

---

## 6.3 Phase Response & Group Delay

The phase $\phi(\omega) = \angle H(e^{j\omega})$ determines how much each
frequency component is delayed.

**Linear phase** means $\phi(\omega) = -\alpha \omega$ for some constant $\alpha$.
This gives **constant group delay** $\tau = \alpha$ — all frequencies are delayed
equally, preserving the waveform shape.

| Filter Type | Phase | Group Delay | Waveform |
|-------------|-------|-------------|----------|
| FIR (symmetric) | Linear | Constant $\tau = (N-1)/2$ | Preserved |
| IIR | Non-linear | Varies with $\omega$ | Distorted |

**Why it matters:** In audio, non-linear phase causes **dispersion** — different
frequencies arrive at different times, smearing transients. FIR filters avoid
this problem.

---

## 6.4 Poles, Zeros, and Their Effect

The transfer function $H(z)$ can be factored as:

$$
H(z) = G \cdot \frac{(z - z_1)(z - z_2) \cdots (z - z_M)}{(z - p_1)(z - p_2) \cdots (z - p_N)}
$$

where $z_k$ are **zeros** and $p_k$ are **poles**.

### Rules of Thumb

| Element | Effect on Frequency Response |
|---------|-----------------------------|
| **Zero on unit circle** | Complete null (−∞ dB) at that frequency |
| **Zero near unit circle** | Deep notch at that frequency |
| **Pole on unit circle** | Infinite gain (unstable!) |
| **Pole near unit circle** | Sharp peak at that frequency |
| **Pole/zero distance from circle** | Width of peak/notch |

### Geometric Interpretation

At frequency $\omega$, the point on the unit circle is $z = e^{j\omega}$.
The magnitude response is:

$$
|H(e^{j\omega})| = |G| \cdot \frac{\prod |e^{j\omega} - z_k|}{\prod |e^{j\omega} - p_k|}
$$

Each factor $|e^{j\omega} - z_k|$ is the **distance** from the evaluation point
to the zero $z_k$ in the z-plane. Zeros close to the unit circle make their
distance factor small → magnitude dips. Poles close to the unit circle make
the denominator small → magnitude peaks.

---

## 6.5 Stability from Poles

**BIBO Stability Criterion:**

> A causal LTI system is BIBO (Bounded-Input Bounded-Output) stable
> if and only if **all poles lie strictly inside the unit circle**: $|p_k| < 1$.

| Pole location | System behaviour |
|---------------|-----------------|
| $|p| < 1$ | Stable — impulse response decays |
| $|p| = 1$ | Marginally stable — impulse response persists |
| $|p| > 1$ | Unstable — impulse response grows without bound |

FIR filters have all their poles at $z = 0$ (inside the unit circle),
so **FIR filters are always stable**. IIR filters have poles elsewhere,
so stability must be verified during design.

---

## 6.6 All-Pass Filters

An all-pass filter has $|H(e^{j\omega})| = 1$ for all $\omega$ but a
non-trivial phase response. The simplest example:

$$
H(z) = \frac{a + z^{-1}}{1 + a z^{-1}}
$$

The zero at $z = -a$ and pole at $z = -a$ are **reciprocals** of each other
with respect to the unit circle ($z_{\text{zero}} = 1/z_{\text{pole}}^*$).
This geometry forces unity magnitude.

**Uses:**
- Phase equalization (correct group delay distortion)
- Building lattice filter structures (see Ch 12)
- Constructing complementary filter pairs

---

## 6.7 Exercises

1. **Notch filter design:** Place a conjugate zero pair on the unit circle at
   $\omega = \pi/4$ and conjugate poles at the same angle but radius 0.9.
   Compute and plot the magnitude response. What happens when you move the
   poles closer to the unit circle?

2. **Phase comparison:** Compare the phase response of a 31-tap FIR lowpass
   (windowed sinc) with a 4th-order Butterworth IIR lowpass at the same cutoff.
   Which has flatter group delay?

3. **Stability experiment:** Start with a stable 2nd-order IIR filter
   (pole radius 0.9). Gradually increase the pole radius past 1.0 and observe
   the impulse response. At what radius does the system become unstable?

---

## 6.8 Summary

| Concept | Key Formula | Notes |
|---------|------------|-------|
| Frequency response | $H(e^{j\omega}) = B(e^{j\omega})/A(e^{j\omega})$ | Evaluate $H(z)$ on unit circle |
| Magnitude (dB) | $20\log_{10}|H|$ | 0 dB = unity gain |
| Group delay | $\tau = -d\phi/d\omega$ | Constant for linear-phase FIR |
| Stability | All $|p_k| < 1$ | FIR always stable; IIR needs checking |
| All-pass | $|H| = 1$, zero = $1/p^*$ | Modifies phase only |

**Next:** [Ch 10 — FIR Filter Design](10-digital-filters.md) (design filters to spec),
then [Ch 11 — IIR Filter Design](11-iir-filter-design.md) (Butterworth, Chebyshev).

# Chapter 19 — Advanced FFT: Goertzel, DTMF & Sliding DFT

| Navigation | |
|---|---|
| **Previous:** [Chapter 18 — Fixed-Point](18-fixed-point.md) | **Next:** [Chapter 16 — Overlap-Add/Save](16-overlap-add-save.md) |
| **Code:** [19-advanced-fft.c](19-advanced-fft.c) | **API:** [advanced_fft.h](../include/advanced_fft.h) |

---

## Beyond Radix-2: When You Don't Need All N Bins

The standard FFT computes **all N** frequency bins in O(N log N).  But many
real-world applications only need **one or a few** specific frequencies:

- **DTMF decoding**: 8 frequencies
- **Vibration monitoring**: 1-3 resonant frequencies
- **Pitch detection**: 1 fundamental
- **Power quality**: 50/60 Hz + harmonics

For these cases, the **Goertzel algorithm** is more efficient.

```
  Algorithm Selection:
  ┌─────────────────────────────────────────────────────────────┐
  │  How many bins do you need?                                │
  │                                                             │
  │    All N bins ──────────► FFT          O(N log N)          │
  │                                                             │
  │    K bins, K < log₂(N) ─► Goertzel    O(K × N)            │
  │                                                             │
  │    1 bin, streaming ────► Sliding DFT  O(1) per sample     │
  └─────────────────────────────────────────────────────────────┘
```

---

## The Goertzel Algorithm

Goertzel computes a **single DFT bin** X[k] using a 2nd-order IIR filter:

```
  Initialise: s₁ = s₂ = 0
  Coefficient: c = 2·cos(2πk/N)

  For each sample x[n]:
      s₀ = x[n] + c·s₁ - s₂
      s₂ = s₁
      s₁ = s₀

  After N samples:
      X[k].re = s₁ - s₂·cos(2πk/N)
      X[k].im =     -s₂·sin(2πk/N)
```

**Cost**: 1 multiply + 2 adds per sample (the IIR loop), plus 4 operations
for the final X[k] extraction.  Total: **~N + 4** operations per bin.

### Comparison Table

| N (samples) | FFT: N·log₂(N) | Goertzel: N/bin | Break-even |
|-------------|-----------------|-----------------|------------|
| 256         | 2,048           | 256             | 8 bins     |
| 1,024       | 10,240          | 1,024           | 10 bins    |
| 4,096       | 49,152          | 4,096           | 12 bins    |

**Rule**: Use Goertzel when you need fewer than log₂(N) bins.

---

## DTMF Detection

DTMF (Dual-Tone Multi-Frequency) is the signalling system used by telephone
keypads.  Each key generates **two simultaneous tones** — one from a row
frequency and one from a column frequency:

```
         1209 Hz  1336 Hz  1477 Hz  1633 Hz
  697 Hz    1        2        3        A
  770 Hz    4        5        6        B
  852 Hz    7        8        9        C
  941 Hz    *        0        #        D
```

To detect a key, we only need to check **8 frequencies** — perfect for
Goertzel.  With N=205 samples at 8 kHz (~25 ms):

- FFT cost: 205 × 8 ≈ 1640 operations (and you'd need zero-padding to 256)
- Goertzel cost: 8 × 205 = 1640 operations (exact, no padding needed)

The advantage of Goertzel: **N doesn't need to be a power of 2**.

---

## Generalised Goertzel

The standard Goertzel computes integer bin index k.  The **generalised** form
can detect any arbitrary frequency — just set k to a fractional value:

$$k = \frac{f_{\text{target}} \cdot N}{f_s}$$

This is useful when:
- N and fs don't give clean bins for your target frequency
- You want sub-bin frequency resolution
- You're scanning non-uniformly spaced frequencies

---

## Sliding DFT

For real-time single-frequency tracking with **per-sample updates**, the
Sliding DFT maintains a running DFT bin:

```
  For each new sample x_new replacing oldest x_old:

      X_new[k] = (X_old[k] + x_new - x_old) × e^{j2πk/N}
```

**Cost**: 1 complex multiply + 2 real adds per sample = **O(1)**.

This is ideal for:
- Vibration monitoring (tracking a resonant frequency)
- Real-time pitch tracking
- Lock-in amplifier applications

```
  Input: chirp sweeping 500 → 2000 Hz
  Track: 1000 Hz bin

  |X[k]|:        ╱╲
                 ╱  ╲
                ╱    ╲           ← peak when chirp hits 1 kHz
         ──────╱      ╲──────
         0    sample →      N
```

---

## Demo Walkthrough

### Demo 1: Goertzel vs FFT
Computes X[k] for a 1 kHz tone using both FFT and Goertzel.
Results match to machine precision (~10⁻¹⁴).

### Demo 2: DTMF Detection
Generates DTMF tones for digits 1,5,9,0,*,#,A,D and attempts detection.
All 8 digits should be correctly identified.

![Goertzel DTMF Spectrum](../plots/ch19/goertzel_dtmf.png)

### Demo 3: Generalised Goertzel
Detects a non-integer frequency (1234.5 Hz) with a coarse+fine scan.

### Demo 4: Sliding DFT
Tracks a 1 kHz bin as a chirp sweeps through — magnitude peaks when
the chirp frequency crosses 1 kHz.

![Sliding DFT Tracking](../plots/ch19/sliding_dft.png)
![Goertzel Spectrum](../plots/ch19/goertzel_spectrum.png)

### Demo 5: Complexity Comparison
Prints the operation counts for FFT vs Goertzel at various N values.

---

## Key Takeaways

1. **Goertzel** = single-bin DFT in O(N) — ideal for DTMF, pitch, vibration
2. **Break-even**: Use Goertzel for < log₂(N) bins, FFT otherwise
3. **Generalised Goertzel** works at arbitrary (non-integer) frequencies
4. **Sliding DFT** gives O(1) per-sample bin updates for real-time tracking
5. N does **not** need to be a power of 2 for Goertzel — unlike radix-2 FFT

---

| Navigation | |
|---|---|
| **Previous:** [Chapter 18 — Fixed-Point](18-fixed-point.md) | **Next:** [Chapter 16 — Overlap-Add/Save](16-overlap-add-save.md) |

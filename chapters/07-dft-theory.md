# Chapter 7: The Discrete Fourier Transform (DFT)

> **Prerequisites:** [Ch 03 — Complex Numbers](03-complex-numbers.md), [Ch 05 — Z-Transform](05-z-transform.md)
> **Library module:** [`fft.h`](../include/fft.h) (DFT via FFT algorithm), [`dsp_utils.h`](../include/dsp_utils.h)
> **Runnable demo:** [`07-dft-theory.c`](07-dft-theory.c) — `make chapters && ./build/bin/ch07`

---

## Why This Chapter Matters

The DFT is the **practical bridge** between time and frequency. While the
Z-transform gives us continuous frequency analysis (all of $|z| = 1$), the DFT
samples the frequency axis at N equally-spaced points, giving us a computable,
finite-length frequency representation.

Every spectral analysis, every frequency-domain filter, every modern
communication system uses the DFT (computed efficiently as the FFT).

> **Textbook references:**
> - Oppenheim & Schafer, *Discrete-Time Signal Processing*, Ch 4–5
> - Proakis & Manolakis, *Digital Signal Processing*, Ch 4–5
> - Lyons, *Understanding DSP*, Ch 3
> - Smith, *Scientist & Engineer's Guide to DSP*, Ch 8–9

---

## 7.1 Definition

The **N-point DFT** of a sequence $x[n]$ of length N is:

$$
X[k] = \sum_{n=0}^{N-1} x[n] \, e^{-j2\pi kn/N}, \quad k = 0, 1, \ldots, N-1
$$

The **inverse DFT (IDFT)**:

$$
x[n] = \frac{1}{N} \sum_{k=0}^{N-1} X[k] \, e^{j2\pi kn/N}, \quad n = 0, 1, \ldots, N-1
$$

The DFT maps N time-domain samples to N frequency-domain samples.
It's a **lossless, invertible** transform — no information is gained or lost.

---

## 7.2 The Twiddle Factor

The complex exponential $W_N = e^{-j2\pi/N}$ is called the **twiddle factor**.
Using it, the DFT is compactly written as:

$$
X[k] = \sum_{n=0}^{N-1} x[n] \, W_N^{kn}
$$

Key properties of $W_N$:
- $W_N^N = 1$ (periodicity)
- $W_N^{k+N} = W_N^k$ (periodicity in k)
- $W_N^{N/2} = -1$ (for even N)
- $(W_N^k)^* = W_N^{-k}$ (conjugate symmetry)

These symmetries are what make the FFT algorithm efficient (see Ch 08).

---

## 7.3 Frequency Bin Interpretation

Each DFT bin $X[k]$ corresponds to the frequency:

$$
f_k = \frac{k \cdot f_s}{N} \quad \text{Hz}
$$

or in normalised terms: $\omega_k = 2\pi k / N$.

| Bin k | Frequency | What it represents |
|-------|-----------|-------------------|
| 0 | 0 Hz | DC component (average) |
| 1 | f_s/N | Lowest resolvable frequency |
| N/2 | f_s/2 | Nyquist frequency |
| N/2+1 to N-1 | Negative frequencies | Mirror (for real signals) |

**Frequency resolution:** $\Delta f = f_s / N$. Want finer resolution? Use more samples (larger N) or zero-pad.

---

## 7.4 DFT as Matrix Operation

The DFT can be written as a matrix multiply: $\mathbf{X} = \mathbf{W} \cdot \mathbf{x}$

The DFT matrix $\mathbf{W}$ is an $N \times N$ matrix where entry $(k,n)$ is $W_N^{kn}$.

For N=4:

$$
\mathbf{W}_4 = \begin{pmatrix} 1 & 1 & 1 & 1 \\ 1 & -j & -1 & j \\ 1 & -1 & 1 & -1 \\ 1 & j & -1 & -j \end{pmatrix}
$$

This has $O(N^2)$ complexity. The FFT reduces it to $O(N \log N)$ by exploiting
the structure of $\mathbf{W}$.

---

## 7.5 Properties of the DFT

| Property | Time Domain | DFT Domain |
|----------|-------------|------------|
| Linearity | $a \, x_1[n] + b \, x_2[n]$ | $a \, X_1[k] + b \, X_2[k]$ |
| Time shift | $x[(n-m) \bmod N]$ | $W_N^{mk} X[k]$ |
| Frequency shift | $W_N^{-\ell n} x[n]$ | $X[(k-\ell) \bmod N]$ |
| **Circular convolution** | $(x \circledast h)[n]$ | $X[k] \cdot H[k]$ |
| Parseval's theorem | $\sum |x[n]|^2$ | $\frac{1}{N} \sum |X[k]|^2$ |
| Conjugate symmetry (real x) | $x[n]$ real | $X[N-k] = X^*[k]$ |

---

## 7.6 Circular vs. Linear Convolution

The DFT multiplication property gives **circular** convolution, not linear:

$$
\text{IDFT}\{X[k] \cdot H[k]\} = (x \circledast h)[n] = \sum_{m=0}^{N-1} x[m] \, h[(n-m) \bmod N]
$$

To get **linear** convolution using the DFT:
1. Zero-pad both sequences to length $\geq L_x + L_h - 1$
2. Compute DFTs, multiply, IDFT

This is the basis of **fast convolution** (overlap-add / overlap-save in Ch 16).

---

## 7.7 Parseval's Theorem

Energy is preserved across the transform:

$$
\sum_{n=0}^{N-1} |x[n]|^2 = \frac{1}{N} \sum_{k=0}^{N-1} |X[k]|^2
$$

This means the total energy computed in time = total energy computed in
frequency. It's essential for power spectral analysis.

---

## 7.8 Symmetry for Real Signals

For **real-valued** signals (which is most practical signals):

$$
X[N-k] = X^*[k]
$$

Consequences:
- The magnitude spectrum is **symmetric** about N/2
- You only need to compute/store bins 0 to N/2
- The real-valued FFT (`fft_real`) exploits this for 2x efficiency

---

## 7.9 Zero-Padding and Interpolation

Adding zeros to the end of a signal before the DFT does NOT increase frequency
resolution — the underlying signal hasn't changed. But it does:

1. **Interpolate** the frequency spectrum (finer grid)
2. Make spectral peaks easier to locate
3. Enable linear convolution via circular (by padding to the right length)

True frequency resolution depends only on the signal duration: $\Delta f = 1/T$.

---

## 7.10 Relationship to Other Transforms

```
DTFT          X(e^{jω}) = Σ x[n] e^{-jωn}     — continuous ω, infinite sum
DFT           X[k] = DTFT sampled at ω_k = 2πk/N — N discrete frequencies
Z-Transform   X(z) = Σ x[n] z^{-n}             — DFT is Z-transform on unit circle
FFT           fast algorithm to compute DFT     — O(N log N) vs O(N²)
```

The DFT is the computable version of the DTFT, and the FFT is the efficient
algorithm for computing the DFT.

---

## Exercises

1. **DFT by hand:** Compute the 4-point DFT of `x = {1, 0, -1, 0}`.
   (Answer: X = {0, 2, 0, 2}.)

2. **Circular vs. linear:** DFT-multiply two length-4 sequences without
   zero-padding. Compare with true linear convolution. Where do they differ?

3. **Parseval's check:** Compute time-domain energy and frequency-domain
   energy of a test signal. Verify they're equal.

---

## What's Next?

- [Ch 08: FFT Algorithms](08-fft-fundamentals.md) — Cooley-Tukey Radix-2 DIT, the O(N log N) breakthrough.
- [Ch 09: Window Functions](09-window-functions.md) — controlling spectral leakage in finite-length DFTs.

---

## Generated Plots

> Regenerate with `make plots` from the project root.

### DFT of a Pure Tone
![DFT magnitude spectrum of a single sinusoid](../plots/ch07/dft_spectrum.png)

### Standard Signals DFT
![DFT of impulse, DC, and alternating signals](../plots/ch07/standard_signals_dft.png)

### Zero-Padding Effect
![Spectral interpolation via zero-padding](../plots/ch07/zero_padding.png)

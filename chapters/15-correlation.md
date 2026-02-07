# Chapter 15: Correlation & Autocorrelation

> **Prerequisites:** Ch 07 (DFT), Ch 08 (FFT), Ch 14 (PSD)
> **Library modules:** `correlation.h` / `correlation.c`

---

## Motivation

Correlation measures **similarity** between signals as a function of
time shift (lag).  It answers questions like:

- *Where is this pulse hidden in this noisy recording?*
- *What is the fundamental frequency (pitch) of this sound?*
- *How many samples apart are these two sensor readings?*

---

## 1. Cross-Correlation

The cross-correlation of $x$ and $y$ is:

$$
R_{xy}[m] = \sum_{n} x[n]\,y[n + m]
$$

where $m$ is the **lag** in samples.

**Interpretation:** $R_{xy}[m]$ is large when $x$ and a shifted copy of $y$
overlap well.  The lag of the peak tells you the time offset.

### FFT-based computation

Direct computation is $O(N^2)$.  We use the correlation theorem:

$$
R_{xy} = \mathcal{F}^{-1}\left\{ X^*[k]\,Y[k] \right\}
$$

This is $O(N\log N)$ — the same efficiency as FFT convolution.

**Important:** zero-pad both signals to length $\geq N_x + N_y - 1$ to get
**linear** (not circular) correlation.

---

## 2. Normalised Cross-Correlation

Raw correlation scales with signal amplitude.  Dividing by
$\sqrt{E_x \cdot E_y}$ gives the **normalised** correlation in $[-1, 1]$:

$$
\hat{R}_{xy}[m] = \frac{R_{xy}[m]}{\sqrt{\sum x[n]^2 \cdot \sum y[n]^2}}
$$

A value of $+1$ means perfect match, $-1$ means perfect anti-match,
$0$ means uncorrelated.

---

## 3. Autocorrelation

Autocorrelation is cross-correlation of a signal with itself:

$$
R_{xx}[m] = \sum_{n} x[n]\,x[n + m]
$$

**Properties:**

| Property | Meaning |
|----------|---------|
| $R_{xx}[0]$ is maximum | Signal is most similar to itself at zero lag |
| $R_{xx}[m] = R_{xx}[-m]$ | Symmetric (even function) |
| Normalised $R_{xx}[0] = 1$ | By definition |
| Periodic signals → periodic $R_{xx}$ | Peaks at multiples of the period |

---

## 4. Applications

### Pitch Estimation

A periodic signal (speech, music) has autocorrelation peaks at lag = period.
The fundamental frequency is simply:

$$
f_0 = \frac{f_s}{\text{lag of first peak}}
$$

Demo 3 estimates the pitch of a 440 Hz tone with harmonics and noise.

### White Noise Fingerprint

White noise is uncorrelated with any shifted version of itself:

$$
R_{nn}[m] = \sigma^2 \delta[m]
$$

The normalised autocorrelation is a **delta function** — a sharp spike at
lag 0 and essentially zero everywhere else (Demo 4).

### Time-Delay Estimation

If $y[n] = \alpha\,x[n - d] + \text{noise}$, the cross-correlation peak
occurs at lag $d$.  This is the basis for:

- Radar / sonar ranging
- Acoustic source localisation
- Audio synchronisation

Demo 5 estimates a 73-sample delay between a chirp and its delayed noisy copy.

---

## 5. API Reference

```c
#include "correlation.h"

/* Cross-correlation (raw) */
int xcorr(const double *x, int nx, const double *y, int ny, double *r);

/* Normalised cross-correlation in [-1, 1] */
int xcorr_normalized(const double *x, int nx,
                     const double *y, int ny, double *r);

/* Autocorrelation */
int autocorr(const double *x, int n, double *r);

/* Normalised autocorrelation (lag 0 = 1.0) */
int autocorr_normalized(const double *x, int n, double *r);

/* Find the peak lag */
int xcorr_peak_lag(const double *r, int r_len, int centre);
```

### Output Layout

| Function | Output length | Lag-0 index |
|----------|--------------|-------------|
| `xcorr` (same length) | 2N − 1 | N − 1 |
| `xcorr` (diff length) | Nx + Ny − 1 | Ny − 1 |
| `autocorr` | 2N − 1 | N − 1 |

---

## 6. Relationship to PSD

The **Wiener-Khinchin theorem** connects autocorrelation and PSD:

$$
R_{xx}[m] \;\xleftrightarrow{\;\text{DFT}\;}\; P_{xx}[k]
$$

The PSD is the DFT of the autocorrelation function.  This means:

- A signal with a sharp autocorrelation peak (white noise) has a flat PSD
- A signal with periodic autocorrelation (tonal) has spectral peaks
- Smoothing the PSD (Welch) is equivalent to windowing the autocorrelation

---

## Exercises

1. Generate two copies of a 200 Hz sine wave, delay one by 15 samples,
   add independent noise (σ = 2).  Use `xcorr_normalized` to recover the
   delay.  How much noise can you add before the estimate breaks?

2. Create a "mystery signal" that is the sum of two harmonically related
   tones (e.g. 330 Hz + 660 Hz).  Use `autocorr_normalized` to find the
   fundamental period.  Which autocorrelation peak corresponds to $f_0$?

3. Verify the Wiener-Khinchin theorem numerically: compute the normalised
   autocorrelation of white noise, then take its DFT.  Compare to the
   periodogram of the same noise.

---

## What's Next?

- [Ch 16 — Overlap-Add/Save](16-overlap-add-save.md) — streaming convolution
  using the same FFT primitives.
- [Ch 17 — Multirate DSP](17-multirate-dsp.md) — decimation and interpolation.

---

## Generated Plots

> Regenerate with `make plots` from the project root.

### Pulse Detection
![Cross-correlation detects a pulse hidden in noise](../plots/ch15/pulse_detection.png)

### Normalised Cross-Correlation
![Normalised xcorr of a sine and its phase-shifted copy](../plots/ch15/normalized_xcorr.png)

### Pitch Estimation
![Autocorrelation reveals the pitch period of A4](../plots/ch15/autocorr_pitch.png)

### White Noise Autocorrelation
![Autocorrelation of white noise — impulse at lag 0](../plots/ch15/noise_autocorr.png)

### Time-Delay Estimation
![Cross-correlation peak at the true delay of 73 samples](../plots/ch15/time_delay.png)

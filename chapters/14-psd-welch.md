# Chapter 14: Power Spectral Density & Welch's Method

> **Prerequisites:** Ch 07 (DFT), Ch 08 (FFT), Ch 09 (Windows)
> **Library modules:** `spectrum.h` / `spectrum.c`

---

## Motivation

A single DFT gives you the spectrum of a signal at one snapshot in time.
For **random** or **noisy** signals the raw spectrum is wildly jagged — every
realisation looks different.  We need a statistically reliable estimate of
*how power is distributed across frequency*.  That's the Power Spectral
Density (PSD).

---

## 1. The Periodogram

The simplest PSD estimate squares the DFT magnitudes:

$$
\hat{P}_{xx}[k] = \frac{1}{N}\left| X[k] \right|^2
$$

where $X[k] = \text{DFT}\{x[n]\}$.

**Problem:** the periodogram is an *inconsistent* estimator — its variance
does **not** decrease as $N$ grows.  A single periodogram of a noisy signal
looks almost as noisy as the signal itself (Demo 2 in the code).

### One-sided PSD

For real signals the spectrum is symmetric.  We double the positive-frequency
bins (except DC and Nyquist) to get the **one-sided** PSD — all the power
is accounted for in $k = 0 \ldots N/2$.

---

## 2. Windowed Periodogram

Applying a window before the DFT reduces **spectral leakage** (see Ch 09):

$$
\hat{P}_{xx}[k] = \frac{1}{U}\left| \sum_{n=0}^{N-1} w[n]\,x[n]\,e^{-j2\pi kn/N} \right|^2
$$

where $U = \sum w[n]^2$ normalises for the window's energy.

The Hann window is a good default — it provides −31 dB sidelobe suppression
at the cost of slightly wider main lobe.

---

## 3. Welch's Method

Welch (1967) dramatically reduces variance:

1. **Segment** the signal into overlapping blocks of length $L$
2. **Window** each segment (Hann, Hamming, …)
3. **Periodogram** each segment
4. **Average** all periodograms

$$
\hat{P}_{Welch}[k] = \frac{1}{K} \sum_{i=1}^{K} \hat{P}_i[k]
$$

With 50% overlap and a Hann window, the variance drops by roughly $1/K$
where $K$ is the number of segments.

### Resolution vs Variance Trade-off

| Longer segments | More frequency resolution, fewer averages (higher variance) |
|-----------------|-------------------------------------------------------------|
| Shorter segments | Less frequency resolution, more averages (smoother PSD) |

Demo 4 in the code shows three different segment lengths on the same signal.

---

## 4. Cross-Power Spectral Density (CPSD)

The cross-PSD measures the shared frequency content between two signals:

$$
\hat{P}_{xy}[k] = \frac{1}{K} \sum_{i=1}^{K} X_i^*[k]\,Y_i[k]
$$

**Key insight:** independent noise in $x$ and $y$ averages to zero in the
cross-PSD, while the common signal accumulates.  This is the spectral
analogue of cross-correlation (Ch 15).

---

## 5. API Reference

```c
#include "spectrum.h"

/* Basic periodogram */
int periodogram(const double *x, int n, double *psd, int nfft);

/* Windowed periodogram */
int periodogram_windowed(const double *x, int n, double *psd, int nfft,
                         window_fn win);

/* Welch's averaged PSD */
int welch_psd(const double *x, int n, double *psd, int nfft,
              int seg_len, int overlap, window_fn win);

/* Cross-PSD */
int cross_psd(const double *x, const double *y, int n,
              Complex *cpsd, int nfft,
              int seg_len, int overlap, window_fn win);

/* Helpers */
void psd_to_db(const double *psd, double *psd_db, int n_bins, double floor_db);
void psd_freq_axis(double *freq, int n_bins, double fs);
```

### Parameter Cheat-Sheet

| Parameter | Typical value | Notes |
|-----------|--------------|-------|
| `nfft` | power-of-2, ≥ seg_len | Zero-pads if > seg_len |
| `seg_len` | 256–4096 | Trades resolution ↔ variance |
| `overlap` | seg_len / 2 | 50% is standard with Hann |
| `win` | `hann_window` | Good general-purpose choice |

---

## 6. Practical Tips

1. **Always window** — rectangular windows leak terribly in PSD estimation.
2. **50% overlap + Hann** is the "Welch default" — hard to beat without
   domain-specific knowledge.
3. **Zero-pad** for interpolation, not resolution — padding can't create
   new frequency content.
4. **Normalise correctly** — psd_to_db applies $10\log_{10}$ (power, not
   amplitude), so peaks occur at expected dB levels.
5. **Use cross-PSD** to find shared components between sensor channels
   in the presence of independent noise.

---

## Exercises

1. Generate a signal with three sinusoids at 500, 1500, and 3000 Hz buried
   in Gaussian noise (σ = 5).  Use Welch PSD with different segment lengths
   to find all three tones.

2. Compare the periodogram variance for N = 1024 vs N = 8192.  Does the
   periodogram get smoother?  (Hint: it shouldn't.)

3. Implement **Bartlett's method** (non-overlapping segments, rectangular
   window) and compare to Welch.  Which gives a smoother PSD for the
   same data length?

---

## What's Next?

- [Ch 15 — Correlation & Autocorrelation](15-correlation.md) — time-domain
  analysis of signal similarity, pitch estimation, and delay detection.
- [Ch 16 — Overlap-Add/Save](16-overlap-add-save.md) — efficient convolution
  for streaming data.

---

## Generated Plots

> Regenerate with `make plots` from the project root.

### Periodogram — Clean Signal
![Periodogram of a clean two-tone signal](../plots/ch14/periodogram_clean.png)

### Periodogram — Noisy Signal
![High-variance periodogram of a noisy signal](../plots/ch14/periodogram_noisy.png)

### Welch PSD
![Welch PSD — lower variance estimate](../plots/ch14/welch_psd.png)
![Periodogram vs Welch comparison](../plots/ch14/periodogram_vs_welch.png)

### Resolution vs Variance Trade-off
![Effect of segment length on Welch resolution](../plots/ch14/welch_resolution.png)
![Resolution Tradeoff](../plots/ch14/resolution_tradeoff.png)

### Cross-PSD
![Cross-PSD reveals common frequency content](../plots/ch14/cross_psd.png)

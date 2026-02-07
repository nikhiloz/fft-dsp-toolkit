# Chapter 4 — Digital Filters

Filters shape a signal by passing some frequencies and attenuating
others. This chapter covers **FIR (Finite Impulse Response)** filters —
the simplest and most predictable type.

---

## 4.1 What Does a Filter Do?

A **lowpass filter** passes low frequencies and blocks high frequencies.
A **bandpass filter** passes a range and blocks the rest.

In the time domain, filtering is **convolution** — a weighted running
average where the weights are called **filter coefficients** $h[k]$.

## 4.2 The FIR Filter Equation

$$y[n] = \sum_{k=0}^{M-1} h[k] \cdot x[n - k]$$

where:
- $x[n]$ is the input signal
- $y[n]$ is the output signal
- $h[k]$ are the filter coefficients (the "impulse response")
- $M$ is the filter order (number of taps)

Each output sample is a weighted sum of the current and $M-1$ previous
input samples. This is called "direct-form convolution".

### Why FIR?

- **Always stable** (no feedback → no poles → can't oscillate)
- **Linear phase** (symmetric coefficients → no phase distortion)
- **Simple to design and understand**
- **Predictable latency**: exactly $M/2$ samples of group delay

## 4.3 Implementation Walk-Through

The core filter is in [`src/filter.c`](../src/filter.c) lines 43–57:

```c
void fir_filter(const double *in, double *out, int n,
                const double *h, int order)
{
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int k = 0; k < order; k++) {
            int idx = i - k;
            if (idx >= 0) {
                sum += h[k] * in[idx];
            }
            /* else: x[idx] = 0 (zero-padding) */
        }
        out[i] = sum;
    }
}
```

**What's happening:**
- Outer loop: iterate over each output sample
- Inner loop: multiply-accumulate with the filter coefficients
- `idx = i - k`: look back $k$ samples into the past
- `if (idx >= 0)`: before the signal started, assume zero

This is $O(N \times M)$ — fine for small filters, but for very long
filters (>1000 taps), FFT-based convolution is faster.

### The API

See [`include/filter.h`](../include/filter.h) for the full interface:

```c
void fir_filter(const double *in, double *out, int n,
                const double *h, int order);
void fir_moving_average(double *h, int taps);
void fir_lowpass(double *h, int taps, double cutoff);
```

## 4.4 Moving Average Filter

The simplest lowpass: every coefficient is $1/M$.

$$h[k] = \frac{1}{M}, \quad k = 0, 1, \ldots, M-1$$

This averages the last $M$ samples, smoothing out high-freq noise.

See [`src/filter.c`](../src/filter.c) lines 72–76:

```c
void fir_moving_average(double *h, int taps) {
    double coeff = 1.0 / taps;
    for (int i = 0; i < taps; i++) {
        h[i] = coeff;
    }
}
```

**Frequency response:** sinc-shaped. Has nulls at multiples of $f_s / M$
and poor side-lobe suppression (−13 dB). Good for quick smoothing, not
for precision filtering.

## 4.5 Windowed-Sinc Lowpass Design

For a proper lowpass filter, we start from the **ideal** impulse
response and make it practical.

### The Ideal Lowpass

The ideal lowpass with cutoff $f_c$ (normalised to $f_s$) has an
impulse response that is a **sinc function**:

$$h_{\text{ideal}}[k] = \frac{\sin(2\pi f_c k)}{\pi k}$$

Problem: sinc extends to $\pm\infty$ — we can't use infinite
coefficients.

### Truncation + Windowing

1. **Truncate** to $M$ taps (keep only $-M/2$ to $+M/2$)
2. **Apply a window** (Hamming) to taper the truncated edges
3. **Normalise** so the DC gain equals 1.0

See [`src/filter.c`](../src/filter.c) lines 95–118:

```c
void fir_lowpass(double *h, int taps, double cutoff) {
    int center = taps / 2;
    double sum = 0.0;

    for (int i = 0; i < taps; i++) {
        int offset = i - center;
        if (offset == 0) {
            h[i] = 2.0 * cutoff;           /* value at centre */
        } else {
            h[i] = sin(2.0 * M_PI * cutoff * offset) / (M_PI * offset);
        }
        h[i] *= hamming_window(taps, i);    /* apply window */
        sum += h[i];
    }

    /* Normalise for unity DC gain */
    if (sum != 0.0) {
        for (int i = 0; i < taps; i++) h[i] /= sum;
    }
}
```

**Note the cross-reference:** The `hamming_window` function comes from
`dsp_utils.c` (Chapter 03). Filter design *requires* window functions
— this is where the two concepts connect.

### Choosing the Number of Taps

More taps = sharper transition band, but:
- Higher computational cost per sample
- Longer settling time (group delay = $M/2$ samples)
- More memory

Rule of thumb: $M \approx 4 / (\text{transition width as fraction of } f_s)$.

## 4.6 The Filter Demo

The demo in [`04-digital-filters.c`](04-digital-filters.c) shows
a complete noise-reduction workflow:

1. Generate a clean 200 Hz sine wave
2. Add high-frequency noise
3. Design a 31-tap lowpass filter at 500 Hz cutoff
4. Filter the noisy signal
5. Compare RMS before and after

```bash
make
./build/bin/filter_demo
```

Output:
```
Signal Analysis:
  Clean signal RMS:    0.7108
  Noisy signal RMS:    0.7556  (noise added 6.3%)
  Filtered signal RMS: 0.6637  (after settling)
```

The filter successfully removes the high-frequency components while
preserving the 200 Hz signal.

## 4.7 Testing

Six tests verify the filter in [`tests/test_filter.c`](../tests/test_filter.c):

| Test | What It Verifies |
|------|-----------------|
| Identity passthrough | 1-tap filter `h={1}` passes signal unchanged |
| Zero input → zero output | No energy from nowhere |
| Impulse response | Delta → output matches coefficients exactly |
| Moving average step | Step function → smooth ramp output |
| Lowpass coeffs sum to 1.0 | DC gain is unity |
| Lowpass attenuates high freq | 3500 Hz signal reduced by >20 dB |

```bash
make test
```

## 4.8 Exercises

1. **Experiment:** Change the filter demo to use 63 taps instead of 31.
   Does the settling time increase? Does the filtering get sharper?

2. **Code exercise:** Implement a **highpass filter** by subtracting
   the lowpass from a delta function:
   $h_{\text{hp}}[k] = \delta[k - M/2] - h_{\text{lp}}[k]$

3. **Code exercise:** Implement a **bandpass filter** by convolving a
   lowpass with a highpass, or by modulating a lowpass to a centre
   frequency.

4. **Thinking question:** Why does the filter demo skip the first
   `TAPS` samples when computing the filtered RMS? *(Hint: transient/
   settling time.)*

5. **Paper exercise:** Calculate the group delay of a 31-tap
   linear-phase FIR filter at 8000 Hz sample rate.

---

**Previous:** [Chapter 03 — Window Functions](03-window-functions.md)
| **Next:** [Chapter 05 — Spectral Analysis →](05-spectral-analysis.md)

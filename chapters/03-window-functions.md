# Chapter 3 — Window Functions

If you take an FFT of a raw chunk of signal, the result will be smeared
and messy. Window functions fix this. Understanding *why* they work is
essential for getting clean spectral analysis.

---

## 3.1 The Problem: Spectral Leakage

When we capture $N$ samples from a continuous signal, we are
implicitly multiplying by a rectangular window (1 inside, 0 outside).
This abrupt truncation causes **spectral leakage** — energy from one
frequency smears into neighbouring bins.

**Example:** A pure 440 Hz sine computed with a 256-point FFT at
8000 Hz sample rate. The ideal frequency bin is $440 / 31.25 = 14.08$
— it falls *between* bins 14 and 15. Without a window, energy spreads
across many bins.

This happens because the DFT assumes the signal is *periodic*. If the
captured chunk doesn't contain an exact integer number of cycles, the
DFT sees a discontinuity at the edges.

## 3.2 The Solution: Windowing

Multiply the signal by a function that tapers smoothly to zero at both
edges. This removes the discontinuity and concentrates energy into
narrower peaks.

The trade-off: **better side-lobe suppression ↔ wider main lobe**
(poorer frequency resolution).

## 3.3 Three Windows Implemented

All window functions live in [`src/dsp_utils.c`](../src/dsp_utils.c) lines 74–105.

### Hann Window

$$w[i] = 0.5 \left(1 - \cos\!\frac{2\pi i}{N-1}\right)$$

```c
double hann_window(int n, int i) {
    return 0.5 * (1.0 - cos(2.0 * M_PI * i / (n - 1)));
}
```

- **Side-lobe level:** −31 dB
- **Main-lobe width:** 4 bins
- **Use when:** General-purpose spectral analysis

### Hamming Window

$$w[i] = 0.54 - 0.46 \cos\!\frac{2\pi i}{N-1}$$

```c
double hamming_window(int n, int i) {
    return 0.54 - 0.46 * cos(2.0 * M_PI * i / (n - 1));
}
```

- **Side-lobe level:** −42 dB (better than Hann)
- **Main-lobe width:** 4 bins
- **Use when:** You need better side-lobe suppression and don't need
  the window to reach exactly zero at the edges

### Blackman Window

$$w[i] = 0.42 - 0.5\cos\!\frac{2\pi i}{N-1} + 0.08\cos\!\frac{4\pi i}{N-1}$$

```c
double blackman_window(int n, int i) {
    double t = 2.0 * M_PI * i / (n - 1);
    return 0.42 - 0.5 * cos(t) + 0.08 * cos(2.0 * t);
}
```

- **Side-lobe level:** −58 dB (excellent)
- **Main-lobe width:** 6 bins (wider — poorer resolution)
- **Use when:** Detecting weak signals near strong ones

### Comparison Table

| Window | Main Lobe (bins) | Side Lobe (dB) | Best For |
|--------|----------------:|---------------:|----------|
| Rectangular (none) | 2 | −13 | Coherent sampling only |
| Hann | 4 | −31 | General analysis |
| Hamming | 4 | −42 | Filter design |
| Blackman | 6 | −58 | Weak signal detection |

## 3.4 Applying a Window

The `apply_window` function uses a function pointer so any window
can be plugged in. See [`include/dsp_utils.h`](../include/dsp_utils.h)
lines 39–40 and [`src/dsp_utils.c`](../src/dsp_utils.c) lines 107–110:

```c
/* Header: typedef and function signature */
typedef double (*window_fn)(int n, int i);
void apply_window(double *signal, int n, window_fn w);

/* Implementation: multiply each sample by the window value */
void apply_window(double *signal, int n, window_fn w) {
    for (int i = 0; i < n; i++) {
        signal[i] *= w(n, i);
    }
}
```

Usage in the chapter demo ([`03-window-functions.c`](03-window-functions.c)):

```c
/* Apply Hann window before FFT to reduce leakage */
apply_window(signal, N, hann_window);
```

> **Design note:** The function pointer `window_fn` avoids code
> duplication. You swap windows by changing one argument, not rewriting
> the loop.

## 3.5 Windows in Filter Design

Windows aren't only for FFT analysis — they appear in filter design too.
The [windowed-sinc lowpass](04-digital-filters.md) in Chapter 04 uses
the Hamming window to taper the ideal sinc impulse response.

See [`src/filter.c`](../src/filter.c) lines 103–105 where the filter
coefficients are multiplied by `hamming_window()`:

```c
/* Apply Hamming window to the sinc function */
h[i] *= hamming_window(taps, i);
```

## 3.6 Try It Yourself

```bash
make
./build/bin/fft_demo
```

The demo uses `hann_window`. Try editing [`03-window-functions.c`](03-window-functions.c)
to use `hamming_window` or `blackman_window` instead, rebuild, and
compare the output. The peaks should change width.

## 3.7 Exercises

1. **Experiment:** Modify `fft_demo.c` to run the FFT *without* any
   window (comment out the `apply_window` line). How much does 440 Hz
   spread into adjacent bins?

2. **Paper exercise:** Plot the Hann window for $N = 16$ by computing
   $w[i]$ for $i = 0 \ldots 15$.

3. **Code exercise:** Implement a **Kaiser window** (requires the
   zeroth-order modified Bessel function $I_0$). The Kaiser window has
   an adjustable $\beta$ parameter that trades main-lobe width for
   side-lobe suppression.

4. **Thinking question:** Why does the Hamming window not reach zero
   at the edges ($w[0] = 0.08$)? What advantage does this give?

---

**Previous:** [Chapter 02 — FFT Fundamentals](02-fft-fundamentals.md)
| **Next:** [Chapter 04 — Digital Filters →](04-digital-filters.md)

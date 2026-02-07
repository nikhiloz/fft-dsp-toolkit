# Chapter 5 — Spectral Analysis

This chapter brings together the FFT (Chapter 02), window functions
(Chapter 03), and utility functions (Chapter 01) to perform real
spectral analysis — the kind used in audio apps, radar, and
scientific instruments.

---

## 5.1 The Complete Analysis Pipeline

```
signal → window → FFT → magnitude → dB → display
```

> **📊 Signal Processing Pipeline** — [View full-size diagram →](../reference/diagrams/signal_flow.png)

Each stage has a purpose:
1. **Acquire** $N$ samples from the signal
2. **apply_window()** — reduce spectral leakage (Ch. 03)
3. **fft()** — transform to frequency domain (Ch. 02)
4. **fft_magnitude()** — extract $|X[k]|$ from complex bins
5. **db_from_magnitude()** — convert to logarithmic dB scale
6. **Display** — print or plot the spectrum

## 5.2 Walk-Through: The FFT Demo

[`05-spectral-analysis.c`](05-spectral-analysis.c) implements this full
pipeline. Let's trace each step.

### Signal Generation (lines 25–35)

```c
#define N     256
#define FS    8000.0
#define F1    440.0
#define F2    1000.0

for (int i = 0; i < N; i++) {
    double t = (double)i / FS;
    signal[i] = sin(2.0 * M_PI * F1 * t)
              + sin(2.0 * M_PI * F2 * t);
}
```

Creates a signal with two pure tones: 440 Hz (musical A4) and 1000 Hz.

### Windowing (line ~40)

```c
apply_window(signal, N, hann_window);
```

The Hann window (−31 dB side lobes) tapers the signal edges. Without
this, 440 Hz energy would smear into bins 12–16 instead of clustering
near bin 14.

### FFT + Magnitude Extraction (lines ~43–50)

```c
Complex spectrum[N];
fft_real(signal, spectrum, N);

double mag[N];
fft_magnitude(spectrum, mag, N);
```

`fft_real()` wraps the real signal into complex format and runs the
in-place Cooley-Tukey FFT. Then `fft_magnitude()` computes
$|X[k]| = \sqrt{\text{re}^2 + \text{im}^2}$ for each bin.

### Magnitude → dB + Display (lines ~55–70)

```c
double bin_hz = FS / N;   /* 31.25 Hz per bin */

for (int k = 0; k <= N / 2; k++) {
    double freq = k * bin_hz;
    double db   = db_from_magnitude(mag[k]);
    if (db > -25.0) {
        printf("%8.1f Hz     |  %6.1f dB", freq, db);
    }
}
```

Only bins with energy above −25 dB are printed. The frequency axis
is `k * fs / N`.

### Expected Output

```
   437.5 Hz     |   -6.1 dB  ◄── peak!
   468.8 Hz     |  -11.0 dB  ◄── peak!
  1000.0 Hz     |  -12.1 dB  ◄── peak!
```

The 440 Hz peak falls between bins 14 and 15 (437.5 and 468.8 Hz)
because $440 / 31.25 = 14.08$ — not exactly on a bin. This is
**spectral leakage in action** — the Hann window keeps it to 2 bins
instead of 10+.

## 5.3 Key Concepts

### Frequency Resolution

$$\Delta f = \frac{f_s}{N}$$

For our demo: $\Delta f = 8000 / 256 = 31.25$ Hz.

To improve resolution, either:
- **Increase $N$** (more samples → finer bins, but more computation)
- **Decrease $f_s$** (only if your signal bandwidth allows it)

### The dB Scale

Human perception is logarithmic, so we use decibels:

$$\text{dB} = 20 \cdot \log_{10}(|X[k]|)$$

See [`src/dsp_utils.c`](../src/dsp_utils.c) lines 122–125:

```c
double db_from_magnitude(double mag) {
    if (mag <= 0.0) return -200.0;  /* floor for silence */
    return 20.0 * log10(mag);
}
```

A −200 dB floor prevents `log10(0)` from returning `-inf`.

### Nyquist Limit

Only bins 0 through $N/2$ contain unique information for real signals.
Bin $N/2$ corresponds to $f_s / 2$ — the **Nyquist frequency**.

Any signal component above $f_s / 2$ will **alias** to a lower
frequency and corrupt the spectrum. Always ensure your sample rate is
at least 2× the highest frequency of interest.

## 5.4 RMS and Signal Power

The **Root Mean Square** measures the "effective amplitude" of a signal:

$$\text{RMS} = \sqrt{\frac{1}{N} \sum_{n=0}^{N-1} x[n]^2}$$

See [`src/dsp_utils.c`](../src/dsp_utils.c) lines 128–133. The filter
demo uses RMS to quantify noise reduction (Chapter 04).

## 5.5 Try It Yourself

```bash
make && ./build/bin/fft_demo
```

Modifications to try:
- Change `N` to 512 or 1024 — observe sharper peaks
- Change `F1` to 500 Hz (exact multiple of bin width) — observe a
  single clean peak with no leakage
- Remove the window — observe leakage spreading everywhere

## 5.6 Exercises

1. **Experiment:** Generate a chirp signal (frequency sweeping linearly
   from 100 Hz to 3000 Hz over $N$ samples). What does the FFT show?

2. **Code exercise:** Add a `find_peak_frequency()` function that
   returns the frequency of the highest-magnitude bin.

3. **Code exercise:** Implement **zero-padding** — append zeros to
   increase $N$ before the FFT. Does it improve resolution or just
   interpolate between bins?

4. **Thinking question:** Why do we only print bins 0 to $N/2$? What
   happens at bin $N/2 + 1$?

5. **Advanced:** Implement a simple **spectrogram** — divide a long
   signal into overlapping frames, FFT each frame, and print the
   time-frequency map.

---

**Previous:** [Chapter 04 — Digital Filters](04-digital-filters.md)
| **Next:** [Chapter 06 — Real-Time Streaming →](06-real-time-streaming.md)

# Chapter 2: Sampling, Aliasing & the Nyquist Theorem

> **Prerequisites:** [Ch 01 — Signals & Sequences](01-signals-and-sequences.md)
> **Library module:** [`signal_gen.h`](../include/signal_gen.h) (signal generators used for demos)
> **Runnable demo:** [`02-sampling-and-aliasing.c`](02-sampling-and-aliasing.c) — `make chapters && ./build/bin/ch02s`

---

## Why This Chapter Matters

Every digital signal begins its life as an analog signal. Sampling is the bridge
between the continuous and discrete worlds. Get it wrong, and information is
**irreversibly destroyed** — frequencies fold on top of each other, and no
amount of clever processing can untangle them.

This chapter answers three questions:
1. How fast must we sample?
2. What happens when we sample too slowly?
3. Can we reconstruct the original continuous signal from samples?

> **Textbook references:**
> - Oppenheim & Willsky, *Signals and Systems*, Ch 7
> - Oppenheim & Schafer, *Discrete-Time Signal Processing*, Ch 4
> - Proakis & Manolakis, *Digital Signal Processing*, Ch 1
> - Lyons, *Understanding DSP*, Ch 1–2
> - Smith, *Scientist & Engineer's Guide to DSP*, Ch 3

---

## 2.1 From Continuous to Discrete: The Sampling Process

A continuous signal `x(t)` is sampled at intervals of `T` seconds to produce
a discrete sequence:

$$
x[n] = x(nT), \quad n = 0, 1, 2, \ldots
$$

The **sample rate** (or sampling frequency) is:

$$
f_s = \frac{1}{T} \quad \text{(samples per second, Hz)}
$$

In C, the sample rate is just a number we track alongside the data:

```c
double signal[256];
double sample_rate = 44100.0;  /* CD-quality audio */
gen_sine(signal, 256, 1.0, 440.0, sample_rate, 0.0);
```

The signal array knows nothing about time — it's our job to remember that
`signal[n]` corresponds to time `t = n / sample_rate`.

---

## 2.2 The Nyquist-Shannon Sampling Theorem

> **Theorem (Shannon, 1949):** A bandlimited signal with no frequency
> component above `f_max` Hz can be perfectly reconstructed from its
> samples if and only if:
>
> $$f_s > 2 f_{\max}$$

The critical frequency $f_N = f_s / 2$ is called the **Nyquist frequency**.
It is the highest frequency that can be represented at sample rate $f_s$.

| Sample Rate | Nyquist Frequency | Can represent up to |
|-------------|-------------------|---------------------|
| 8,000 Hz | 4,000 Hz | Telephone speech |
| 44,100 Hz | 22,050 Hz | CD audio (full hearing range) |
| 48,000 Hz | 24,000 Hz | Professional audio |
| 96,000 Hz | 48,000 Hz | High-res audio |

**Key insight:** The sample rate doesn't just limit resolution — it defines a
**hard ceiling** on what frequencies can exist in the digital domain.

---

## 2.3 Aliasing: When Sampling Goes Wrong

If a signal contains frequencies above `f_s/2`, those frequencies **fold back**
into the range `[0, f_s/2]` and become indistinguishable from lower frequencies.
This is **aliasing**.

$$
\text{A signal at frequency } f \text{ aliases to } |f - k \cdot f_s|
$$

where `k` is the integer that brings the result into `[0, f_s/2]`.

### The Wagon Wheel Effect

Think of a movie showing a spinning wheel. At 24 frames per second, a wheel
spinning at 23 revolutions per second looks like it's going forward at 23 rps.
But at 25 rps, it looks like it's going **backward** at 1 rps. That's aliasing —
the sample rate (24 fps) can't distinguish 23 from 25.

### In DSP terms:

At sample rate 1000 Hz (Nyquist = 500 Hz):
- A 300 Hz tone appears correctly at 300 Hz
- A 700 Hz tone aliases to **300 Hz** (700 = 1000 - 300)
- A 1300 Hz tone also aliases to **300 Hz** (1300 = 1000 + 300)

**They all produce identical samples.** The demo below proves this.

---

## 2.4 Normalised Frequency

In discrete-time DSP, we often work with **normalised frequency**:

$$
\omega = 2\pi \frac{f}{f_s} \quad \text{(radians per sample)}
$$

Or equivalently, $f/f_s$ (cycles per sample, range [0, 0.5] for real signals).

The advantage: algorithms become independent of sample rate. An FFT doesn't
know or care whether `f_s` is 8 kHz or 96 kHz — it works in normalised
frequency internally.

| Physical Frequency | Normalised (ω) | Meaning |
|-------------------|-----------------|---------|
| 0 Hz | 0 | DC (constant) |
| f_s/4 | π/2 | Quarter of Nyquist |
| f_s/2 | π | Nyquist frequency |
| f_s | 2π | Wraps back to DC |

---

## 2.5 Anti-Aliasing Filters

In real hardware, an **anti-aliasing filter** (a lowpass analog filter) is placed
before the ADC to remove frequencies above `f_s/2`. Without it, those
frequencies alias down and corrupt the signal irreversibly.

```
Analog signal  →  [Anti-alias LPF]  →  [ADC: sample at f_s]  →  x[n]
                   cutoff = f_s/2
```

In this tutorial, our signals are already discrete (generated in C), so we
don't need a hardware anti-aliasing filter. But the concept matters when you
connect DSP code to real-world sensors or audio interfaces.

---

## 2.6 Quantization

Sampling in time (discretisation) is one axis. Sampling in **amplitude**
(quantization) is the other. An ADC maps continuous voltage to a fixed set
of integer levels:

$$
x_q[n] = \text{round}(x[n] \times 2^{B-1})
$$

where `B` is the number of bits.

| Bits | Levels | Dynamic Range (dB) | Typical Use |
|------|--------|---------------------|-------------|
| 8 | 256 | 48 dB | Low-quality audio, sensors |
| 16 | 65,536 | 96 dB | CD audio |
| 24 | 16.7M | 144 dB | Professional audio |
| 32 | float | ~150 dB | Scientific, internal processing |

**Quantization noise** — the rounding error — acts like additive white noise
with an RMS level of approximately:

$$
\sigma_q = \frac{\Delta}{\sqrt{12}}
$$

where $\Delta = 2 / 2^B$ is the step size (for input normalised to [-1, +1]).

We'll explore quantization effects deeply in [Ch 18: Fixed-Point Arithmetic](18-fixed-point.md).

---

## 2.7 Reconstruction (D/A Conversion)

The sampling theorem guarantees perfect reconstruction using **sinc interpolation**:

$$
x(t) = \sum_{n=-\infty}^{\infty} x[n] \cdot \text{sinc}\left(\frac{t - nT}{T}\right)
$$

where $\text{sinc}(u) = \sin(\pi u) / (\pi u)$.

In practice, DACs use simpler methods (zero-order hold, then analog smoothing
filter), but the theoretical result is important: **no information is lost by
sampling, as long as the Nyquist condition is met.**

---

## Exercises

1. **Predict the alias:** At `f_s = 8000 Hz`, what frequency does a 5000 Hz
   tone alias to? (Answer: 3000 Hz, because 5000 = 8000 - 3000.)

2. **Critical sampling:** Generate a sine at exactly `f_s/2`. What do you see?
   (The signal alternates +1, -1, +1, -1... — the maximum oscillation possible.)

3. **Quantization experiment:** Generate a sine wave, quantize it to 4 bits
   (16 levels), and measure the resulting noise floor. Compare with the
   theoretical $6.02B + 1.76$ dB formula.

---

## What's Next?

- [Ch 03: Complex Numbers](01-complex-numbers.md) — the mathematical engine that makes frequency analysis possible.
- [Ch 04: LTI Systems & Convolution](04-lti-systems.md) — how filters process signals sample by sample.

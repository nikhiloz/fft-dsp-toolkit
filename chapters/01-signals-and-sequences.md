# Chapter 1: Discrete-Time Signals & Sequences

> **Prerequisites:** Basic calculus, familiarity with `sin`/`cos`.
> **Library module:** [`signal_gen.h`](../include/signal_gen.h) / [`signal_gen.c`](../src/signal_gen.c)
> **Runnable demo:** [`01-signals-and-sequences.c`](01-signals-and-sequences.c) — `make chapters && ./build/bin/ch01`

---

## Why Start Here?

Every DSP system begins with a **signal** — a sequence of numbers. Before we can
analyse frequencies (FFT), design filters (FIR/IIR), or estimate spectra, we need
to understand the raw material we're working with.

This chapter introduces the five fundamental discrete-time signals that appear
throughout DSP theory and practice. Master these, and every later chapter
becomes easier to follow.

> **Textbook references:**
> - Oppenheim & Willsky, *Signals and Systems*, Ch 1–2
> - Oppenheim & Schafer, *Discrete-Time Signal Processing*, Ch 2
> - Lyons, *Understanding DSP*, Ch 1
> - Smith, *Scientist & Engineer's Guide to DSP*, Ch 2–4

---

## 1.1 What Is a Discrete-Time Signal?

A **continuous-time** signal `x(t)` is defined for every real value of `t`.
We cannot store it in a computer — we'd need infinite memory.

A **discrete-time** signal `x[n]` is defined only at integer indices `n`. It is
a **sequence** of numbers. This is what computers actually work with.

```
x[n]:   ... 0  0  0.5  1.0  0.8  0.3  0  0 ...
             ↑        ↑              ↑
            n=-2     n=0            n=3
```

In C, we represent `x[n]` as a `double` array. The index `n` maps directly to
the array index. Simple.

---

## 1.2 The Unit Impulse (Kronecker Delta)

$$
\delta[n] = \begin{cases} 1 & \text{if } n = 0 \\ 0 & \text{otherwise} \end{cases}
$$

The impulse is the simplest non-trivial signal — a single `1` surrounded by zeros.
Despite its simplicity, it is **the most important signal in DSP**:

- The **impulse response** `h[n]` of a system completely characterises that system.
- The FFT of an impulse is flat (all frequencies equal) — see [Ch 08](02-fft-fundamentals.md).
- Convolution with an impulse gives back the original signal: `x[n] * δ[n] = x[n]`.

```c
/* From signal_gen.c */
void gen_impulse(double *out, int n, int delay)
{
    memset(out, 0, n * sizeof(double));
    if (delay >= 0 && delay < n)
        out[delay] = 1.0;
}
```

A **shifted impulse** `δ[n - d]` places the `1` at index `d`. This is used to
model delayed copies of signals.

---

## 1.3 The Unit Step

$$
u[n] = \begin{cases} 1 & \text{if } n \geq 0 \\ 0 & \text{if } n < 0 \end{cases}
$$

The step represents a signal that "turns on" and stays on. It is the running
sum of the impulse:

$$
u[n] = \sum_{k=-\infty}^{n} \delta[k]
$$

Conversely, the impulse is the **first difference** of the step:

$$
\delta[n] = u[n] - u[n-1]
$$

**In DSP practice:** The step function models system startup transients. When you
apply a step to a filter, the output shows how the filter *settles* — this is
the **step response**, closely related to the impulse response.

---

## 1.4 The Real Exponential

$$
x[n] = A \cdot r^n
$$

- If `|r| < 1`: signal **decays** (damped oscillation, charging capacitor)
- If `|r| > 1`: signal **grows** (unstable system — a red flag in filter design)
- If `r = 1`: constant signal (DC)
- If `r = -1`: alternating `+1, -1, +1, -1, …` (Nyquist frequency!)

```c
void gen_exponential(double *out, int n, double amplitude, double base)
{
    double val = amplitude;
    for (int i = 0; i < n; i++) {
        out[i] = val;
        val *= base;
    }
}
```

**Why it matters:** Exponential signals are the **eigenfunctions of LTI systems**.
If you put an exponential into a linear system, you get the same exponential
back, scaled. This property is the foundation of Z-transform analysis
(see [Ch 05](05-z-transform.md)).

---

## 1.5 The Sinusoid

$$
x[n] = A \cos(2\pi f n / f_s + \phi)
$$

Where:
- `A` = amplitude (peak value)
- `f` = frequency in Hz (cycles per second)
- `f_s` = sample rate (samples per second)
- `φ` = initial phase in radians
- `ω = 2πf/f_s` = **normalised angular frequency** (radians per sample)

The sinusoid is the **atomic building block of frequency analysis**. Fourier's key
insight: any signal can be decomposed into a sum of sinusoids at different
frequencies. The FFT (Ch 08) computes exactly this decomposition.

**Key insight for discrete time:** Unlike continuous sinusoids, discrete-time
sinusoids are **periodic only if `f/f_s` is rational**. If `f/f_s = k/N`, the
signal repeats every `N` samples.

---

## 1.6 The Complex Exponential

$$
x[n] = A \cdot e^{j(\omega n + \phi)} = A[\cos(\omega n + \phi) + j\sin(\omega n + \phi)]
$$

This is Euler's formula in action (see [Ch 03](01-complex-numbers.md)). The
complex exponential wraps a cosine (real part) and sine (imaginary part) into
one object. It traces a circle in the complex plane.

**Why use complex exponentials instead of sines?**
- The math is cleaner — multiplication instead of trig identities.
- The DFT/FFT work with complex exponentials natively.
- Phase information is preserved automatically.

```c
void gen_complex_exp(Complex *out, int n, double amplitude,
                     double freq_hz, double sample_rate, double phase_rad)
{
    double omega = 2.0 * M_PI * freq_hz / sample_rate;
    for (int i = 0; i < n; i++) {
        double angle = omega * i + phase_rad;
        out[i] = complex_from_polar(amplitude, angle);
    }
}
```

---

## 1.7 Composite Signals

### Multi-tone (sum of sinusoids)

Any periodic signal can be built from sinusoids (**Fourier series**):

$$
x[n] = \sum_{k=0}^{K-1} A_k \sin(2\pi f_k n / f_s)
$$

This is what `gen_multi_tone()` does. The FFT will later reveal these individual
components — that's the entire point of spectral analysis.

### Chirp (swept frequency)

A chirp sweeps from frequency `f₀` to `f₁` over time:

$$
x[n] = A \sin\left(2\pi \left[f_0 t + \frac{(f_1 - f_0) t^2}{2T}\right]\right)
$$

Chirps are used in:
- Radar (pulse compression)
- Speaker/room measurement (sine sweep)
- System identification (finding the impulse response)

### Noise

White noise has equal power at all frequencies — a flat spectrum. It is the DSP
equivalent of "everything at once." We use it to:
- Test filter frequency response (apply filter to noise, measure output spectrum)
- Model real-world interference (sensor noise, communication channel noise)
- Stress-test algorithms

`gen_gaussian_noise()` uses the **Box-Muller transform** to generate Gaussian
samples from uniform random numbers — a classic numerical technique.

---

## 1.8 Signal Operations

Two basic operations you'll see throughout this tutorial:

**Addition:** `signal_add(a, b, n)` — adds `b[i]` into `a[i]`. Used to combine
tones, add noise, or sum filter outputs.

**Scaling:** `signal_scale(x, n, scale)` — multiplies every sample by `scale`.
Used for normalisation, gain control, and amplitude adjustment.

These operations are **linear** — they satisfy superposition. This is the key
property that makes Fourier analysis, filtering, and convolution work.
(See [Ch 04](04-lti-systems.md) for the full LTI story.)

---

## Exercises

1. **Frequency aliasing preview:** Generate a 440 Hz sine at sample rate 1000 Hz.
   Now generate a 1440 Hz sine at the same rate. Print both — they're identical!
   Why? (Answer: Ch 02 on sampling and aliasing.)

2. **Build a square wave** by summing odd harmonics: `sin(ωn) + sin(3ωn)/3 + sin(5ωn)/5 + …`.
   Use `gen_multi_tone()` with frequencies `f, 3f, 5f, ...` and amplitudes `1, 1/3, 1/5, ...`.
   How many harmonics do you need before it looks "square"?

3. **Exponential decay + sinusoid:** Generate `x[n] = 0.95^n * sin(2π·100·n/1000)`.
   This is what a plucked guitar string sounds like — a decaying oscillation.

---

## What's Next?

- [Ch 02: Sampling, Aliasing & the Nyquist Theorem](02-sampling-and-aliasing.md) — why discrete-time signals exist, and the rules that govern them.
- [Ch 03: Complex Numbers](01-complex-numbers.md) — the mathematical engine behind frequency analysis.

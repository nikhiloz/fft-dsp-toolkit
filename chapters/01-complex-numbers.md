# Chapter 1 — Complex Numbers

Complex numbers are the language of frequency-domain analysis.
Every FFT input and output is a complex number. Before we can understand
the FFT, we need to be fluent in complex arithmetic.

---

## 1.1 What Is a Complex Number?

A complex number $z$ has a **real** part and an **imaginary** part:

$$z = a + bi$$

where $i = \sqrt{-1}$.

In C we represent this as a struct (see [`include/dsp_utils.h`](../include/dsp_utils.h) line 17):

```c
typedef struct {
    double re;  /* Real part      */
    double im;  /* Imaginary part */
} Complex;
```

### Why not use C99 `<complex.h>`?

The standard library provides `double _Complex`, but:
- The struct approach is more explicit and educational
- We can see exactly what each operation does
- Portable across compilers with varying `<complex.h>` support

## 1.2 Polar Form

Any complex number can also be written in **polar form**:

$$z = r \cdot e^{i\theta} = r(\cos\theta + i\sin\theta)$$

where:
- $r = |z| = \sqrt{a^2 + b^2}$ is the **magnitude**
- $\theta = \text{atan2}(b, a)$ is the **phase angle**

This is crucial for the FFT because the "twiddle factors" are complex
exponentials in polar form: $W_N^k = e^{-j \cdot 2\pi k / N}$.

### Implementation

See [`src/dsp_utils.c`](../src/dsp_utils.c) lines 46–57:

```c
double complex_mag(Complex z) {
    return sqrt(z.re * z.re + z.im * z.im);
}

double complex_phase(Complex z) {
    return atan2(z.im, z.re);
}

Complex complex_from_polar(double mag, double phase) {
    return (Complex){ mag * cos(phase), mag * sin(phase) };
}
```

The `complex_from_polar` function converts $(r, \theta) \to (a, b)$ using
Euler's formula. This is exactly what the FFT does when it creates
twiddle factors.

## 1.3 Complex Arithmetic

### Addition & Subtraction

Add the real parts together, imaginary parts together:

$$(a + bi) + (c + di) = (a+c) + (b+d)i$$

See [`src/dsp_utils.c`](../src/dsp_utils.c) lines 25–30:

```c
Complex complex_add(Complex a, Complex b) {
    return (Complex){ a.re + b.re, a.im + b.im };
}

Complex complex_sub(Complex a, Complex b) {
    return (Complex){ a.re - b.re, a.im - b.im };
}
```

### Multiplication (The Important One)

$$(a + bi)(c + di) = (ac - bd) + (ad + bc)i$$

This is the **most performance-critical operation** in the FFT.
Every butterfly computes one complex multiply with the twiddle factor.
A 1024-point FFT does $512 \times 10 = 5120$ complex multiplies.

See [`src/dsp_utils.c`](../src/dsp_utils.c) lines 38–43:

```c
Complex complex_mul(Complex a, Complex b) {
    return (Complex){
        a.re * b.re - a.im * b.im,   /* real */
        a.re * b.im + a.im * b.re    /* imaginary */
    };
}
```

> **Key insight:** Each complex multiply requires 4 real multiplies and
> 2 real adds. This is why SIMD optimization (Chapter 07) targets this
> operation first.

## 1.4 The Complex Exponential

Euler's formula connects exponentials to trigonometry:

$$e^{i\theta} = \cos\theta + i\sin\theta$$

In the FFT, the "twiddle factor" for stage of size $N$ is:

$$W_N = e^{-j \cdot 2\pi / N} = \cos\!\left(\frac{2\pi}{N}\right) - i\sin\!\left(\frac{2\pi}{N}\right)$$

Look at how this appears in [`src/fft.c`](../src/fft.c) line 97:

```c
double angle = -2.0 * M_PI / stage_size;
Complex w_base = { cos(angle), sin(angle) };
```

This creates $W_N$ in rectangular form, ready for butterfly multiplies.

## 1.5 Try It Yourself

Build and run the FFT demo — it creates complex numbers, computes
magnitudes, and prints them in dB:

```bash
make
./build/bin/fft_demo
```

Look at the output and trace back to [`01-complex-numbers.c`](01-complex-numbers.c)
to see how `Complex`, `complex_mag`, and `db_from_magnitude` work together.

## 1.6 Exercises

1. **Paper exercise:** Multiply $(3 + 4i)$ by $(1 - 2i)$ by hand.
   Verify with `complex_mul`.

2. **Code exercise:** Write a small `main()` that creates
   $z = e^{i\pi/4}$ using `complex_from_polar(1.0, M_PI/4)` and prints
   the real and imaginary parts. What do you expect?

3. **Thinking question:** Why does the FFT use $e^{-j\theta}$ (negative
   sign) instead of $e^{+j\theta}$? *(Hint: convention for forward vs.
   inverse transform.)*

---

**Previous:** [Chapter 00 — Overview](00-overview.md)
| **Next:** [Chapter 02 — FFT Fundamentals →](02-fft-fundamentals.md)

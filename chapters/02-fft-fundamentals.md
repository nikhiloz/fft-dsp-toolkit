# Chapter 2 — FFT Fundamentals

The Fast Fourier Transform converts a time-domain signal into its
frequency components. It's the single most important algorithm in
digital signal processing.

---

## 2.1 From DFT to FFT

The **Discrete Fourier Transform** of an N-point signal $x[n]$ is:

$$X[k] = \sum_{n=0}^{N-1} x[n] \cdot e^{-j 2\pi kn / N}, \quad k = 0, 1, \ldots, N{-}1$$

Computing this directly requires $N^2$ complex multiplies. The **FFT**
exploits symmetry and periodicity to do it in $O(N \log N)$.

For $N = 1024$: DFT needs ~1,000,000 multiplies. FFT needs ~10,000.
That's a **100× speedup**.

## 2.2 The Cooley-Tukey Algorithm

Our implementation uses the **Radix-2 Decimation-In-Time (DIT)** variant,
the most common form. It requires $N$ to be a power of 2.

> **📊 FFT Processing Sequence** — [View full-size diagram →](../reference/diagrams/fft_sequence.png)

### The Idea

Split the DFT into two half-size DFTs — one for even-indexed samples,
one for odd-indexed samples:

$$X[k] = \underbrace{\sum_{m=0}^{N/2-1} x[2m] \cdot W_{N/2}^{mk}}_{\text{Even: } E[k]}
       + W_N^k \cdot \underbrace{\sum_{m=0}^{N/2-1} x[2m+1] \cdot W_{N/2}^{mk}}_{\text{Odd: } O[k]}$$

This gives us the **butterfly equation**:

$$X[k] = E[k] + W_N^k \cdot O[k]$$
$$X[k + N/2] = E[k] - W_N^k \cdot O[k]$$

Two outputs from one complex multiply. Repeat recursively → $O(N \log N)$.

## 2.3 Implementation Walk-Through

The full implementation is in [`src/fft.c`](../src/fft.c) (~185 lines).

### Step 1: Bit-Reversal Permutation

Before the butterfly stages, we must reorder the input so that the
output ends up in natural order. The reordering follows a bit-reversal
pattern.

For $N = 8$:

```
Index:    0  1  2  3  4  5  6  7
Binary:  000 001 010 011 100 101 110 111
Reverse: 000 100 010 110 001 101 011 111
Result:   0  4  2  6  1  5  3  7
```

See [`src/fft.c`](../src/fft.c) lines 43–59:

```c
static void bit_reverse_permute(Complex *x, int n) {
    int j = 0;
    for (int i = 0; i < n - 1; i++) {
        if (i < j) {
            Complex tmp = x[i];
            x[i] = x[j];
            x[j] = tmp;
        }
        int m = n >> 1;
        while (m >= 1 && j >= m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
}
```

This is an in-place algorithm — no extra memory needed. The `if (i < j)`
guard ensures each pair is swapped only once.

### Step 2: Butterfly Stages

The core of the FFT: $\log_2(N)$ stages, each with $N/2$ butterflies.

```
Butterfly diagram (one pair):

    a ──────┬──── a + W·b     ("top wing")
            ×  W
    b ──────┴──── a − W·b     ("bottom wing")
```

See [`src/fft.c`](../src/fft.c) lines 83–118 for the main `fft()` function:

```c
void fft(Complex *x, int n) {
    if (n <= 1) return;
    bit_reverse_permute(x, n);

    for (int stage_size = 2; stage_size <= n; stage_size <<= 1) {
        int half = stage_size >> 1;
        double angle = -2.0 * M_PI / stage_size;
        Complex w_base = { cos(angle), sin(angle) };

        for (int group = 0; group < n; group += stage_size) {
            Complex w = { 1.0, 0.0 };
            for (int k = 0; k < half; k++) {
                int top = group + k;
                int bot = group + k + half;

                Complex t = complex_mul(w, x[bot]);  /* twiddle */
                Complex u = x[top];

                x[top] = complex_add(u, t);  /* top wing */
                x[bot] = complex_sub(u, t);  /* bottom wing */

                w = complex_mul(w, w_base);  /* advance twiddle */
            }
        }
    }
}
```

Key observations:
- **In-place**: Output overwrites input — no allocation needed
- **Three nested loops**: stages → groups → butterflies within group
- **Twiddle accumulation**: `w = complex_mul(w, w_base)` avoids
  recomputing `cos`/`sin` for every butterfly

### Step 3: Extracting Results

After `fft()` completes, use the helper functions:

```c
/* Get magnitude of each frequency bin */
fft_magnitude(result, magnitudes, N);

/* Get phase angle of each bin */
fft_phase(result, phases, N);
```

See [`src/fft.c`](../src/fft.c) lines 173–185.

## 2.4 Inverse FFT

The IFFT recovers the original signal from its frequency representation.
The trick: conjugate, forward-FFT, conjugate again, scale by $1/N$.

$$x[n] = \frac{1}{N} \sum_{k=0}^{N-1} X[k] \cdot e^{+j 2\pi kn / N}$$

See [`src/fft.c`](../src/fft.c) lines 138–153:

```c
void ifft(Complex *x, int n) {
    /* Conjugate input */
    for (int i = 0; i < n; i++) x[i].im = -x[i].im;

    /* Forward FFT */
    fft(x, n);

    /* Conjugate again + scale */
    double scale = 1.0 / n;
    for (int i = 0; i < n; i++) {
        x[i].re *= scale;
        x[i].im = -x[i].im * scale;
    }
}
```

This elegantly reuses the forward FFT — no need for a separate
implementation.

## 2.5 Frequency Resolution

For a signal sampled at $f_s$ Hz with $N$ samples:

- **Frequency resolution**: $\Delta f = f_s / N$
- **Bin $k$ corresponds to**: $f = k \cdot f_s / N$
- **Nyquist frequency**: $f_s / 2$ (bin $N/2$)

Only bins $0$ to $N/2$ contain unique information for real-valued signals
(the upper half is a mirror image — the complex conjugate).

## 2.6 Try It Yourself

```bash
make
./build/bin/fft_demo
```

The demo generates a signal with 440 Hz and 1000 Hz components, applies
a Hann window, computes a 256-point FFT, and prints the magnitude
spectrum. You should see peaks near those two frequencies.

Trace the code in [`02-fft-fundamentals.c`](02-fft-fundamentals.c) to see
each step.

## 2.7 Testing

The FFT is verified by 6 tests in [`tests/test_fft.c`](../tests/test_fft.c):

| Test | What It Verifies |
|------|-----------------|
| DC component | Constant signal → all energy in bin 0 |
| Impulse → flat | Delta function has equal magnitude across all bins |
| Alternating → Nyquist | $[1,-1,1,-1,\ldots]$ → peak at bin $N/2$ |
| FFT ↔ IFFT roundtrip | Forward then inverse recovers original within $10^{-9}$ |
| Pure sine peaks | Known-frequency sine → correct bin |
| `fft_real` wrapper | Matches manual complex FFT |

```bash
make test
```

## 2.8 Exercises

1. **Modify the demo** to use 512 points instead of 256. What happens
   to the frequency resolution? Do the peaks sharpen?

2. **Add a third frequency** (e.g., 2500 Hz) to the demo signal. Can
   you see it in the output? Is it aliased?

3. **Paper exercise:** For an 8-point FFT, draw the complete butterfly
   diagram (3 stages, 4 butterflies each). Label each twiddle factor.

4. **Code exercise:** Create a test that verifies Parseval's theorem:
   $$\sum|x[n]|^2 = \frac{1}{N}\sum|X[k]|^2$$

---

**Previous:** [Chapter 01 — Complex Numbers](01-complex-numbers.md)
| **Next:** [Chapter 03 — Window Functions →](03-window-functions.md)

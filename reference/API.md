# FFT-DSP Toolkit: API Reference

Public API for the three implemented modules. All functions are C99,
thread-safe (no global state), and have zero external dependencies
beyond `<math.h>`.

> **📊 API Overview** — [View full-size API diagram →](diagrams/api_reference.png)
>
> **📊 Module Dependencies** — [View full-size diagram →](diagrams/modules.png)

**Dependency rule:** `fft.h` and `filter.h` both include `dsp_utils.h`.
They do not depend on each other.

---

## Module 1: dsp_utils.h — Core Utilities

**Header:** [`include/dsp_utils.h`](../include/dsp_utils.h)
| **Source:** [`src/dsp_utils.c`](../src/dsp_utils.c)
| **Tutorial:** [Chapter 01](../chapters/01-complex-numbers.md), [Chapter 03](../chapters/03-window-functions.md)

### Data Types

```c
typedef struct {
    double re;  /* Real part      */
    double im;  /* Imaginary part */
} Complex;

typedef double (*window_fn)(int n, int i);
```

`Complex` is the fundamental type used throughout the toolkit.
`window_fn` is a function pointer for any window function with the
signature `double f(int window_length, int sample_index)`.

### Complex Arithmetic

| Function | Signature | Description |
|----------|-----------|-------------|
| `complex_add` | `Complex complex_add(Complex a, Complex b)` | $(a + b)$ |
| `complex_sub` | `Complex complex_sub(Complex a, Complex b)` | $(a - b)$ |
| `complex_mul` | `Complex complex_mul(Complex a, Complex b)` | $(a \times b)$ using $(ac-bd) + (ad+bc)i$ |
| `complex_mag` | `double complex_mag(Complex z)` | $\|z\| = \sqrt{re^2 + im^2}$ |
| `complex_phase` | `double complex_phase(Complex z)` | $\text{atan2}(im, re)$ in radians |
| `complex_from_polar` | `Complex complex_from_polar(double mag, double phase)` | $(r, \theta) \to (a, b)$ |

### Window Functions

| Function | Signature | Side Lobes | Main Lobe |
|----------|-----------|-----------|-----------|
| `hann_window` | `double hann_window(int n, int i)` | −31 dB | 4 bins |
| `hamming_window` | `double hamming_window(int n, int i)` | −42 dB | 4 bins |
| `blackman_window` | `double blackman_window(int n, int i)` | −58 dB | 6 bins |

#### `void apply_window(double *signal, int n, window_fn w)`

Multiplies each sample in-place: `signal[i] *= w(n, i)`.

**Parameters:**
- `signal` — Array of `n` doubles (modified in-place)
- `n` — Number of samples
- `w` — Window function pointer (e.g., `hann_window`)

### Helpers

#### `int next_power_of_2(int n)`

Returns the smallest power of 2 ≥ `n`. Used to determine zero-padding
length for FFT input.

#### `double db_from_magnitude(double mag)`

Converts a magnitude value to decibels: $20 \cdot \log_{10}(\text{mag})$.
Returns −200 dB for zero or negative input (avoids `log(0)`).

#### `double rms(const double *signal, int n)`

Root-mean-square: $\sqrt{\frac{1}{N}\sum x[i]^2}$.

---

## Module 2: fft.h — Fast Fourier Transform

**Header:** [`include/fft.h`](../include/fft.h)
| **Source:** [`src/fft.c`](../src/fft.c)
| **Tutorial:** [Chapter 02](../chapters/02-fft-fundamentals.md)

**Algorithm:** Cooley-Tukey Radix-2 Decimation-In-Time.
**Complexity:** $O(N \log N)$.
**Constraint:** `n` must be a power of 2.

### Functions

#### `void fft(Complex *x, int n)`

In-place forward FFT. Transforms `x` from time domain to frequency
domain.

**Parameters:**
- `x` — Array of `n` complex values (overwritten with result)
- `n` — Transform size (must be power of 2)

**Example:**
```c
Complex x[8] = { {1,0}, {0,0}, {0,0}, {0,0},
                 {0,0}, {0,0}, {0,0}, {0,0} };
fft(x, 8);
/* x now contains the 8-point DFT of an impulse */
```

#### `void fft_real(const double *in, Complex *out, int n)`

Convenience wrapper: copies real-valued input into complex format,
then runs `fft()`.

**Parameters:**
- `in` — Real-valued input array [n]
- `out` — Complex output array [n] (receives FFT result)
- `n` — Transform size (power of 2)

#### `void ifft(Complex *x, int n)`

In-place inverse FFT. Recovers the time-domain signal.
Uses the conjugate trick: conjugate → fft → conjugate → scale by 1/N.

**Parameters:** Same as `fft()`.

#### `void fft_magnitude(const Complex *x, double *mag, int n)`

Extracts magnitude from each complex bin: `mag[i] = |x[i]|`.

#### `void fft_phase(const Complex *x, double *phase, int n)`

Extracts phase angle from each complex bin: `phase[i] = atan2(im, re)`.

---

## Module 3: filter.h — FIR Digital Filters

**Header:** [`include/filter.h`](../include/filter.h)
| **Source:** [`src/filter.c`](../src/filter.c)
| **Tutorial:** [Chapter 04](../chapters/04-digital-filters.md)

**Method:** Direct-form FIR convolution.
**Complexity:** $O(N \times M)$ where N = signal length, M = filter order.

### Functions

#### `void fir_filter(const double *in, double *out, int n, const double *h, int order)`

Applies an FIR filter using direct convolution.

**Parameters:**
- `in` — Input signal [n]
- `out` — Output signal [n] (may not alias `in`)
- `n` — Signal length
- `h` — Filter coefficients [order]
- `order` — Number of taps

**Behaviour:** For samples before the signal start (`i - k < 0`),
zero-padding is assumed.

**Example:**
```c
double in[100], out[100];
double h[5] = {0.2, 0.2, 0.2, 0.2, 0.2};  /* 5-tap average */

/* ... fill in[] ... */
fir_filter(in, out, 100, h, 5);
```

#### `void fir_moving_average(double *h, int taps)`

Generates moving-average coefficients: `h[i] = 1.0 / taps` for all i.

**Parameters:**
- `h` — Output coefficient array [taps]
- `taps` — Number of filter taps

#### `void fir_lowpass(double *h, int taps, double cutoff)`

Designs a lowpass FIR filter using the windowed-sinc method with a
Hamming window. Coefficients are normalised for unity DC gain.

**Parameters:**
- `h` — Output coefficient array [taps]
- `taps` — Number of filter taps (odd recommended)
- `cutoff` — Normalised cutoff frequency ($f_c / f_s$, range 0.0–0.5)

**Example:**
```c
#define TAPS 31
double h[TAPS];
double fs = 8000.0;

/* 500 Hz lowpass at 8000 Hz sample rate */
fir_lowpass(h, TAPS, 500.0 / fs);

/* Apply to signal */
fir_filter(input, output, signal_len, h, TAPS);
```

---

## Compilation & Linking

### Build with Make

```bash
make              # Debug build (-g -Wall -Wextra -Werror)
make release      # Optimised build (-O3 -DNDEBUG)
make test         # Build + run all tests
make clean        # Remove build artefacts
```

### Link Your Application

```bash
cc -Iinclude -o my_app my_app.c src/fft.c src/filter.c src/dsp_utils.c -lm
```

Or link against the built object files:

```bash
cc -Iinclude -o my_app my_app.c build/obj/fft.o build/obj/filter.o build/obj/dsp_utils.o -lm
```

### CMake

```cmake
add_subdirectory(fft-dsp-toolkit)
target_link_libraries(my_app PRIVATE fft_dsp_toolkit)
```

---

## See Also

- [ARCHITECTURE.md](ARCHITECTURE.md) — System design and module relationships
- [chapters/](../chapters/00-overview.md) — Progressive learning chapters
- [diagrams/](diagrams/) — All PlantUML diagrams

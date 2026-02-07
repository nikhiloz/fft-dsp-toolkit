# FFT-DSP Toolkit: API Reference

Complete public API reference for the FFT-DSP Toolkit library.

## Overview

The FFT-DSP Toolkit provides a modular C99 API organized into 8 core modules:

![Public API Reference](diagrams/api_reference.png)

![Module Dependencies](diagrams/modules.png)

- **fft.h** - Fast Fourier Transform operations
- **filter.h** - Digital filtering (FIR/IIR)
- **window.h** - Window functions for spectral analysis
- **spectrum.h** - Spectral analysis and feature extraction
- **convolution.h** - Convolution and correlation
- **ring_buffer.h** - Real-time circular buffers
- **signal_gen.h** - Signal generation and synthesis
- **dsp_utils.h** - Utility functions and complex arithmetic

---

## Module 1: fft.h - Fast Fourier Transform

### Overview
Core FFT implementation with support for real and complex signals, multiple precision levels, and streaming mode.

### Data Types

```c
typedef struct {
    double re, im;  // Real and imaginary parts
} complex_t;

typedef struct {
    int N;                  // Transform size
    double *twiddle;        // Pre-computed twiddle factors
    int *bit_reverse;       // Bit-reversal permutation
    double norm_factor;     // Normalization for IFFT
} fft_plan_t;
```

### Functions

#### `fft_plan_t* fft_init(int N)`
Initialize FFT plan with transform size N.
- **Parameters**: N - Transform size (power of 2)
- **Returns**: FFT plan pointer (must be freed with fft_free)
- **Time Complexity**: O(N log N) preprocessing
- **Space Complexity**: O(N) for twiddle factors

**Example**:
```c
fft_plan_t *plan = fft_init(1024);
// Use plan for FFT operations
fft_free(plan);
```

#### `void fft_transform(fft_plan_t *plan, complex_t *input, complex_t *output)`
Compute forward FFT.
- **Parameters**: 
  - plan: FFT plan (from fft_init)
  - input: Complex input array [N]
  - output: Complex output array [N]
- **In-place**: Can use input==output
- **Time Complexity**: O(N log N)
- **Notes**: Input assumed to be in time domain; output in frequency domain

#### `void ifft_transform(fft_plan_t *plan, complex_t *input, complex_t *output)`
Compute inverse FFT with automatic scaling.
- **Parameters**: Same as fft_transform
- **Notes**: Properly scaled (divide by N not required)
- **Time Complexity**: O(N log N)

#### `void fft_real(fft_plan_t *plan, double *input, complex_t *output)`
Real-valued FFT (optimized for real signals).
- **Parameters**:
  - input: Real-valued input [N]
  - output: Complex output [N/2+1] (only unique frequencies)
- **Performance**: ~2x faster than complex FFT
- **Memory**: Output size N/2+1 vs N for complex

#### `void fft_free(fft_plan_t *plan)`
Free FFT plan resources.
- **Parameters**: plan pointer
- **Behavior**: Invalidates plan; further use is undefined

### Example Usage

```c
#include "fft.h"

int main() {
    int N = 512;
    fft_plan_t *plan = fft_init(N);
    
    complex_t input[512], output[512];
    
    // Fill input
    for (int i = 0; i < N; i++) {
        input[i].re = sin(2 * M_PI * i / N);
        input[i].im = 0.0;
    }
    
    // Transform
    fft_transform(plan, input, output);
    
    // Use output (frequency domain)
    for (int k = 0; k < N; k++) {
        double mag = sqrt(output[k].re*output[k].re + 
                         output[k].im*output[k].im);
        printf("Bin %d: %.6f\n", k, mag);
    }
    
    fft_free(plan);
    return 0;
}
```

---

## Module 2: filter.h - Digital Filtering

### Overview
Time-domain FIR and IIR digital filtering for signal processing.

### Data Types

```c
typedef struct {
    double *coeffs;         // Filter coefficients
    double *state;          // Internal state (IIR only)
    int order;              // Filter order
    int type;               // FIR or IIR
} filter_t;

typedef enum {
    FILTER_FIR = 0,
    FILTER_IIR = 1
} filter_type_t;
```

### Functions

#### `filter_t* fir_init(double *coeffs, int order)`
Initialize FIR filter.
- **Parameters**:
  - coeffs: Filter coefficients [order]
  - order: Number of coefficients (M)
- **Returns**: Filter structure
- **Notes**: Coefficients are copied (not aliased)

#### `void fir_process(filter_t *filter, double *input, double *output, int N)`
Process signal through FIR filter.
- **Parameters**:
  - input: Input samples [N]
  - output: Output samples [N]
  - N: Number of samples
- **Output Delay**: order-1 samples latency
- **Complexity**: O(N * order)

#### `filter_t* iir_init(double *b, int b_len, double *a, int a_len)`
Initialize IIR filter.
- **Parameters**:
  - b: Feedforward coefficients (numerator)
  - a: Feedback coefficients (denominator)
  - b_len, a_len: Coefficient lengths
- **Note**: a[0] typically normalized to 1.0

#### `void iir_process(filter_t *filter, double *input, double *output, int N)`
Process signal through IIR filter.
- **Complexity**: O(N * max(a_len, b_len))
- **Stability**: Check pole locations before use

#### `void filter_reset(filter_t *filter)`
Reset filter state to zero.
- **Effect**: Clears internal buffers and history
- **Use Case**: Between discontinuous signal segments

#### `void filter_free(filter_t *filter)`
Release filter resources.

### Example Usage

```c
#include "filter.h"

int main() {
    // Simple FIR low-pass (3-tap moving average)
    double coeffs[] = {0.333, 0.334, 0.333};
    filter_t *fir = fir_init(coeffs, 3);
    
    double input[] = {1, 0, 1, 0, 1, 0};
    double output[6];
    
    fir_process(fir, input, output, 6);
    
    filter_free(fir);
    return 0;
}
```

---

## Module 3: window.h - Window Functions

### Overview
Spectral window functions for reducing spectral leakage in FFT analysis.

### Functions

#### `double* window_hann(int N)`
Hann (Hanning) window.
- **Sidelobe Level**: -32 dB
- **Main Lobe Width**: 4.0 bins
- **Use Case**: General-purpose spectral analysis

#### `double* window_hamming(int N)`
Hamming window.
- **Sidelobe Level**: -43 dB
- **Main Lobe Width**: 4.0 bins
- **Use Case**: Standard spectral analysis

#### `double* window_blackman(int N)`
Blackman window (3-term).
- **Sidelobe Level**: -58 dB
- **Main Lobe Width**: 6.0 bins
- **Use Case**: High dynamic range applications

#### `double* window_kaiser(int N, double beta)`
Kaiser window with adjustable shape.
- **Parameters**: 
  - beta: Shape parameter (0-20+)
  - beta ~ 5: Near Hamming response
  - beta ~ 8: Near Blackman response
- **Use Case**: Design for specific sidelobe requirements

### Properties Table

| Window | Main Lobe | Sidelobe | Rolloff | Processing |
|--------|-----------|----------|--------|------------|
| Rectangular | 1.0 | -13 dB | -6 dB/oct | None |
| Hann | 4.0 | -32 dB | -18 dB/oct | Moderate |
| Hamming | 4.0 | -43 dB | -6 dB/oct | Moderate |
| Blackman | 6.0 | -58 dB | -18 dB/oct | High |
| Kaiser (β=8.6) | 5.0 | -50 dB | -6 dB/oct | Tunable |

### Example Usage

```c
#include "window.h"
#include "fft.h"

int main() {
    int N = 1024;
    double *window = window_hann(N);
    
    double *signal = malloc(N * sizeof(double));
    // ... fill signal ...
    
    // Apply window
    for (int i = 0; i < N; i++) {
        signal[i] *= window[i];
    }
    
    free(window);
    free(signal);
    return 0;
}
```

---

## Module 4: spectrum.h - Spectral Analysis

### Overview
Frequency domain analysis tools (PSD, peak detection, spectral features).

### Functions

#### `double* spectrum_psd(double *signal, int N, double sr)`
Compute power spectral density using periodogram.
- **Parameters**:
  - signal: Time-domain signal [N]
  - N: Signal length
  - sr: Sample rate (Hz)
- **Returns**: PSD [N/2] in dBFS
- **Algorithm**: Welch's method with default settings

#### `double* spectrum_welch(double *signal, int N, double sr, int n_segments)`
Compute PSD using Welch's method (lower variance).
- **Parameters**: 
  - n_segments: Number of overlapping windows
- **Returns**: Averaged PSD
- **Overlap**: 50% (standard)
- **Variance**: Reduced by ~n_segments

#### `int* spectrum_peaks(double *magnitude, int N, double threshold)`
Find spectral peaks above threshold.
- **Parameters**:
  - magnitude: Magnitude spectrum [N]
  - threshold: Minimum magnitude (0.0-1.0)
- **Returns**: Array of peak bin indices (NULL-terminated)

#### `double* spectrum_freqs(double sr, int N)`
Generate frequency axis for spectrum display.
- **Returns**: Frequency array [N/2] in Hz
- **Spacing**: sr/N Hz between bins

### Example Usage

```c
#include "spectrum.h"

int main() {
    int N = 1024;
    double sr = 44100.0;  // 44.1 kHz sample rate
    
    double *signal = malloc(N * sizeof(double));
    // ... fill signal ...
    
    double *psd = spectrum_psd(signal, N, sr);
    int *peaks = spectrum_peaks(psd, N/2, 0.1);
    
    printf("Peaks found at bins: ");
    for (int i = 0; peaks[i] != -1; i++) {
        double freq = peaks[i] * sr / N;
        printf("%d (%.1f Hz) ", peaks[i], freq);
    }
    printf("\n");
    
    free(psd);
    free(peaks);
    free(signal);
    return 0;
}
```

---

## Module 5: convolution.h - Convolution & Correlation

### Overview
Linear and circular convolution, auto/cross-correlation.

### Functions

#### `double* convolve(double *x, int x_len, double *h, int h_len)`
Linear convolution (zero-padded).
- **Parameters**:
  - x: Signal [x_len]
  - h: Filter/kernel [h_len]
- **Returns**: Output [x_len + h_len - 1]
- **Algorithm**: FFT-based for efficiency
- **Complexity**: O((x_len + h_len) log(x_len + h_len))

#### `double* convolve_circular(double *x, int N, double *h, int M)`
Circular convolution (periodic).
- **Complexity**: O(N log N)
- **Use Case**: Frequency-domain filtering

#### `double* correlate_auto(double *x, int N)`
Autocorrelation (lag correlation with self).
- **Returns**: Autocorrelation [N]
- **Peak**: At lag 0
- **Use Case**: Periodicity detection

#### `double* correlate_cross(double *x, int x_len, double *y, int y_len)`
Cross-correlation between two signals.
- **Returns**: Correlation [x_len + y_len - 1]
- **Use Case**: Signal alignment, matching

---

## Module 6: ring_buffer.h - Circular Buffers

### Overview
Lock-free ring buffers for real-time streaming applications.

### Functions

#### `ring_buffer_t* rbuf_create(int capacity)`
Create circular ring buffer.
- **Parameters**: capacity - Buffer size (power of 2 recommended)
- **Returns**: Ring buffer handle
- **Thread-Safe**: Yes (lock-free)

#### `int rbuf_write(ring_buffer_t *buf, double *data, int count)`
Write samples to buffer (producer).
- **Returns**: Number of samples actually written
- **Behavior**: Non-blocking; may write fewer samples if full

#### `int rbuf_read(ring_buffer_t *buf, double *data, int count)`
Read samples from buffer (consumer).
- **Returns**: Number of samples read
- **Behavior**: Non-blocking; may read fewer if empty

#### `int rbuf_capacity(ring_buffer_t *buf)`
Get current buffer utilization.
- **Returns**: Number of samples in buffer

#### `void rbuf_destroy(ring_buffer_t *buf)`
Free ring buffer resources.

### Example Usage

```c
#include "ring_buffer.h"

int main() {
    ring_buffer_t *buf = rbuf_create(4096);
    
    double samples_in[512];
    double samples_out[512];
    
    // Write to buffer
    rbuf_write(buf, samples_in, 512);
    
    // Read when ready
    int read = rbuf_read(buf, samples_out, 512);
    printf("Read %d samples\n", read);
    
    rbuf_destroy(buf);
    return 0;
}
```

---

## Module 7: signal_gen.h - Signal Generation

### Overview
Synthesize test and reference signals.

### Functions

#### `double* gen_sine(double freq, double duration, double sr)`
Generate sine wave.
- **Duration**: In seconds
- **sr**: Sample rate (Hz)
- **Returns**: Signal array [duration * sr]

#### `double* gen_chirp(double f0, double f1, double duration, double sr)`
Generate linear frequency sweep (chirp).
- **f0, f1**: Start and end frequencies (Hz)
- **Sweep Type**: Linear frequency increase

#### `double* gen_noise_white(int N)`
White Gaussian noise.
- **Power**: Normalized to unit variance
- **Returns**: Noise array [N]

#### `double* gen_noise_pink(int N)`
Pink 1/f noise (1/f^0.5 power spectrum).
- **Use Case**: Natural signal simulation

---

## Module 8: dsp_utils.h - Utility Functions

### Overview
Core mathematical operations and complex number arithmetic.

### Functions

#### `complex_t complex_add(complex_t a, complex_t b)`
Complex number addition.

#### `complex_t complex_mult(complex_t a, complex_t b)`
Complex number multiplication.

#### `complex_t complex_conj(complex_t a)`
Complex conjugate.

#### `double magnitude(complex_t c)`
Magnitude/absolute value.

#### `double phase(complex_t c)`
Phase angle (-π to +π).

#### `double db_magnitude(complex_t c)`
Magnitude in dB (20*log10(|c|)).

### Example Usage

```c
#include "dsp_utils.h"
#include <math.h>

int main() {
    complex_t a = {3.0, 4.0};
    double mag = magnitude(a);  // 5.0
    double phase_val = phase(a);  // atan2(4, 3) 
    double db = db_magnitude(a);  // 20*log10(5) ≈ 13.98
    
    return 0;
}
```

---

## Compilation & Linking

### Build with Make
```bash
make release                 # Build optimized library
make debug                   # Build with debug symbols
make install                 # Install to /usr/local
```

### Link Your Application
```bash
gcc -Iinclude myapp.c -o myapp -Lbuild/lib -lfft_dsp -lm
```

### CMake
```cmake
find_package(FFT_DSP REQUIRED)
target_link_libraries(myapp FFT_DSP::fft_dsp)
```

---

## Error Handling

All functions return NULL pointers or -1 on error. Check return values:

```c
fft_plan_t *plan = fft_init(512);
if (plan == NULL) {
    perror("FFT initialization failed");
    return -1;
}
```

---

## Performance Tips

1. **FFT Size**: Use power-of-2 sizes for best performance
2. **Windowing**: Always window before FFT to reduce spectral leakage
3. **Streaming**: Use ring buffers with overlap-add for real-time
4. **SIMD**: Compile with `-O3 -march=native` to enable SIMD
5. **Memory**: Pre-allocate and reuse buffers in loops

---

## See Also

- [ARCHITECTURE.md](ARCHITECTURE.md) - System design
- [PERFORMANCE.md](PERFORMANCE.md) - Optimization guide
- [REALTIME.md](REALTIME.md) - Real-time programming
- [docs/examples/](../examples/) - Complete code examples
- [docs/diagrams/](diagrams/) - Architecture diagrams

---

**Version**: 1.0 | **Last Updated**: 2026-02-07 | **License**: MIT

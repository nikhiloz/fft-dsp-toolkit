# Chapter 8 — Putting It Together

This final chapter walks through two complete end-to-end examples that
combine everything from Chapters 01–07.

---

## 8.1 Example 1: Spectral Analysis of a Dual-Tone Signal

**Goal:** Detect the frequencies present in a signal.

**Chapters used:** Complex numbers (01), FFT (02), Windows (03),
Spectral analysis (05).

### The Code

Full source: [`08-putting-it-together.c`](08-putting-it-together.c).

### Step-by-Step Trace

| Step | Code | Concept | Chapter |
|------|------|---------|---------|
| 1. Define signal | `sin(2π·440·t) + sin(2π·1000·t)` | Signal generation | — |
| 2. Apply window | `apply_window(signal, N, hann_window)` | Reduce leakage | 03 |
| 3. FFT | `fft_real(signal, spectrum, N)` | Time → frequency | 02 |
| 4. Magnitude | `fft_magnitude(spectrum, mag, N)` | Complex → real | 01 |
| 5. dB scale | `db_from_magnitude(mag[k])` | Human-readable | 05 |
| 6. Frequency axis | `k * FS / N` | Resolution | 05 |

### Run It

```bash
make && ./build/bin/fft_demo
```

### What to Observe

- Peaks appear near 437.5 Hz and 1000 Hz
- The 440 Hz peak splits across 2 bins (it doesn't fall exactly on a
  bin boundary)
- Side lobes are suppressed below −25 dB thanks to the Hann window
- Only bins 0–128 (0 to 4000 Hz) are unique for real input

## 8.2 Example 2: Noise Reduction with FIR Filtering

**Goal:** Remove high-frequency noise from a signal while preserving
the desired low-frequency content.

**Chapters used:** Complex numbers (01), Windows (03), Filters (04).

### The Code

Full source: [`08-putting-it-together.c`](08-putting-it-together.c) (combines both demos).

### Step-by-Step Trace

| Step | Code | Concept | Chapter |
|------|------|---------|---------|
| 1. Clean signal | `sin(2π·200·t)` | 200 Hz tone | — |
| 2. Add noise | `+= 0.3 · sin(2π·3500·t) + ...` | High-freq contamination | — |
| 3. Design filter | `fir_lowpass(h, 31, 500/FS)` | Windowed-sinc design | 04, 03 |
| 4. Filter | `fir_filter(noisy, filtered, N, h, TAPS)` | Direct convolution | 04 |
| 5. Measure | `rms(filtered + TAPS, N - TAPS)` | Skip transient | 05 |

### Run It

```bash
make && ./build/bin/filter_demo
```

### What to Observe

- **Clean RMS** ~0.71 (pure sine)
- **Noisy RMS** ~0.76 (added ~6% noise energy)
- **Filtered RMS** ~0.66 (noise removed, slight amplitude loss
  during transient settling)
- The filter has 31 taps → 15-sample group delay → first 31 outputs
  are ignored as transient

## 8.3 How the Modules Connect

> **📊 Module Dependencies** — [View full-size diagram →](../reference/diagrams/modules.png)

```
         ┌──────────────┐
         │  dsp_utils.h  │  ← Complex type, windows, helpers
         └──────┬───────┘
                │
        ┌───────┴───────┐
        ↓               ↓
   ┌─────────┐    ┌──────────┐
   │  fft.h   │    │ filter.h  │
   └────┬────┘    └─────┬────┘
        │               │
        ↓               ↓
   fft_demo.c     filter_demo.c
```

- `dsp_utils` is the foundation — everyone includes it
- `fft` and `filter` are independent — they don't include each other
- The demos include everything they need via the headers

## 8.4 Testing the Full Stack

All 12 tests can be run with:

```bash
make test
```

The test structure mirrors the module structure:

| Test File | Module Tested | Tests |
|-----------|--------------|------:|
| [`tests/test_fft.c`](../tests/test_fft.c) | `fft.c` + `dsp_utils.c` | 6 |
| [`tests/test_filter.c`](../tests/test_filter.c) | `filter.c` + `dsp_utils.c` | 6 |

Each test:
1. Creates a known input signal
2. Processes it through the module under test
3. Checks the output against expected mathematical properties
4. Reports PASS/FAIL

## 8.5 Where to Go from Here

### Extend the Current Code

- Add more window functions (Kaiser, flat-top, Gaussian)
- Implement `fir_highpass()` and `fir_bandpass()`
- Add power spectral density estimation
- Implement auto-correlation and cross-correlation

### Build Phase 3: Real-Time Streaming

- Implement ring buffer (Chapter 06)
- Add ALSA audio input
- Create a live spectrum analyser

### Build Phase 4: Optimisation

- Radix-4 FFT (Chapter 07)
- SIMD kernels for x86 and ARM
- Benchmark against FFTW

> **📊 Development Roadmap** — [View full-size diagram →](../reference/diagrams/roadmap.png)

## 8.6 Final Exercises

1. **Integration project:** Combine the FFT and filter: generate a
   noisy signal, filter it, then FFT both the noisy and filtered
   versions. Compare the spectra side by side.

2. **Build a tuner:** Create a program that reads PCM audio from stdin,
   computes the FFT, and prints the dominant frequency continuously.

3. **Spectrogram:** Process a long signal in overlapping frames (50%
   overlap, Hann window). Print the magnitude of each frame as one
   row — this creates a time-frequency representation.

4. **Benchmark:** Time `fft()` for powers of 2 from 64 to 65536.
   Plot the results and verify $O(N \log N)$ scaling.

5. **Cross-platform:** Build with CMake instead of Make. Verify
   identical test results on a different compiler or OS.

---

**Previous:** [Chapter 07 — Optimisation](07-optimisation.md)
| **Back to start:** [Chapter 00 — Overview](00-overview.md)

---

*Congratulations — you've completed the FFT-DSP tutorial. You now
understand complex arithmetic, the Cooley-Tukey FFT, window functions,
FIR filter design, spectral analysis, and have a foundation for
real-time processing and optimisation.*

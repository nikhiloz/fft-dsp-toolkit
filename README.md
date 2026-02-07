# DSP Tutorial Suite

**A comprehensive C tutorial for Digital Signal Processing — from first principles to postgraduate topics.**

This repository is a progressive learning resource covering 30 chapters of DSP, from
discrete-time signals through adaptive filters. Every source file is written to be
*read*, with detailed comments that cross-reference tutorial chapters, textbook
sections, and diagrams. Zero external dependencies — just C99 and `math.h`.

---

## What You'll Learn

### Part I — Foundations

| Ch | Topic | Demo | Key Files |
|----|-------|------|-----------|
| [01](chapters/01-signals-and-sequences.md) | Discrete-time signals & sequences | `ch01` | [`signal_gen.h`](include/signal_gen.h) |
| [02](chapters/02-sampling-and-aliasing.md) | Sampling, aliasing & Nyquist theorem | `ch02` | — |
| [03](chapters/03-complex-numbers.md) | Complex numbers & Euler's formula | `ch03` | [`dsp_utils.h`](include/dsp_utils.h) |
| [04](chapters/04-lti-systems.md) | LTI systems & discrete convolution | `ch04` | [`convolution.h`](include/convolution.h) |

### Part II — Transform Domain

| Ch | Topic | Demo | Key Files |
|----|-------|------|-----------|
| [05](chapters/05-z-transform.md) | The Z-Transform | `ch05` | — |
| [06](chapters/06-frequency-response.md) | Frequency response, poles & zeros | `ch06` | [`iir.h`](include/iir.h) |
| [07](chapters/07-dft-theory.md) | The DFT — theory & properties | `ch07` | — |
| [08](chapters/08-fft-fundamentals.md) | FFT algorithms (Cooley-Tukey Radix-2) | `ch08` | [`fft.h`](include/fft.h), [`fft.c`](src/fft.c) |
| [09](chapters/09-window-functions.md) | Window functions & spectral leakage | `ch09` | [`dsp_utils.c`](src/dsp_utils.c) |

### Part III — Filter Design

| Ch | Topic | Demo | Key Files |
|----|-------|------|-----------|
| [10](chapters/10-digital-filters.md) | FIR filter design | `ch10` | [`filter.h`](include/filter.h), [`filter.c`](src/filter.c) |
| [11](chapters/11-iir-filter-design.md) | IIR filter design (Butterworth, Chebyshev) | `ch11` | [`iir.h`](include/iir.h), [`iir.c`](src/iir.c) |
| [12](chapters/12-filter-structures.md) | Filter structures (biquads, SOS cascades) | `ch12` | — |

### Part IV — Analysis

| Ch | Topic | Demo | Key Files |
|----|-------|------|-----------|
| [13](chapters/13-spectral-analysis.md) | Spectral analysis | `ch13` | [`spectrum.h`](include/spectrum.h) |
| [14](chapters/14-psd-welch.md) | Power spectral density (Welch's method) | `ch14` | [`spectrum.c`](src/spectrum.c) |
| [15](chapters/15-correlation.md) | Correlation & autocorrelation | `ch15` | [`correlation.h`](include/correlation.h) |

### Part V — C-Specific & Advanced UG

| Ch | Topic | Demo | Key Files |
|----|-------|------|-----------|
| [16](chapters/16-overlap-add-save.md) | Overlap-Add/Save streaming convolution | `ch16` | [`streaming.h`](include/streaming.h) |
| [17](chapters/17-multirate-dsp.md) | Multirate DSP (decimation, interpolation, polyphase) | `ch17` | [`multirate.h`](include/multirate.h) |
| [18](chapters/18-fixed-point.md) | Fixed-point arithmetic (Q15/Q31, SQNR) | `ch18` | [`fixed_point.h`](include/fixed_point.h) |
| [19](chapters/19-advanced-fft.md) | Advanced FFT (Goertzel, DTMF, Sliding DFT) | `ch19` | [`advanced_fft.h`](include/advanced_fft.h) |
| [20](chapters/20-hilbert-transform.md) | Quadrature signals & Hilbert transform | `ch20` | [`hilbert.h`](include/hilbert.h) |
| [21](chapters/21-signal-averaging.md) | Signal averaging & noise reduction | `ch21` | [`averaging.h`](include/averaging.h) |
| [22](chapters/22-advanced-fir.md) | Advanced FIR design (Parks-McClellan / IRLS) | `ch22` | [`remez.h`](include/remez.h) |

### Part VI — Postgraduate

| Ch | Topic | Status |
|----|-------|--------|
| 23 | Adaptive filters (LMS / RLS) | 🔜 |
| 24 | Linear prediction & parametric modelling | 🔜 |
| 25 | Parametric spectral estimation (MUSIC, ESPRIT) | 🔜 |
| 26 | Cepstrum analysis & MFCCs | 🔜 |

### Part VII — Applied / Capstone

| Ch | Topic | Demo | Status |
|----|-------|------|--------|
| 27 | 2D DSP & image processing | — | 🔜 |
| [28](chapters/28-real-time-streaming.md) | Real-time system design | — | 📋 design |
| [29](chapters/29-optimisation.md) | SIMD & hardware optimisation | — | 📋 design |
| [30](chapters/30-putting-it-together.md) | End-to-end capstone project | `ch30` | ✅ |

## Quick Start

```bash
# Clone
git clone git@github.com:nikhiloz/dsp-tutorial-suite.git
cd dsp-tutorial-suite

# Build everything (C99, no external deps)
make release

# Run a specific chapter demo
./build/bin/ch01    # Signals & sequences
./build/bin/ch08    # FFT fundamentals
./build/bin/ch18    # Fixed-point arithmetic

# Run the test suite (61 tests across 6 suites)
make test

# Run all chapter demos
make run

# Generate all gnuplot visualisations
make plots
```

### Requirements

- GCC or Clang with C99 support
- GNU Make
- Linux / macOS (POSIX)
- *Optional*: gnuplot >= 5.0 for plot generation (`apt install gnuplot`)
- *Optional*: Java 11+ and [PlantUML](https://plantuml.com) to regenerate diagrams

## Project Layout

```
dsp-tutorial-suite/
├── include/          ← Public headers (start reading here)
│   ├── dsp_utils.h       Complex type, windows, helpers
│   ├── fft.h             FFT / IFFT API
│   ├── filter.h          FIR filter API
│   ├── iir.h             IIR filter design (Butterworth, Chebyshev)
│   ├── signal_gen.h      Signal generation (sine, noise, chirp)
│   ├── convolution.h     Convolution & correlation
│   ├── spectrum.h        PSD, Welch's method, cross-PSD
│   ├── correlation.h     Cross/auto-correlation (FFT-based)
│   ├── gnuplot.h         Pipe-based gnuplot plotting helpers
│   ├── fixed_point.h     Q15/Q31 fixed-point arithmetic
│   ├── advanced_fft.h    Goertzel, DTMF detection, sliding DFT
│   ├── streaming.h       Overlap-Add/Save block convolution
│   ├── multirate.h       Decimation, interpolation, polyphase
│   ├── hilbert.h         Hilbert transform, analytic signal
│   ├── averaging.h       Coherent averaging, EMA, median filter
│   └── remez.h           Parks-McClellan / IRLS equiripple FIR
├── src/              ← Reusable library (builds to libdsp_core.a)
│   ├── fft.c             Cooley-Tukey Radix-2 DIT
│   ├── filter.c          Direct convolution + sinc design
│   ├── dsp_utils.c       Complex arithmetic, 3 window functions
│   ├── signal_gen.c      Signal generators (Box-Muller noise, chirp)
│   ├── convolution.c     Linear/causal conv, cross/auto-correlation
│   ├── iir.c             IIR design + SOS cascade processing
│   ├── spectrum.c        Periodogram, Welch PSD, cross-PSD
│   ├── correlation.c     FFT-based xcorr/autocorr
│   ├── gnuplot.c         Gnuplot pipe helpers (pngcairo)
│   ├── fixed_point.c     Q15/Q31 conversion, saturating ops, FIR
│   ├── advanced_fft.c    Goertzel, DTMF, sliding DFT
│   ├── streaming.c       OLA/OLS block FFT convolution
│   ├── multirate.c       Multirate processing & polyphase
│   ├── hilbert.c         Hilbert FIR design, envelope, inst freq
│   ├── averaging.c       Coherent avg, EMA, MA, median filter
│   └── remez.c           IRLS-based equiripple FIR design
├── tests/            ← Unit tests (zero-dependency framework)
│   ├── test_framework.h  Lightweight test macros
│   ├── test_fft.c        6 FFT tests
│   ├── test_filter.c     6 FIR filter tests
│   ├── test_iir.c        10 IIR filter tests
│   ├── test_spectrum_corr.c  12 spectrum & correlation tests
│   ├── test_phase4.c     12 fixed-point, Goertzel, streaming tests
│   └── test_phase5.c     15 multirate, Hilbert, averaging, Remez tests
├── chapters/         ← START HERE — progressive tutorial (30 chapters)
│   ├── NN-topic.c        Demo code with ASCII art & rich comments
│   └── NN-topic.md       Theory tutorial with equations & exercises
├── tools/            ← Utilities
│   └── generate_plots.c  Generates 60+ gnuplot PNGs for all chapters
├── plots/            ← Generated visualisations (by chapter)
├── reference/        ← Architecture, API reference, diagrams
│   ├── ARCHITECTURE.md
│   ├── API.md
│   └── diagrams/     PlantUML sources + rendered PNGs
├── Makefile          ← Primary build (30+ targets)
└── CMakeLists.txt    ← Cross-platform alternative
```

## Diagrams

All architecture diagrams live in [`reference/diagrams/`](reference/diagrams/). They are
rendered from PlantUML source as high-resolution PNGs. Click any link below to
view the full-size image:

| Diagram | Description |
|---------|-------------|
| [System Architecture](reference/diagrams/architecture.png) | Four-layer system design |
| [Signal Flow](reference/diagrams/signal_flow.png) | Time → Frequency domain pipeline |
| [Module Dependencies](reference/diagrams/modules.png) | Source file dependency graph |
| [FFT Sequence](reference/diagrams/fft_sequence.png) | Runtime call sequence |
| [Real-Time Architecture](reference/diagrams/realtime_architecture.png) | Streaming pipeline |
| [Optimisation Roadmap](reference/diagrams/optimization_roadmap.png) | 5-stage speedup plan |
| [API Reference](reference/diagrams/api_reference.png) | Public function map |
| [Project Roadmap](reference/diagrams/roadmap.png) | 6-phase development plan |
| [Benchmarks](reference/diagrams/benchmarks.png) | Latency comparison |
| [Use Cases](reference/diagrams/use_cases.png) | Target applications |

To regenerate PNGs after editing `.puml` files:

```bash
java -jar ~/tools/plantuml.jar -tpng reference/diagrams/*.puml
```

## Test Output (61 tests)

```
=== Test Suite: FFT Functions ===
  Total: 6, Passed: 6, Failed: 0   (100%)

=== Test Suite: Filter Functions ===
  Total: 6, Passed: 6, Failed: 0   (100%)

=== Test Suite: IIR Filter Functions ===
  Total: 10, Passed: 10, Failed: 0 (100%)

=== Test Suite: Spectrum & Correlation ===
  Total: 12, Passed: 12, Failed: 0 (100%)

=== Test Suite: Phase 4: Fixed-Point, Advanced FFT, Streaming ===
  Total: 12, Passed: 12, Failed: 0 (100%)

=== Test Suite: Phase 5: Multirate, Hilbert, Averaging, Remez ===
  Total: 15, Passed: 15, Failed: 0 (100%)
```

## License

MIT — see [LICENSE](LICENSE).

## References

- Oppenheim & Willsky, *Signals and Systems* (3rd ed.)
- Oppenheim & Schafer, *Discrete-Time Signal Processing* (3rd ed.)
- Proakis & Manolakis, *Digital Signal Processing* (4th ed.)
- Haykin, *Adaptive Filter Theory* (5th ed.)
- Lyons, *Understanding DSP* (3rd ed.)
- Smith, *The Scientist and Engineer's Guide to DSP* (free online)

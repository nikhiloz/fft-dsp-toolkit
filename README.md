# FFT-DSP Toolkit

**A hands-on C tutorial for Digital Signal Processing — from complex numbers to real-time FFT.**

This repository is a progressive learning resource. Every source file is
written to be *read*, with detailed comments that cross-reference the
tutorial chapters and diagrams. You can study the theory, look at the
diagram, then open the matching `.c` file and see exactly how the math
becomes code.

---

## What You'll Learn

| Chapter | Topic | Key Source Files |
|---------|-------|-----------------|
| [00 – Overview](chapters/00-overview.md) | Project structure, build system, how to navigate | — |
| [01 – Complex Numbers](chapters/01-complex-numbers.md) | The building block of all frequency-domain work | [`dsp_utils.h`](include/dsp_utils.h), [`dsp_utils.c`](src/dsp_utils.c) |
| [02 – FFT Fundamentals](chapters/02-fft-fundamentals.md) | Cooley-Tukey Radix-2 DIT, butterfly operations | [`fft.h`](include/fft.h), [`fft.c`](src/fft.c) |
| [03 – Window Functions](chapters/03-window-functions.md) | Why windowing matters, Hann / Hamming / Blackman | [`dsp_utils.c`](src/dsp_utils.c) |
| [04 – Digital Filters](chapters/04-digital-filters.md) | FIR convolution, moving average, lowpass design | [`filter.h`](include/filter.h), [`filter.c`](src/filter.c) |
| [05 – Spectral Analysis](chapters/05-spectral-analysis.md) | Putting FFT + windows together to analyse signals | [`05-spectral-analysis.c`](chapters/05-spectral-analysis.c) |
| [06 – Real-Time Streaming](chapters/06-real-time-streaming.md) | Overlap-add, ring buffers, latency budgets | *(planned)* |
| [07 – Optimisation](chapters/07-optimisation.md) | SIMD, cache layout, multithreading | *(planned)* |
| [08 – Putting It Together](chapters/08-putting-it-together.md) | End-to-end project walkthrough | [`08-putting-it-together.c`](chapters/08-putting-it-together.c) |

## Quick Start

```bash
# Clone
git clone git@github.com:nikhiloz/fft-dsp-toolkit.git
cd fft-dsp-toolkit

# Build everything (C99, no external deps)
make

# Run the demos
./build/bin/fft_demo
./build/bin/filter_demo

# Run the test suite (12 tests)
make test
```

### Requirements

- GCC or Clang with C99 support
- GNU Make
- Linux / macOS (POSIX)
- *Optional*: Java 11+ and [PlantUML](https://plantuml.com) to regenerate diagrams

## Project Layout

```
fft-dsp-toolkit/
├── include/          ← Public headers (start reading here)
│   ├── dsp_utils.h       Complex type, windows, helpers
│   ├── fft.h             FFT / IFFT API
│   └── filter.h          FIR filter API
├── src/              ← Implementations (heavily commented)
│   ├── dsp_utils.c
│   ├── fft.c             Cooley-Tukey Radix-2 DIT
│   └── filter.c          Direct convolution + sinc design
├── tests/            ← Unit tests (zero-dependency framework)
│   ├── test_framework.h
│   ├── test_fft.c        6 FFT tests
│   └── test_filter.c     6 FIR filter tests
├── chapters/         ← START HERE — progressive tutorial
├── reference/        ← Architecture, API reference, diagrams
│   ├── ARCHITECTURE.md
│   ├── API.md
│   └── diagrams/     PlantUML sources + rendered PNGs
├── Makefile          ← Primary build (15+ targets)
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

## Test Output

```
=== Test Suite: FFT Functions ===
  [TEST] DC component (constant signal) .............. PASS
  [TEST] Impulse gives flat spectrum .................. PASS
  [TEST] Alternating signal → Nyquist bin ............. PASS
  [TEST] FFT then IFFT recovers original ............. PASS
  [TEST] Pure sine wave: peaks at correct bin ......... PASS
  [TEST] fft_real matches manual complex FFT .......... PASS
Total: 6, Passed: 6, Failed: 0

=== Test Suite: Filter Functions ===
  [TEST] Identity filter (1-tap passthrough) .......... PASS
  [TEST] Zero input gives zero output ................. PASS
  [TEST] Impulse response matches coefficients ........ PASS
  [TEST] Moving average smooths step input ............ PASS
  [TEST] Lowpass filter coefficients sum to 1.0 ....... PASS
  [TEST] Lowpass attenuates high frequency ............ PASS
Total: 6, Passed: 6, Failed: 0
```

## License

MIT — see [LICENSE](LICENSE).

## References

- Oppenheim & Schafer, *Discrete-Time Signal Processing*
- Haykin, *Signals and Systems*
- Proakis & Manolakis, *Digital Signal Processing*
- Smith, *The Scientist and Engineer's Guide to DSP* (free online)

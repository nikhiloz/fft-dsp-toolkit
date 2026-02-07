# Chapter 0 — Project Overview

Welcome to the **DSP Tutorial Suite**. This is a from-scratch
Digital Signal Processing library written in C99 with no external
dependencies. Every source file is designed to be *read* as a learning
resource.

---

## How to Navigate This Tutorial

Each chapter follows a consistent pattern:

1. **Theory** — The math and intuition behind the concept
2. **Diagram** — A visual reference (click to view full-size)
3. **Implementation Walk-Through** — Line-by-line analysis of the C code
4. **Try It Yourself** — Commands to build and run the related demo
5. **Exercises** — Practice problems to deepen understanding
6. **Next Chapter** — Link to continue the learning path

## Reading Order

| # | Chapter | Prerequisite Knowledge |
|---|---------|----------------------|
| **01** | [Complex Numbers](01-complex-numbers.md) | Basic algebra |
| **02** | [FFT Fundamentals](02-fft-fundamentals.md) | Chapter 01 |
| **03** | [Window Functions](03-window-functions.md) | Chapters 01–02 |
| **04** | [Digital Filters](04-digital-filters.md) | Chapter 01 |
| **05** | [Spectral Analysis](05-spectral-analysis.md) | Chapters 02–03 |
| **06** | [Real-Time Streaming](06-real-time-streaming.md) | Chapters 02–04 |
| **07** | [Optimisation](07-optimisation.md) | All previous |
| **08** | [Putting It Together](08-putting-it-together.md) | All previous |

## Project Structure

```
dsp-tutorial-suite/
├── include/            ← Start reading the PUBLIC API here
│   ├── dsp_utils.h         Complex type + windows + helpers
│   ├── fft.h               FFT/IFFT API
│   └── filter.h            FIR filter API
│
├── src/                ← Then dig into the IMPLEMENTATIONS
│   ├── dsp_utils.c         Complex arithmetic, 3 window functions
│   ├── fft.c               Cooley-Tukey Radix-2 DIT (~185 lines)
│   └── filter.c            FIR convolution + windowed-sinc design
│
├── tests/              ← Verify correctness
│   ├── test_framework.h    Zero-dependency test macros
│   ├── test_fft.c          6 FFT tests
│   └── test_filter.c       6 FIR filter tests
│
├── chapters/           ← YOU ARE HERE — start with 00-overview.md
├── reference/          ← Architecture docs + diagrams
│   ├── ARCHITECTURE.md
│   ├── API.md
│   └── diagrams/       PlantUML sources + PNG renders
│
├── Makefile            ← GNU Make (primary build)
└── CMakeLists.txt      ← CMake (cross-platform)
```

## Building & Running

```bash
# Build everything (debug mode, with warnings-as-errors)
make

# Build with optimisations
make release

# Run demos
./build/bin/fft_demo
./build/bin/filter_demo

# Run all 12 tests
make test

# Clean build artefacts
make clean
```

## Architecture at a Glance

> **📊 System Architecture** — [View full-size diagram →](../reference/diagrams/architecture.png)

The toolkit is organised in layers:

- **Application Layer** — Your code (demos, tests, and custom programs)
- **Core Library** — `fft.c`, `filter.c`, `dsp_utils.c`
- **System Interface** — File I/O, future ALSA audio
- **Platform Abstraction** — POSIX, math library

> **📊 Module Dependencies** — [View full-size diagram →](../reference/diagrams/modules.png)

Dependency rule: everything depends on `dsp_utils` (the `Complex` type lives
there). `fft.c` and `filter.c` are independent of each other.

---

**Next:** [Chapter 01 — Complex Numbers →](01-complex-numbers.md)

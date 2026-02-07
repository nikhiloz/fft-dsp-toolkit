# Chapter 18: Fixed-Point Arithmetic

Q-format, saturation, and fixed-point DSP implementation.

## Concept Diagram

![Fixed-Point Arithmetic](fixed_point_format.png)

## Contents

| File | Description |
|------|------------|
| [tutorial.md](tutorial.md) | Full theory tutorial with equations and exercises |
| [demo.c](demo.c) | Self-contained runnable demo |
| [`fixed_point.h`](../../include/fixed_point.h) | Library API |

## What You'll Learn

- Convert between floating-point and Q15/Q31 formats
- Perform addition with saturation to prevent overflow
- Implement Q15 × Q15 multiplication with correct shift
- Measure SQNR to evaluate fixed-point filter quality

---

[← Ch 17](../17-multirate-dsp/README.md) | [Index](../../reference/CHAPTER_INDEX.md) | [Ch 19 →](../19-advanced-fft/README.md)

# Chapter 16: Overlap-Add & Overlap-Save

Efficient block convolution for long FIR filtering.

## Concept Diagram

![Overlap-Add & Overlap-Save](ola_ols_pipeline.png)

## Contents

| File | Description |
|------|------------|
| [tutorial.md](tutorial.md) | Full theory tutorial with equations and exercises |
| [demo.c](demo.c) | Self-contained runnable demo |
| [`streaming.h`](../../include/streaming.h) | Library API |

## What You'll Learn

- Implement overlap-add block convolution
- Implement overlap-save block convolution
- Choose block size for optimal FFT performance
- Apply block filtering to streaming audio data

---

[← Ch 15](../15-correlation/README.md) | [Index](../../reference/CHAPTER_INDEX.md) | [Ch 17 →](../17-multirate-dsp/README.md)

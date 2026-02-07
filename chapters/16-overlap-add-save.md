# Chapter 16 — Overlap-Add/Save Streaming Convolution

| Navigation | |
|---|---|
| **Previous:** [Chapter 15 — Correlation](15-correlation.md) | **Next:** [Chapter 18 — Fixed-Point](18-fixed-point.md) |
| **Code:** [16-overlap-add-save.c](16-overlap-add-save.c) | **API:** [streaming.h](../include/streaming.h) |

---

## The Streaming Problem

When processing a long audio signal with an FIR filter:

- **Direct convolution**: O(L × M) per block — too slow for long filters
- **Single FFT**: Need the entire signal in memory — impossible for real-time

**Solution**: Process in blocks using **FFT-based convolution**, but handle
block boundaries correctly to avoid artefacts.

```
  Direct:  y[n] = Σ h[k]·x[n-k]    O(L·M) per L samples

  OLA:     Block → Pad → FFT → ×H → IFFT → overlap-add
           O(N·log N) per block, where N ≈ L + M
```

---

## Overlap-Add (OLA)

### Algorithm

1. **Partition** input into blocks of L samples
2. **Zero-pad** each block to N = L + M - 1 (round up to power-of-2)
3. **FFT** the zero-padded block → X[k]
4. **Multiply** X[k] · H[k] (pre-computed once)
5. **IFFT** → N-sample result
6. **Overlap-add**: first L samples go to output; last M-1 samples (the "tail")
   are added to the start of the next block's output

```
  Block 0: [===x0===|0000000]  → FFT×H→IFFT → [====y0=====|tail0]
  Block 1: [===x1===|0000000]  → FFT×H→IFFT → [====y1=====|tail1]

  Output:  [====y0=====]
                 [tail0 + y1====]
                          [tail1 + y2====]

  ✓ No boundary artefacts — output = exact linear convolution
```

### Why Zero-Padding?

Linear convolution of L-sample input with M-tap filter produces L+M-1 output
samples.  Without zero-padding, the FFT computes **circular** convolution,
which wraps around and corrupts the output.  Padding to N ≥ L+M-1 makes the
circular convolution equivalent to linear.

---

## Overlap-Save (OLS)

### Algorithm

1. Maintain a buffer of N samples: last M-1 from previous block + L new samples
2. **FFT** the entire N-sample buffer → X[k]
3. **Multiply** X[k] · H[k]
4. **IFFT** → N samples
5. **Discard** first M-1 samples (corrupted by circular wrap-around)
6. **Output** the remaining L = N - M + 1 valid samples

```
  Segment:  [prev M-1 | ===new L samples===]
             └overlap─┘

  → FFT → ×H → IFFT → [discard | ===valid L samples===]
                        └─M-1──┘
```

### OLA vs OLS

| Feature | Overlap-Add | Overlap-Save |
|---------|------------|--------------|
| Block output | L samples | L = N-M+1 samples |
| Accumulation | Yes (tail addition) | No |
| Memory | Extra tail buffer | Extra overlap buffer |
| Padding | Zero-pad each block | Maintain overlap region |
| Simpler to implement | ✓ | |

Both produce **identical results** (exact linear convolution).

---

## Efficiency Analysis

| Signal Len | Filter Taps | Direct (L×M) | OLA (≈5N·log₂N per block) | Speedup |
|------------|-------------|---------------|---------------------------|---------|
| 16,384     | 31          | 508K          | 205K                      | ×2.5    |
| 16,384     | 101         | 1.65M         | 287K                      | ×5.7    |
| 65,536     | 255         | 16.7M         | 1.15M                     | ×14.5   |

**Rule of thumb**: OLA/OLS becomes worthwhile when M > ~32 taps.

---

## Demo Walkthrough

### Demo 1: OLA Basic
Filters a 512-sample signal with a 31-tap lowpass using OLA.
Compares to direct FIR — max error should be ~10⁻¹⁴ (machine epsilon).

![OLA vs Direct](../plots/ch16/ola_vs_direct.png)
![OLA Reconstruction Error](../plots/ch16/ola_error.png)

### Demo 2: OLS Basic
Same test using Overlap-Save.  After warm-up block, matches direct FIR.

### Demo 3: OLA vs OLS
Processes a chirp through both methods with a 63-tap filter.
Both outputs match to machine precision.

### Demo 4: Streaming 16K Samples
Simulates real-time processing: 16384 samples of multi-tone audio
streamed through OLA in 128-sample blocks with a 101-tap lowpass.

![Streaming OLA](../plots/ch16/streaming_ola.png)

### Demo 5: Efficiency Table
Prints operation counts comparing direct FIR vs OLA across various
signal lengths and filter sizes.

---

## Key Takeaways

1. **OLA/OLS** enable real-time FFT-based FIR filtering in fixed-size blocks
2. **Zero-padding** to N ≥ L+M-1 converts circular convolution to linear
3. **OLA** adds tails between blocks; **OLS** discards corrupted prefix
4. Both produce **identical results** to direct convolution
5. Speedup grows with filter length — essential for long FIR filters (>32 taps)
6. **Pre-compute H[k]** once — amortised across all blocks

---

| Navigation | |
|---|---|
| **Previous:** [Chapter 15 — Correlation](15-correlation.md) | **Next:** [Chapter 18 — Fixed-Point](18-fixed-point.md) |

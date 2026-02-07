# Chapter Index — DSP Tutorial Suite

Quick-reference linking each chapter to its tutorial, demo, library modules, plots,
and test coverage.

> Regenerate plot PNGs with `make plots`. Rebuild all demos with `make release`.

---

## Part I — Foundations

| Ch | Tutorial | Demo | Library | Plots | Tests |
|----|----------|------|---------|-------|-------|
| [01](../chapters/01-signals-and-sequences.md) | Discrete-time signals | [`ch01`](../chapters/01-signals-and-sequences.c) | [`signal_gen.h`](../include/signal_gen.h) | [plots/ch01/](../plots/ch01/) (5 PNGs) | — |
| [02](../chapters/02-sampling-and-aliasing.md) | Sampling & aliasing | [`ch02`](../chapters/02-sampling-and-aliasing.c) | — | [plots/ch02/](../plots/ch02/) (3 PNGs) | — |
| [03](../chapters/03-complex-numbers.md) | Complex numbers | [`ch03`](../chapters/03-complex-numbers.c) | [`dsp_utils.h`](../include/dsp_utils.h) | [plots/ch03/](../plots/ch03/) (1 PNG) | — |
| [04](../chapters/04-lti-systems.md) | LTI systems | [`ch04`](../chapters/04-lti-systems.c) | [`convolution.h`](../include/convolution.h) | [plots/ch04/](../plots/ch04/) (2 PNGs) | — |

## Part II — Transform Domain

| Ch | Tutorial | Demo | Library | Plots | Tests |
|----|----------|------|---------|-------|-------|
| [05](../chapters/05-z-transform.md) | Z-Transform | [`ch05`](../chapters/05-z-transform.c) | — | [plots/ch05/](../plots/ch05/) (2 PNGs) | — |
| [06](../chapters/06-frequency-response.md) | Frequency response | [`ch06`](../chapters/06-frequency-response.c) | [`iir.h`](../include/iir.h) | [plots/ch06/](../plots/ch06/) (4 PNGs) | — |
| [07](../chapters/07-dft-theory.md) | DFT theory | [`ch07`](../chapters/07-dft-theory.c) | — | [plots/ch07/](../plots/ch07/) (3 PNGs) | — |
| [08](../chapters/08-fft-fundamentals.md) | FFT algorithms | [`ch08`](../chapters/08-fft-fundamentals.c) | [`fft.h`](../include/fft.h) | [plots/ch08/](../plots/ch08/) (1 PNG) | `test_fft` (6) |
| [09](../chapters/09-window-functions.md) | Window functions | [`ch09`](../chapters/09-window-functions.c) | [`dsp_utils.h`](../include/dsp_utils.h) | [plots/ch09/](../plots/ch09/) (2 PNGs) | — |

## Part III — Filter Design

| Ch | Tutorial | Demo | Library | Plots | Tests |
|----|----------|------|---------|-------|-------|
| [10](../chapters/10-digital-filters.md) | FIR filter design | [`ch10`](../chapters/10-digital-filters.c) | [`filter.h`](../include/filter.h) | [plots/ch10/](../plots/ch10/) (2 PNGs) | `test_filter` (6) |
| [11](../chapters/11-iir-filter-design.md) | IIR filter design | [`ch11`](../chapters/11-iir-filter-design.c) | [`iir.h`](../include/iir.h) | [plots/ch11/](../plots/ch11/) (3 PNGs) | `test_iir` (10) |
| [12](../chapters/12-filter-structures.md) | Filter structures | [`ch12`](../chapters/12-filter-structures.c) | — | [plots/ch12/](../plots/ch12/) (2 PNGs) | — |

## Part IV — Analysis

| Ch | Tutorial | Demo | Library | Plots | Tests |
|----|----------|------|---------|-------|-------|
| [13](../chapters/13-spectral-analysis.md) | Spectral analysis | [`ch13`](../chapters/13-spectral-analysis.c) | [`spectrum.h`](../include/spectrum.h) | [plots/ch13/](../plots/ch13/) (1 PNG) | — |
| [14](../chapters/14-psd-welch.md) | PSD & Welch's method | [`ch14`](../chapters/14-psd-welch.c) | [`spectrum.h`](../include/spectrum.h) | [plots/ch14/](../plots/ch14/) (7 PNGs) | `test_spectrum_corr` (12) |
| [15](../chapters/15-correlation.md) | Correlation | [`ch15`](../chapters/15-correlation.c) | [`correlation.h`](../include/correlation.h) | [plots/ch15/](../plots/ch15/) (5 PNGs) | `test_spectrum_corr` (12) |

## Part V — Advanced UG

| Ch | Tutorial | Demo | Library | Plots | Tests |
|----|----------|------|---------|-------|-------|
| [16](../chapters/16-overlap-add-save.md) | Overlap-Add/Save | [`ch16`](../chapters/16-overlap-add-save.c) | [`streaming.h`](../include/streaming.h) | [plots/ch16/](../plots/ch16/) (3 PNGs) | `test_phase4` (12) |
| [17](../chapters/17-multirate-dsp.md) | Multirate DSP | [`ch17`](../chapters/17-multirate-dsp.c) | [`multirate.h`](../include/multirate.h) | [plots/ch17/](../plots/ch17/) (3 PNGs) | `test_phase5` (15) |
| [18](../chapters/18-fixed-point.md) | Fixed-point arithmetic | [`ch18`](../chapters/18-fixed-point.c) | [`fixed_point.h`](../include/fixed_point.h) | [plots/ch18/](../plots/ch18/) (5 PNGs) | `test_phase4` (12) |
| [19](../chapters/19-advanced-fft.md) | Advanced FFT | [`ch19`](../chapters/19-advanced-fft.c) | [`advanced_fft.h`](../include/advanced_fft.h) | [plots/ch19/](../plots/ch19/) (3 PNGs) | `test_phase4` (12) |
| [20](../chapters/20-hilbert-transform.md) | Hilbert transform | [`ch20`](../chapters/20-hilbert-transform.c) | [`hilbert.h`](../include/hilbert.h) | [plots/ch20/](../plots/ch20/) (2 PNGs) | `test_phase5` (15) |
| [21](../chapters/21-signal-averaging.md) | Signal averaging | [`ch21`](../chapters/21-signal-averaging.c) | [`averaging.h`](../include/averaging.h) | [plots/ch21/](../plots/ch21/) (3 PNGs) | `test_phase5` (15) |
| [22](../chapters/22-advanced-fir.md) | Advanced FIR (Remez) | [`ch22`](../chapters/22-advanced-fir.c) | [`remez.h`](../include/remez.h) | [plots/ch22/](../plots/ch22/) (2 PNGs) | `test_phase5` (15) |

## Part VI — Postgraduate

| Ch | Tutorial | Demo | Library | Plots | Tests |
|----|----------|------|---------|-------|-------|
| [23](../chapters/23-adaptive-filters.md) | Adaptive filters | [`ch23`](../chapters/23-adaptive-filters.c) | [`adaptive.h`](../include/adaptive.h) | [plots/ch23/](../plots/ch23/) (2 PNGs) | `test_phase6` (19) |
| [24](../chapters/24-linear-prediction.md) | Linear prediction | [`ch24`](../chapters/24-linear-prediction.c) | [`lpc.h`](../include/lpc.h) | [plots/ch24/](../plots/ch24/) (2 PNGs) | `test_phase6` (19) |
| [25](../chapters/25-parametric-spectral.md) | Parametric spectral | [`ch25`](../chapters/25-parametric-spectral.c) | [`spectral_est.h`](../include/spectral_est.h) | [plots/ch25/](../plots/ch25/) (2 PNGs) | `test_phase6` (19) |
| [26](../chapters/26-cepstrum-mfcc.md) | Cepstrum & MFCC | [`ch26`](../chapters/26-cepstrum-mfcc.c) | [`cepstrum.h`](../include/cepstrum.h) | [plots/ch26/](../plots/ch26/) (3 PNGs) | `test_phase6` (19) |
| [27](../chapters/27-2d-dsp.md) | 2-D DSP | [`ch27`](../chapters/27-2d-dsp.c) | [`dsp2d.h`](../include/dsp2d.h) | [plots/ch27/](../plots/ch27/) (2 PNGs) | `test_phase6` (19) |

## Part VII — Applied / Capstone

| Ch | Tutorial | Demo | Library | Plots | Tests |
|----|----------|------|---------|-------|-------|
| [30](../chapters/30-putting-it-together.md) | Capstone pipeline | [`ch30`](../chapters/30-putting-it-together.c) | All | [plots/ch30/](../plots/ch30/) (2 PNGs) | — |

| Ch | Topic | Status |
|----|-------|--------|
| 28 | Real-time system design | Phase 7 |
| 29 | SIMD & hardware optimisation | Phase 7 |

---

## System-Level Documentation

| Document | Description |
|----------|-------------|
| [ARCHITECTURE.md](ARCHITECTURE.md) | Layered system design, module dependencies, roadmap |
| [API.md](API.md) | Complete public function reference |
| [diagrams/](diagrams/) | 10 PlantUML system diagrams (architecture, signal flow, modules, etc.) |

## Test Summary (80 tests)

| Suite | File | Count |
|-------|------|-------|
| FFT | [`test_fft.c`](../tests/test_fft.c) | 6 |
| Filter | [`test_filter.c`](../tests/test_filter.c) | 6 |
| IIR | [`test_iir.c`](../tests/test_iir.c) | 10 |
| Spectrum & Correlation | [`test_spectrum_corr.c`](../tests/test_spectrum_corr.c) | 12 |
| Phase 4 | [`test_phase4.c`](../tests/test_phase4.c) | 12 |
| Phase 5 | [`test_phase5.c`](../tests/test_phase5.c) | 15 |
| Phase 6 | [`test_phase6.c`](../tests/test_phase6.c) | 19 |
| **Total** | | **80** |

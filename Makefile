# Makefile - Simplified build for C-only DSP Tutorial Suite
# Requirements: gcc/clang, make

CC ?= gcc
CFLAGS := -Wall -Wextra -Werror -std=c99 -Iinclude -fPIC
CFLAGS_DEBUG := $(CFLAGS) -g -O0 -DDEBUG
CFLAGS_RELEASE := $(CFLAGS) -O3 -DNDEBUG
LDFLAGS := -lm

# Build directories
BUILD_DIR := build
BIN_DIR := $(BUILD_DIR)/bin
LIB_DIR := $(BUILD_DIR)/lib
OBJ_DIR := $(BUILD_DIR)/obj

# Source files
SOURCES := src/fft.c src/filter.c src/dsp_utils.c src/signal_gen.c src/convolution.c src/iir.c src/gnuplot.c src/spectrum.c src/correlation.c src/fixed_point.c src/advanced_fft.c src/streaming.c src/multirate.c src/hilbert.c src/averaging.c src/remez.c src/adaptive.c src/lpc.c src/spectral_est.c src/cepstrum.c src/dsp2d.c
OBJECTS := $(patsubst src/%.c, $(OBJ_DIR)/%.o, $(SOURCES))

TESTS := tests/test_fft.c tests/test_filter.c tests/test_iir.c tests/test_spectrum_corr.c tests/test_phase4.c tests/test_phase5.c tests/test_phase6.c

# Chapter demos
CHAPTER_DEMOS := chapters/01-signals-and-sequences.c \
	chapters/02-sampling-and-aliasing.c \
	chapters/03-complex-numbers.c \
	chapters/04-lti-systems.c \
	chapters/05-z-transform.c \
	chapters/06-frequency-response.c \
	chapters/07-dft-theory.c \
	chapters/08-fft-fundamentals.c \
	chapters/09-window-functions.c \
	chapters/10-digital-filters.c \
	chapters/11-iir-filter-design.c \
	chapters/12-filter-structures.c \
	chapters/13-spectral-analysis.c \
	chapters/14-psd-welch.c \
	chapters/15-correlation.c \
	chapters/16-overlap-add-save.c \
	chapters/18-fixed-point.c \
	chapters/19-advanced-fft.c \
	chapters/20-hilbert-transform.c \
	chapters/21-signal-averaging.c \
	chapters/22-advanced-fir.c \
	chapters/17-multirate-dsp.c \
	chapters/23-adaptive-filters.c \
	chapters/24-linear-prediction.c \
	chapters/25-parametric-spectral.c \
	chapters/26-cepstrum-mfcc.c \
	chapters/27-2d-dsp.c \
	chapters/30-putting-it-together.c

# Targets
all: release

# Create directories
$(BUILD_DIR) $(BIN_DIR) $(LIB_DIR) $(OBJ_DIR):
	mkdir -p $@

# Object compilation
$(OBJ_DIR)/%.o: src/%.c | $(OBJ_DIR)
	$(CC) $(CFLAGS_RELEASE) -c $< -o $@

# Debug build
debug: CFLAGS_RELEASE = $(CFLAGS_DEBUG)
debug: $(OBJ_DIR) $(BIN_DIR) $(LIB_DIR) \
	$(BIN_DIR)/ch01 \
	$(BIN_DIR)/ch02 \
	$(BIN_DIR)/ch03 \
	$(BIN_DIR)/ch04 \
	$(BIN_DIR)/ch05 \
	$(BIN_DIR)/ch06 \
	$(BIN_DIR)/ch07 \
	$(BIN_DIR)/ch08 \
	$(BIN_DIR)/ch09 \
	$(BIN_DIR)/ch10 \
	$(BIN_DIR)/ch11 \
	$(BIN_DIR)/ch12 \
	$(BIN_DIR)/ch13 \
	$(BIN_DIR)/ch14 \
	$(BIN_DIR)/ch15 \
	$(BIN_DIR)/ch16 \
	$(BIN_DIR)/ch18 \
	$(BIN_DIR)/ch19 \
	$(BIN_DIR)/ch20 \
	$(BIN_DIR)/ch21 \
	$(BIN_DIR)/ch22 \
	$(BIN_DIR)/ch17 \
	$(BIN_DIR)/ch23 \
	$(BIN_DIR)/ch24 \
	$(BIN_DIR)/ch25 \
	$(BIN_DIR)/ch26 \
	$(BIN_DIR)/ch27 \
	$(BIN_DIR)/ch30 \
	$(BIN_DIR)/test_fft \
	$(BIN_DIR)/test_filter \
	$(BIN_DIR)/test_iir \
	$(BIN_DIR)/test_spectrum_corr \
	$(BIN_DIR)/test_phase4 \
	$(BIN_DIR)/test_phase5 \
	$(BIN_DIR)/test_phase6 \
	$(BIN_DIR)/generate_plots

# Release build
release: $(OBJ_DIR) $(BIN_DIR) $(LIB_DIR) \
	$(BIN_DIR)/ch01 \
	$(BIN_DIR)/ch02 \
	$(BIN_DIR)/ch03 \
	$(BIN_DIR)/ch04 \
	$(BIN_DIR)/ch05 \
	$(BIN_DIR)/ch06 \
	$(BIN_DIR)/ch07 \
	$(BIN_DIR)/ch08 \
	$(BIN_DIR)/ch09 \
	$(BIN_DIR)/ch10 \
	$(BIN_DIR)/ch11 \
	$(BIN_DIR)/ch12 \
	$(BIN_DIR)/ch13 \
	$(BIN_DIR)/ch14 \
	$(BIN_DIR)/ch15 \
	$(BIN_DIR)/ch16 \
	$(BIN_DIR)/ch18 \
	$(BIN_DIR)/ch19 \
	$(BIN_DIR)/ch20 \
	$(BIN_DIR)/ch21 \
	$(BIN_DIR)/ch22 \
	$(BIN_DIR)/ch17 \
	$(BIN_DIR)/ch23 \
	$(BIN_DIR)/ch24 \
	$(BIN_DIR)/ch25 \
	$(BIN_DIR)/ch26 \
	$(BIN_DIR)/ch27 \
	$(BIN_DIR)/ch30 \
	$(BIN_DIR)/test_fft \
	$(BIN_DIR)/test_filter \
	$(BIN_DIR)/test_iir \
	$(BIN_DIR)/test_spectrum_corr \
	$(BIN_DIR)/test_phase4 \
	$(BIN_DIR)/test_phase5 \
	$(BIN_DIR)/test_phase6 \
	$(BIN_DIR)/generate_plots

# Static library
$(LIB_DIR)/libdsp_core.a: $(OBJECTS)
	ar rcs $@ $^

# Shared library
$(LIB_DIR)/libdsp_core.so: $(OBJECTS)
	$(CC) -shared -fPIC $(OBJECTS) $(LDFLAGS) -o $@

# Chapter demos
$(BIN_DIR)/ch01: chapters/01-signals-and-sequences.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) $< $(OBJECTS) $(LDFLAGS) -o $@

$(BIN_DIR)/ch02: chapters/02-sampling-and-aliasing.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) $< $(OBJECTS) $(LDFLAGS) -o $@

$(BIN_DIR)/ch03: chapters/03-complex-numbers.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) $< $(OBJECTS) $(LDFLAGS) -o $@

$(BIN_DIR)/ch04: chapters/04-lti-systems.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) $< $(OBJECTS) $(LDFLAGS) -o $@

$(BIN_DIR)/ch05: chapters/05-z-transform.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) $< $(OBJECTS) $(LDFLAGS) -o $@

$(BIN_DIR)/ch06: chapters/06-frequency-response.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) $< $(OBJECTS) $(LDFLAGS) -o $@

$(BIN_DIR)/ch07: chapters/07-dft-theory.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) $< $(OBJECTS) $(LDFLAGS) -o $@

$(BIN_DIR)/ch08: chapters/08-fft-fundamentals.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) $< $(OBJECTS) $(LDFLAGS) -o $@

$(BIN_DIR)/ch09: chapters/09-window-functions.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) $< $(OBJECTS) $(LDFLAGS) -o $@

$(BIN_DIR)/ch10: chapters/10-digital-filters.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) $< $(OBJECTS) $(LDFLAGS) -o $@

$(BIN_DIR)/ch11: chapters/11-iir-filter-design.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) $< $(OBJECTS) $(LDFLAGS) -o $@

$(BIN_DIR)/ch12: chapters/12-filter-structures.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) $< $(OBJECTS) $(LDFLAGS) -o $@

$(BIN_DIR)/ch13: chapters/13-spectral-analysis.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) $< $(OBJECTS) $(LDFLAGS) -o $@

$(BIN_DIR)/ch14: chapters/14-psd-welch.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) $< $(OBJECTS) $(LDFLAGS) -o $@

$(BIN_DIR)/ch15: chapters/15-correlation.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) $< $(OBJECTS) $(LDFLAGS) -o $@

$(BIN_DIR)/ch16: chapters/16-overlap-add-save.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) $< $(OBJECTS) $(LDFLAGS) -o $@

$(BIN_DIR)/ch18: chapters/18-fixed-point.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) $< $(OBJECTS) $(LDFLAGS) -o $@

$(BIN_DIR)/ch19: chapters/19-advanced-fft.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) $< $(OBJECTS) $(LDFLAGS) -o $@

$(BIN_DIR)/ch17: chapters/17-multirate-dsp.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) $< $(OBJECTS) $(LDFLAGS) -o $@

$(BIN_DIR)/ch20: chapters/20-hilbert-transform.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) $< $(OBJECTS) $(LDFLAGS) -o $@

$(BIN_DIR)/ch21: chapters/21-signal-averaging.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) $< $(OBJECTS) $(LDFLAGS) -o $@

$(BIN_DIR)/ch22: chapters/22-advanced-fir.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) $< $(OBJECTS) $(LDFLAGS) -o $@

$(BIN_DIR)/ch23: chapters/23-adaptive-filters.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) $< $(OBJECTS) $(LDFLAGS) -o $@

$(BIN_DIR)/ch24: chapters/24-linear-prediction.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) $< $(OBJECTS) $(LDFLAGS) -o $@

$(BIN_DIR)/ch25: chapters/25-parametric-spectral.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) $< $(OBJECTS) $(LDFLAGS) -o $@

$(BIN_DIR)/ch26: chapters/26-cepstrum-mfcc.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) $< $(OBJECTS) $(LDFLAGS) -o $@

$(BIN_DIR)/ch27: chapters/27-2d-dsp.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) $< $(OBJECTS) $(LDFLAGS) -o $@

$(BIN_DIR)/ch30: chapters/30-putting-it-together.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) $< $(OBJECTS) $(LDFLAGS) -o $@

# Plot generator
$(BIN_DIR)/generate_plots: tools/generate_plots.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) $< $(OBJECTS) $(LDFLAGS) -o $@

# Build only chapter demos
chapters: $(BIN_DIR)/ch01 $(BIN_DIR)/ch02 $(BIN_DIR)/ch03 $(BIN_DIR)/ch04 $(BIN_DIR)/ch05 $(BIN_DIR)/ch06 $(BIN_DIR)/ch07 $(BIN_DIR)/ch08 $(BIN_DIR)/ch09 $(BIN_DIR)/ch10 $(BIN_DIR)/ch11 $(BIN_DIR)/ch12 $(BIN_DIR)/ch13 $(BIN_DIR)/ch14 $(BIN_DIR)/ch15 $(BIN_DIR)/ch16 $(BIN_DIR)/ch17 $(BIN_DIR)/ch18 $(BIN_DIR)/ch19 $(BIN_DIR)/ch20 $(BIN_DIR)/ch21 $(BIN_DIR)/ch22 $(BIN_DIR)/ch23 $(BIN_DIR)/ch24 $(BIN_DIR)/ch25 $(BIN_DIR)/ch26 $(BIN_DIR)/ch27 $(BIN_DIR)/ch30

# Tests
$(BIN_DIR)/test_fft: tests/test_fft.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) -Itests $< $(OBJECTS) $(LDFLAGS) -o $@

$(BIN_DIR)/test_filter: tests/test_filter.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) -Itests $< $(OBJECTS) $(LDFLAGS) -o $@

$(BIN_DIR)/test_iir: tests/test_iir.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) -Itests $< $(OBJECTS) $(LDFLAGS) -o $@

$(BIN_DIR)/test_spectrum_corr: tests/test_spectrum_corr.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) -Itests $< $(OBJECTS) $(LDFLAGS) -o $@

$(BIN_DIR)/test_phase4: tests/test_phase4.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) -Itests $< $(OBJECTS) $(LDFLAGS) -o $@

$(BIN_DIR)/test_phase5: tests/test_phase5.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) -Itests $< $(OBJECTS) $(LDFLAGS) -o $@

$(BIN_DIR)/test_phase6: tests/test_phase6.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) -Itests $< $(OBJECTS) $(LDFLAGS) -o $@

# Run tests
test: $(BIN_DIR)/test_fft $(BIN_DIR)/test_filter $(BIN_DIR)/test_iir $(BIN_DIR)/test_spectrum_corr $(BIN_DIR)/test_phase4 $(BIN_DIR)/test_phase5 $(BIN_DIR)/test_phase6
	@echo "=== Running FFT tests ==="
	$(BIN_DIR)/test_fft
	@echo "\n=== Running Filter tests ==="
	$(BIN_DIR)/test_filter
	@echo "\n=== Running IIR tests ==="
	$(BIN_DIR)/test_iir
	@echo "\n=== Running Spectrum & Correlation tests ==="
	$(BIN_DIR)/test_spectrum_corr
	@echo "\n=== Running Phase 4 tests ==="
	$(BIN_DIR)/test_phase4
	@echo "\n=== Running Phase 5 tests ==="
	$(BIN_DIR)/test_phase5
	@echo "\n=== Running Phase 6 tests ==="
	$(BIN_DIR)/test_phase6

# Run chapter demos
run: chapters
	@echo "=== Ch01: Signals & Sequences ==="
	$(BIN_DIR)/ch01
	@echo "\n=== Ch02: Sampling & Aliasing ==="
	$(BIN_DIR)/ch02
	@echo "\n=== Ch03: Complex Numbers ==="
	$(BIN_DIR)/ch03
	@echo "\n=== Ch04: LTI Systems & Convolution ==="
	$(BIN_DIR)/ch04
	@echo "\n=== Ch05: Z-Transform ==="
	$(BIN_DIR)/ch05
	@echo "\n=== Ch07: DFT Theory ==="
	$(BIN_DIR)/ch07
	@echo "\n=== Ch08: FFT Fundamentals ==="
	$(BIN_DIR)/ch08
	@echo "\n=== Ch09: Window Functions ==="
	$(BIN_DIR)/ch09
	@echo "\n=== Ch10: Digital Filters ==="
	$(BIN_DIR)/ch10
	@echo "\n=== Ch11: IIR Filter Design ==="
	$(BIN_DIR)/ch11
	@echo "\n=== Ch12: Filter Structures ==="
	$(BIN_DIR)/ch12
	@echo "\n=== Ch13: Spectral Analysis ==="
	$(BIN_DIR)/ch13
	@echo "\n=== Ch14: PSD & Welch ==="
	$(BIN_DIR)/ch14
	@echo "\n=== Ch15: Correlation ==="
	$(BIN_DIR)/ch15
	@echo "\n=== Ch16: Overlap-Add/Save ==="
	$(BIN_DIR)/ch16
	@echo "\n=== Ch18: Fixed-Point ==="
	$(BIN_DIR)/ch18
	@echo "\n=== Ch17: Multirate DSP ==="
	$(BIN_DIR)/ch17
	@echo "\n=== Ch19: Advanced FFT ==="
	$(BIN_DIR)/ch19
	@echo "\n=== Ch20: Hilbert Transform ==="
	$(BIN_DIR)/ch20
	@echo "\n=== Ch21: Signal Averaging ==="
	$(BIN_DIR)/ch21
	@echo "\n=== Ch22: Advanced FIR ==="
	$(BIN_DIR)/ch22
	@echo "\n=== Ch23: Adaptive Filters ==="
	$(BIN_DIR)/ch23
	@echo "\n=== Ch24: Linear Prediction ==="
	$(BIN_DIR)/ch24
	@echo "\n=== Ch25: Parametric Spectral ==="
	$(BIN_DIR)/ch25
	@echo "\n=== Ch26: Cepstrum & MFCC ==="
	$(BIN_DIR)/ch26
	@echo "\n=== Ch27: 2D DSP ==="
	$(BIN_DIR)/ch27
	@echo "\n=== Ch30: Putting It Together ==="
	$(BIN_DIR)/ch30

# Generate all gnuplot PNG visualizations (requires gnuplot)
plots: $(BIN_DIR)/generate_plots
	@echo "=== Generating all plots ==="
	$(BIN_DIR)/generate_plots

# Code formatting
format:
	clang-format -i src/*.c include/*.h chapters/*.c tests/test_*.c

# Static analysis
lint:
	clang-tidy src/*.c -- -Iinclude

# Memory checking
memcheck: debug
	valgrind --leak-check=full --error-exitcode=1 $(BIN_DIR)/test_fft
	valgrind --leak-check=full --error-exitcode=1 $(BIN_DIR)/test_filter

# Profiling (Linux only)
profile: $(BIN_DIR)/ch08
	perf record -g $(BIN_DIR)/ch08
	perf report

# Clean
clean:
	rm -rf $(BUILD_DIR)

# Deep clean
distclean: clean
	find . -name "*.o" -o -name "*.a" -o -name "*.so" -o -name "perf.data*" | xargs rm -f

# Install
install: release
	@echo "Installing to /usr/local..."
	mkdir -p /usr/local/include/dsp_core /usr/local/lib
	cp include/*.h /usr/local/include/dsp_core/
	cp $(LIB_DIR)/* /usr/local/lib/
	ldconfig

# Help
help:
	@echo "DSP Tutorial Suite Makefile"
	@echo ""
	@echo "Targets:"
	@echo "  make release     - Build release version (default)"
	@echo "  make debug       - Build debug version with symbols"
	@echo "  make test        - Run unit tests"
	@echo "  make run         - Run all chapter demos"
	@echo "  make chapters    - Build chapter demos only"
	@echo "  make plots       - Generate all gnuplot visualizations"
	@echo "  make memcheck    - Run tests with valgrind"
	@echo "  make profile     - Profile ch08 (FFT) with perf"
	@echo "  make format      - Format code with clang-format"
	@echo "  make lint        - Static analysis with clang-tidy"
	@echo "  make install     - Install headers & libraries to /usr/local"
	@echo "  make clean       - Remove build directory"
	@echo "  make distclean   - Remove all generated files"
	@echo "  make help        - Show this help message"

.PHONY: all debug release test run chapters plots memcheck profile format lint clean distclean install help

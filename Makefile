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
SOURCES := src/fft.c src/filter.c src/dsp_utils.c src/signal_gen.c src/convolution.c src/iir.c
OBJECTS := $(patsubst src/%.c, $(OBJ_DIR)/%.o, $(SOURCES))

TESTS := tests/test_fft.c tests/test_filter.c tests/test_iir.c

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
	$(BIN_DIR)/ch30 \
	$(BIN_DIR)/test_fft \
	$(BIN_DIR)/test_filter \
	$(BIN_DIR)/test_iir

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
	$(BIN_DIR)/ch30 \
	$(BIN_DIR)/test_fft \
	$(BIN_DIR)/test_filter \
	$(BIN_DIR)/test_iir

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

$(BIN_DIR)/ch30: chapters/30-putting-it-together.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) $< $(OBJECTS) $(LDFLAGS) -o $@

# Build only chapter demos
chapters: $(BIN_DIR)/ch01 $(BIN_DIR)/ch02 $(BIN_DIR)/ch03 $(BIN_DIR)/ch04 $(BIN_DIR)/ch05 $(BIN_DIR)/ch06 $(BIN_DIR)/ch07 $(BIN_DIR)/ch08 $(BIN_DIR)/ch09 $(BIN_DIR)/ch10 $(BIN_DIR)/ch11 $(BIN_DIR)/ch12 $(BIN_DIR)/ch13 $(BIN_DIR)/ch30

# Tests
$(BIN_DIR)/test_fft: tests/test_fft.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) -Itests $< $(OBJECTS) $(LDFLAGS) -o $@

$(BIN_DIR)/test_filter: tests/test_filter.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) -Itests $< $(OBJECTS) $(LDFLAGS) -o $@

$(BIN_DIR)/test_iir: tests/test_iir.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) -Itests $< $(OBJECTS) $(LDFLAGS) -o $@

# Run tests
test: $(BIN_DIR)/test_fft $(BIN_DIR)/test_filter $(BIN_DIR)/test_iir
	@echo "=== Running FFT tests ==="
	$(BIN_DIR)/test_fft
	@echo "\n=== Running Filter tests ==="
	$(BIN_DIR)/test_filter
	@echo "\n=== Running IIR tests ==="
	$(BIN_DIR)/test_iir

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
	@echo "\n=== Ch30: Putting It Together ==="
	$(BIN_DIR)/ch30

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
	@echo "  make memcheck    - Run tests with valgrind"
	@echo "  make profile     - Profile ch08 (FFT) with perf"
	@echo "  make format      - Format code with clang-format"
	@echo "  make lint        - Static analysis with clang-tidy"
	@echo "  make install     - Install headers & libraries to /usr/local"
	@echo "  make clean       - Remove build directory"
	@echo "  make distclean   - Remove all generated files"
	@echo "  make help        - Show this help message"

.PHONY: all debug release test run chapters memcheck profile format lint clean distclean install help

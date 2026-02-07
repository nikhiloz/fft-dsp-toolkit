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
SOURCES := src/fft.c src/filter.c src/dsp_utils.c src/signal_gen.c
OBJECTS := $(patsubst src/%.c, $(OBJ_DIR)/%.o, $(SOURCES))

TESTS := tests/test_fft.c tests/test_filter.c

# Chapter demos
CHAPTER_DEMOS := chapters/01-signals-and-sequences.c \
	chapters/01-complex-numbers.c \
	chapters/02-fft-fundamentals.c \
	chapters/03-window-functions.c \
	chapters/04-digital-filters.c \
	chapters/05-spectral-analysis.c \
	chapters/08-putting-it-together.c

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
	$(BIN_DIR)/ch01s \
	$(BIN_DIR)/ch01 \
	$(BIN_DIR)/ch02 \
	$(BIN_DIR)/ch03 \
	$(BIN_DIR)/ch04 \
	$(BIN_DIR)/ch05 \
	$(BIN_DIR)/ch08 \
	$(BIN_DIR)/test_fft \
	$(BIN_DIR)/test_filter

# Release build
release: $(OBJ_DIR) $(BIN_DIR) $(LIB_DIR) \
	$(BIN_DIR)/ch01s \
	$(BIN_DIR)/ch01 \
	$(BIN_DIR)/ch02 \
	$(BIN_DIR)/ch03 \
	$(BIN_DIR)/ch04 \
	$(BIN_DIR)/ch05 \
	$(BIN_DIR)/ch08 \
	$(BIN_DIR)/test_fft \
	$(BIN_DIR)/test_filter

# Static library
$(LIB_DIR)/libfft_dsp.a: $(OBJECTS)
	ar rcs $@ $^

# Shared library
$(LIB_DIR)/libfft_dsp.so: $(OBJECTS)
	$(CC) -shared -fPIC $(OBJECTS) $(LDFLAGS) -o $@

# Chapter demos
$(BIN_DIR)/ch01s: chapters/01-signals-and-sequences.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) $< $(OBJECTS) $(LDFLAGS) -o $@

$(BIN_DIR)/ch01: chapters/01-complex-numbers.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) $< $(OBJECTS) $(LDFLAGS) -o $@

$(BIN_DIR)/ch02: chapters/02-fft-fundamentals.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) $< $(OBJECTS) $(LDFLAGS) -o $@

$(BIN_DIR)/ch03: chapters/03-window-functions.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) $< $(OBJECTS) $(LDFLAGS) -o $@

$(BIN_DIR)/ch04: chapters/04-digital-filters.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) $< $(OBJECTS) $(LDFLAGS) -o $@

$(BIN_DIR)/ch05: chapters/05-spectral-analysis.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) $< $(OBJECTS) $(LDFLAGS) -o $@

$(BIN_DIR)/ch08: chapters/08-putting-it-together.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) $< $(OBJECTS) $(LDFLAGS) -o $@

# Build only chapter demos
chapters: $(BIN_DIR)/ch01s $(BIN_DIR)/ch01 $(BIN_DIR)/ch02 $(BIN_DIR)/ch03 $(BIN_DIR)/ch04 $(BIN_DIR)/ch05 $(BIN_DIR)/ch08

# Tests
$(BIN_DIR)/test_fft: tests/test_fft.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) -Itests $< $(OBJECTS) $(LDFLAGS) -o $@

$(BIN_DIR)/test_filter: tests/test_filter.c $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS_RELEASE) -Itests $< $(OBJECTS) $(LDFLAGS) -o $@

# Run tests
test: $(BIN_DIR)/test_fft $(BIN_DIR)/test_filter
	@echo "=== Running FFT tests ==="
	$(BIN_DIR)/test_fft
	@echo "\n=== Running Filter tests ==="
	$(BIN_DIR)/test_filter

# Run chapter demos
run: chapters
	@echo "=== Ch01s: Signals & Sequences ==="
	$(BIN_DIR)/ch01s
	@echo "\n=== Ch01: Complex Numbers ==="
	$(BIN_DIR)/ch01
	@echo "\n=== Ch02: FFT Fundamentals ==="
	$(BIN_DIR)/ch02
	@echo "\n=== Ch03: Window Functions ==="
	$(BIN_DIR)/ch03
	@echo "\n=== Ch04: Digital Filters ==="
	$(BIN_DIR)/ch04
	@echo "\n=== Ch05: Spectral Analysis ==="
	$(BIN_DIR)/ch05
	@echo "\n=== Ch08: Putting It Together ==="
	$(BIN_DIR)/ch08

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
profile: $(BIN_DIR)/ch02
	perf record -g $(BIN_DIR)/ch02
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
	mkdir -p /usr/local/include/fft_dsp /usr/local/lib
	cp include/*.h /usr/local/include/fft_dsp/
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
	@echo "  make profile     - Profile ch02 (FFT) with perf"
	@echo "  make format      - Format code with clang-format"
	@echo "  make lint        - Static analysis with clang-tidy"
	@echo "  make install     - Install headers & libraries to /usr/local"
	@echo "  make clean       - Remove build directory"
	@echo "  make distclean   - Remove all generated files"
	@echo "  make help        - Show this help message"

.PHONY: all debug release test run chapters memcheck profile format lint clean distclean install help

/**
 * @file fft.c
 * @brief Cooley-Tukey Radix-2 Decimation-In-Time FFT.
 *
 * TUTORIAL CROSS-REFERENCES:
 *   Theory & butterfly diagram  → chapters/02-fft-fundamentals.md
 *   Complex multiplication used → chapters/01-complex-numbers.md
 *   Spectral leakage & windows  → chapters/03-window-functions.md
 *
 * ALGORITHM OVERVIEW:
 *   1. Bit-reversal permutation  (reorder indices)
 *   2. log₂(N) stages of butterfly operations
 *   3. Each butterfly: X[k] = E[k] + W·O[k],  X[k+N/2] = E[k] - W·O[k]
 *      where W = exp(-j·2π·k/N) is the "twiddle factor"
 *
 * Complexity: O(N log N) — compared to O(N²) for the naive DFT.
 */

#define _GNU_SOURCE
#include "fft.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ════════════════════════════════════════════════════════════════════
 *  STEP 1: Bit-reversal permutation
 *
 *  The Cooley-Tukey algorithm processes data in "bit-reversed" order.
 *  For N=8, the mapping is:
 *     index:    0  1  2  3  4  5  6  7
 *     binary:  000 001 010 011 100 101 110 111
 *     reverse: 000 100 010 110 001 101 011 111
 *     result:   0  4  2  6  1  5  3  7
 *
 *  This swap ensures that after all butterfly stages, the output
 *  is in natural (sequential) order.
 * ════════════════════════════════════════════════════════════════════ */

static void bit_reverse_permute(Complex *x, int n) {
    int j = 0;
    for (int i = 0; i < n - 1; i++) {
        if (i < j) {
            /* Swap x[i] and x[j] */
            Complex tmp = x[i];
            x[i] = x[j];
            x[j] = tmp;
        }
        int m = n >> 1;
        while (m >= 1 && j >= m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
}

/* ════════════════════════════════════════════════════════════════════
 *  STEP 2: Butterfly operations  (the core of the FFT)
 *
 *  For each stage s (1 to log₂N):
 *    - Group size = 2^s
 *    - Half group = 2^(s-1)
 *    - Twiddle factor W = exp(-j · 2π / group_size)
 *
 *  The "butterfly" combines two values:
 *
 *      a ──────┬──── a + W·b     ("top wing")
 *              ×  W
 *      b ──────┴──── a − W·b     ("bottom wing")
 *
 *  This is where all the computational work happens.
 *  With N samples and log₂(N) stages, there are N/2 butterflies per
 *  stage → total N/2 · log₂(N) complex multiplications.
 *
 *  See chapters/02-fft-fundamentals.md § "Butterfly Diagram" for
 *  a visual walkthrough of an 8-point FFT.
 * ════════════════════════════════════════════════════════════════════ */

void fft(Complex *x, int n) {
    if (n <= 1) return;

    /* Step 1: reorder data by bit-reversal */
    bit_reverse_permute(x, n);

    /* Step 2: butterfly stages */
    for (int stage_size = 2; stage_size <= n; stage_size <<= 1) {
        int half = stage_size >> 1;

        /*
         * Twiddle factor for this stage:
         *   W_base = exp(-j · 2π / stage_size)
         * We start with W = 1 and multiply by W_base at each step
         * within the group.
         */
        double angle = -2.0 * M_PI / stage_size;
        Complex w_base = { cos(angle), sin(angle) };

        /* Process each group of 'stage_size' elements */
        for (int group = 0; group < n; group += stage_size) {
            Complex w = { 1.0, 0.0 };  /* W^0 = 1 */

            for (int k = 0; k < half; k++) {
                int top = group + k;
                int bot = group + k + half;

                /*
                 * THE BUTTERFLY:
                 *   t       = W · x[bot]        (twiddle multiply)
                 *   x[top]  = x[top] + t        (top wing)
                 *   x[bot]  = x[top_old] - t    (bottom wing)
                 */
                Complex t = complex_mul(w, x[bot]);
                Complex u = x[top];

                x[top] = complex_add(u, t);
                x[bot] = complex_sub(u, t);

                /* Advance twiddle: W = W · W_base */
                w = complex_mul(w, w_base);
            }
        }
    }
}

/* ════════════════════════════════════════════════════════════════════
 *  Real-valued FFT wrapper
 *  Copies real input into complex array, then runs FFT.
 * ════════════════════════════════════════════════════════════════════ */

void fft_real(const double *in, Complex *out, int n) {
    for (int i = 0; i < n; i++) {
        out[i].re = in[i];
        out[i].im = 0.0;
    }
    fft(out, n);
}

/* ════════════════════════════════════════════════════════════════════
 *  Inverse FFT
 *
 *  IFFT is just FFT with:
 *    1. Conjugated input  (negate imaginary parts)
 *    2. Run FFT
 *    3. Conjugate output again
 *    4. Divide by N
 *
 *  This works because the DFT matrix is unitary up to a scale factor.
 * ════════════════════════════════════════════════════════════════════ */

void ifft(Complex *x, int n) {
    /* Step 1: conjugate */
    for (int i = 0; i < n; i++) {
        x[i].im = -x[i].im;
    }

    /* Step 2: forward FFT */
    fft(x, n);

    /* Step 3: conjugate again and scale by 1/N */
    double scale = 1.0 / n;
    for (int i = 0; i < n; i++) {
        x[i].re *= scale;
        x[i].im = -x[i].im * scale;
    }
}

/* ════════════════════════════════════════════════════════════════════
 *  Feature extraction helpers
 * ════════════════════════════════════════════════════════════════════ */

void fft_magnitude(const Complex *x, double *mag, int n) {
    for (int i = 0; i < n; i++) {
        mag[i] = complex_mag(x[i]);
    }
}

void fft_phase(const Complex *x, double *phase, int n) {
    for (int i = 0; i < n; i++) {
        phase[i] = complex_phase(x[i]);
    }
}

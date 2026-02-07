/**
 * @file 01-complex-numbers.c
 * @brief Chapter 1 demo — Complex number arithmetic in action.
 *
 * Demonstrates:
 *   - Creating complex numbers (rectangular & polar)
 *   - Addition, subtraction, multiplication
 *   - Magnitude and phase extraction
 *   - The twiddle factor W_N used in the FFT
 *
 * Build & run:
 *   make chapters && ./build/bin/ch01
 *
 * Read alongside: chapters/01-complex-numbers.md
 */

#include <stdio.h>
#include <math.h>
#include "dsp_utils.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static void print_complex(const char *label, Complex z) {
    printf("  %-20s = %+.4f %+.4fi\n", label, z.re, z.im);
}

int main(void) {
    printf("=== Chapter 1: Complex Numbers ===\n\n");

    /* ── Basic arithmetic ────────────────────────────────────────── */
    Complex a = { 3.0, 4.0 };   /* 3 + 4i */
    Complex b = { 1.0, -2.0 };  /* 1 - 2i */

    printf("Given:\n");
    print_complex("a", a);
    print_complex("b", b);
    printf("\nArithmetic:\n");
    print_complex("a + b", complex_add(a, b));
    print_complex("a - b", complex_sub(a, b));
    print_complex("a * b", complex_mul(a, b));

    printf("\n  |a| = %.4f   (magnitude)\n", complex_mag(a));
    printf("  ∠a = %.4f rad = %.1f°\n",
           complex_phase(a), complex_phase(a) * 180.0 / M_PI);

    /* ── Polar form ──────────────────────────────────────────────── */
    printf("\nPolar form:\n");
    Complex unit = complex_from_polar(1.0, M_PI / 4.0);
    printf("  e^(iπ/4) = ");
    print_complex("", unit);
    printf("  Expected:  0.7071 + 0.7071i (unit circle at 45°)\n");

    /* ── Twiddle factors (preview of FFT) ────────────────────────── */
    printf("\nTwiddle factors (W_8^k = e^(-j·2π·k/8)):\n");
    for (int k = 0; k < 8; k++) {
        double angle = -2.0 * M_PI * k / 8.0;
        Complex w = complex_from_polar(1.0, angle);
        printf("  W_8^%d = %+.4f %+.4fi   |W| = %.4f\n",
               k, w.re, w.im, complex_mag(w));
    }
    printf("\n  Note: all magnitudes are 1.0 — twiddle factors\n");
    printf("  are rotations on the unit circle.\n");

    /* ── Verify multiplication identity ──────────────────────────── */
    printf("\nVerification: (3+4i)(1-2i) = ?\n");
    Complex product = complex_mul(a, b);
    printf("  By hand:   (3·1 - 4·(-2)) + (3·(-2) + 4·1)i = 11 - 2i\n");
    printf("  By code:   %.0f %+.0fi  ✓\n", product.re, product.im);

    return 0;
}

/**
 * @file test_phase7.c
 * @brief Unit tests for Phase 7 modules: realtime, optimization.
 *
 * Tests:
 *   1.  Ring buffer create/destroy
 *   2.  Ring buffer write/read round-trip
 *   3.  Ring buffer peek does not consume
 *   4.  Ring buffer skip discards correctly
 *   5.  Ring buffer wraps around correctly
 *   6.  Ring buffer full/empty detection
 *   7.  Frame processor create/destroy
 *   8.  Frame processor detects sine frequency
 *   9.  Latency timer returns increasing values
 *  10.  Latency stats min/max/avg
 *  11.  Radix-4 FFT matches radix-2
 *  12.  Radix-4 FFT round-trip (fft → ifft)
 *  13.  Radix-4 fallback for non-power-of-4
 *  14.  Twiddle table FFT matches direct FFT
 *  15.  Aligned memory is properly aligned
 *  16.  Aligned memory read/write works
 *  17.  Benchmark produces valid results
 *  18.  Bench print does not crash
 *
 * Run: make test
 */

#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "test_framework.h"
#include "realtime.h"
#include "optimization.h"
#include "fft.h"
#include "dsp_utils.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void)
{
    TEST_SUITE("Phase 7: Real-Time & Optimisation");

    /* ─── Ring Buffer Tests ─────────────────────────────── */

    TEST_CASE_BEGIN("Ring buffer create/destroy");
    {
        RingBuffer *rb = ring_buffer_create(100);
        int ok = (rb != NULL) && (rb->cap >= 100);
        /* Power of 2: 128 >= 100 */
        ok = ok && (rb->cap == 128);
        ring_buffer_destroy(rb);
        if (ok) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("bad capacity rounding"); }
    }

    TEST_CASE_BEGIN("Ring buffer write/read round-trip");
    {
        RingBuffer *rb = ring_buffer_create(64);
        double data[] = {1.0, 2.0, 3.0, 4.0, 5.0};
        int written = ring_buffer_write(rb, data, 5);
        double out[5] = {0};
        int nread = ring_buffer_read(rb, out, 5);

        int ok = (written == 5) && (nread == 5);
        for (int i = 0; i < 5; i++)
            ok = ok && (fabs(out[i] - data[i]) < 1e-15);

        ring_buffer_destroy(rb);
        if (ok) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("data mismatch"); }
    }

    TEST_CASE_BEGIN("Ring buffer peek does not consume");
    {
        RingBuffer *rb = ring_buffer_create(16);
        double data[] = {10.0, 20.0, 30.0};
        ring_buffer_write(rb, data, 3);

        double peek_buf[3] = {0};
        ring_buffer_peek(rb, peek_buf, 3);
        int avail_after_peek = ring_buffer_available(rb);

        double read_buf[3] = {0};
        ring_buffer_read(rb, read_buf, 3);

        int ok = (avail_after_peek == 3);
        for (int i = 0; i < 3; i++)
            ok = ok && (fabs(peek_buf[i] - data[i]) < 1e-15)
                    && (fabs(read_buf[i] - data[i]) < 1e-15);

        ring_buffer_destroy(rb);
        if (ok) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("peek consumed data"); }
    }

    TEST_CASE_BEGIN("Ring buffer skip discards correctly");
    {
        RingBuffer *rb = ring_buffer_create(16);
        double data[] = {1.0, 2.0, 3.0, 4.0, 5.0};
        ring_buffer_write(rb, data, 5);
        ring_buffer_skip(rb, 2);

        int avail = ring_buffer_available(rb);
        double out[3] = {0};
        ring_buffer_read(rb, out, 3);

        int ok = (avail == 3) &&
                 fabs(out[0] - 3.0) < 1e-15 &&
                 fabs(out[1] - 4.0) < 1e-15 &&
                 fabs(out[2] - 5.0) < 1e-15;

        ring_buffer_destroy(rb);
        if (ok) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("skip failed"); }
    }

    TEST_CASE_BEGIN("Ring buffer wraps around correctly");
    {
        RingBuffer *rb = ring_buffer_create(8);  /* cap = 8 */
        /* Fill nearly full, then read, then write more to force wrap */
        double d1[] = {1, 2, 3, 4, 5, 6};
        ring_buffer_write(rb, d1, 6);  /* 6 out of 7 usable */
        double tmp[4];
        ring_buffer_read(rb, tmp, 4);  /* read 4, avail=2 */

        double d2[] = {7, 8, 9, 10};
        ring_buffer_write(rb, d2, 4);  /* wraps around */

        int avail = ring_buffer_available(rb);
        double out[6] = {0};
        int nread = ring_buffer_read(rb, out, 6);

        int ok = (avail == 6) && (nread == 6);
        double expected[] = {5, 6, 7, 8, 9, 10};
        for (int i = 0; i < 6; i++)
            ok = ok && (fabs(out[i] - expected[i]) < 1e-15);

        ring_buffer_destroy(rb);
        if (ok) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("wrap error"); }
    }

    TEST_CASE_BEGIN("Ring buffer full/empty detection");
    {
        RingBuffer *rb = ring_buffer_create(4);  /* cap = 4, usable = 3 */
        int empty_ok = (ring_buffer_available(rb) == 0);
        int space_ok = (ring_buffer_space(rb) == 3);

        double d[] = {1, 2, 3, 4};
        int w = ring_buffer_write(rb, d, 4);  /* only 3 fit */
        int full_avail = ring_buffer_available(rb);
        int full_space = ring_buffer_space(rb);

        ring_buffer_destroy(rb);
        int ok = empty_ok && space_ok && (w == 3) &&
                 (full_avail == 3) && (full_space == 0);
        if (ok) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("full/empty detection failed"); }
    }

    /* ─── Frame Processor Tests ─────────────────────────── */

    TEST_CASE_BEGIN("Frame processor create/destroy");
    {
        FrameProcessor *fp = frame_processor_create(256, 128);
        int ok = (fp != NULL) && (fp->frame_size == 256) && (fp->hop_size == 128);
        frame_processor_destroy(fp);
        if (ok) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("bad frame processor"); }
    }

    TEST_CASE_BEGIN("Frame processor detects sine frequency");
    {
        int N = 512;
        double fs = 8000.0;
        double freq = 1000.0;  /* easy to detect */
        FrameProcessor *fp = frame_processor_create(N, N / 2);

        /* Generate a pure sine */
        double *sig = (double *)malloc((size_t)N * sizeof(double));
        for (int i = 0; i < N; i++)
            sig[i] = sin(2.0 * M_PI * freq * i / fs);

        frame_processor_feed(fp, sig, N);
        double detected = frame_processor_peak_freq(fp, fs);

        /* Should be within one FFT bin of 1000 Hz */
        double bin_width = fs / N;
        int ok = fabs(detected - freq) <= bin_width;

        free(sig);
        frame_processor_destroy(fp);
        if (ok) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("wrong peak frequency"); }
    }

    /* ─── Latency Timer Tests ───────────────────────────── */

    TEST_CASE_BEGIN("Latency timer returns increasing values");
    {
        double t0 = timer_usec();
        /* Busy wait a tiny bit */
        volatile double sum = 0;
        for (int i = 0; i < 10000; i++) sum += 0.001;
        (void)sum;
        double t1 = timer_usec();

        int ok = (t1 > t0) && ((t1 - t0) > 0.1);  /* at least 0.1 µs */
        if (ok) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("timer not monotonic"); }
    }

    TEST_CASE_BEGIN("Latency stats min/max/avg");
    {
        LatencyStats ls;
        latency_init(&ls);
        latency_record(&ls, 10.0);
        latency_record(&ls, 20.0);
        latency_record(&ls, 30.0);

        int ok = (fabs(ls.min_us - 10.0) < 1e-10) &&
                 (fabs(ls.max_us - 30.0) < 1e-10) &&
                 (fabs(latency_avg(&ls) - 20.0) < 1e-10) &&
                 (ls.count == 3);

        if (ok) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("latency stats wrong"); }
    }

    /* ─── Radix-4 FFT Tests ─────────────────────────────── */

    TEST_CASE_BEGIN("Radix-4 FFT matches radix-2");
    {
        int N = 256;
        Complex *x2 = (Complex *)malloc((size_t)N * sizeof(Complex));
        Complex *x4 = (Complex *)malloc((size_t)N * sizeof(Complex));

        for (int i = 0; i < N; i++) {
            double t = (double)i / N;
            x2[i].re = sin(2.0 * M_PI * 5.0 * t) + cos(2.0 * M_PI * 17.0 * t);
            x2[i].im = 0.0;
            x4[i] = x2[i];
        }

        fft(x2, N);
        fft_radix4(x4, N);

        double max_err = 0;
        for (int i = 0; i < N; i++) {
            double err = sqrt((x2[i].re - x4[i].re) * (x2[i].re - x4[i].re) +
                              (x2[i].im - x4[i].im) * (x2[i].im - x4[i].im));
            if (err > max_err) max_err = err;
        }

        free(x2);
        free(x4);
        if (max_err < 1e-6) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("radix-4 mismatch"); }
    }

    TEST_CASE_BEGIN("Radix-4 FFT round-trip");
    {
        int N = 64;
        Complex *x = (Complex *)malloc((size_t)N * sizeof(Complex));
        double *orig = (double *)malloc((size_t)N * sizeof(double));

        for (int i = 0; i < N; i++) {
            orig[i] = sin(2.0 * M_PI * 3.0 * (double)i / N);
            x[i].re = orig[i];
            x[i].im = 0.0;
        }

        fft_radix4(x, N);
        ifft_radix4(x, N);

        double max_err = 0;
        for (int i = 0; i < N; i++) {
            double err = fabs(x[i].re - orig[i]);
            if (err > max_err) max_err = err;
        }

        free(x);
        free(orig);
        if (max_err < 1e-12) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("round-trip error too large"); }
    }

    TEST_CASE_BEGIN("Radix-4 fallback for non-power-of-4");
    {
        /* N=128 is power of 2 but not power of 4 */
        int N = 128;
        Complex *x2 = (Complex *)malloc((size_t)N * sizeof(Complex));
        Complex *x4 = (Complex *)malloc((size_t)N * sizeof(Complex));

        for (int i = 0; i < N; i++) {
            x2[i].re = cos(2.0 * M_PI * 10.0 * (double)i / N);
            x2[i].im = 0.0;
            x4[i] = x2[i];
        }

        fft(x2, N);
        fft_radix4(x4, N);  /* should fall back to radix-2 */

        double max_err = 0;
        for (int i = 0; i < N; i++) {
            double err = sqrt((x2[i].re - x4[i].re) * (x2[i].re - x4[i].re) +
                              (x2[i].im - x4[i].im) * (x2[i].im - x4[i].im));
            if (err > max_err) max_err = err;
        }

        free(x2);
        free(x4);
        if (max_err < 1e-12) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("radix-4 fallback failed"); }
    }

    TEST_CASE_BEGIN("Twiddle table FFT matches direct FFT");
    {
        int N = 512;
        TwiddleTable *tt = twiddle_create(N);
        Complex *xd = (Complex *)malloc((size_t)N * sizeof(Complex));
        Complex *xt = (Complex *)malloc((size_t)N * sizeof(Complex));

        for (int i = 0; i < N; i++) {
            double t = (double)i / N;
            xd[i].re = sin(2.0 * M_PI * 7.0 * t) + 0.5 * sin(2.0 * M_PI * 23.0 * t);
            xd[i].im = 0.0;
            xt[i] = xd[i];
        }

        fft(xd, N);
        fft_with_twiddles(xt, N, tt);

        double max_err = 0;
        for (int i = 0; i < N; i++) {
            double err = sqrt((xd[i].re - xt[i].re) * (xd[i].re - xt[i].re) +
                              (xd[i].im - xt[i].im) * (xd[i].im - xt[i].im));
            if (err > max_err) max_err = err;
        }

        free(xd);
        free(xt);
        twiddle_destroy(tt);
        if (max_err < 1e-10) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("twiddle FFT mismatch"); }
    }

    /* ─── Aligned Memory Tests ──────────────────────────── */

    TEST_CASE_BEGIN("Aligned memory is properly aligned");
    {
        int ok = 1;
        int aligns[] = {16, 32, 64};
        for (int a = 0; a < 3; a++) {
            void *ptr = aligned_alloc_dsp(aligns[a], 4096);
            if (!ptr || ((size_t)ptr % (size_t)aligns[a]) != 0)
                ok = 0;
            aligned_free_dsp(ptr);
        }
        if (ok) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("alignment failure"); }
    }

    TEST_CASE_BEGIN("Aligned memory read/write works");
    {
        double *ptr = (double *)aligned_alloc_dsp(32, 256 * (int)sizeof(double));
        int ok = (ptr != NULL);
        if (ok) {
            for (int i = 0; i < 256; i++)
                ptr[i] = (double)i * 1.5;
            for (int i = 0; i < 256; i++) {
                if (fabs(ptr[i] - (double)i * 1.5) > 1e-15) {
                    ok = 0;
                    break;
                }
            }
        }
        aligned_free_dsp(ptr);
        if (ok) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("read/write failure"); }
    }

    /* ─── Benchmark Tests ───────────────────────────────── */

    TEST_CASE_BEGIN("Benchmark produces valid results");
    {
        BenchResult r = bench_fft_radix2(64, 10);
        int ok = (r.n == 64) && (r.runs == 10) &&
                 (r.min_us > 0) && (r.avg_us >= r.min_us) &&
                 (r.max_us >= r.avg_us) && (r.mflops > 0);
        if (ok) { TEST_PASS_STMT; }
        else { TEST_FAIL_STMT("invalid bench result"); }
    }

    TEST_CASE_BEGIN("Bench print does not crash");
    {
        BenchResult r = bench_fft_radix4(64, 5);
        bench_print("test", &r);
        /* If we get here, it didn't crash */
        TEST_PASS_STMT;
    }

    /* ── Summary ──────────────────────────────────────────── */
    printf("\n  ────────────────────────────\n");
    printf("  Results: %d/%d passed", test_passed, test_count);
    if (test_failed > 0)
        printf(", %d FAILED", test_failed);
    printf("\n\n");

    return test_failed > 0 ? 1 : 0;
}

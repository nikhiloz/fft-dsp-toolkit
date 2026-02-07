/**
 * @file 28-real-time-streaming.c
 * @brief Chapter 28 demo — Real-Time Streaming DSP.
 *
 * ── Streaming Pipeline ──────────────────────────────────────────
 *
 *   ┌──────────┐   ┌──────────┐   ┌───────────┐   ┌──────────┐
 *   │  Signal  │──►│  Ring    │──►│  Frame    │──►│ Spectral │
 *   │Generator │   │  Buffer  │   │ Processor │   │ Display  │
 *   └──────────┘   └──────────┘   └───────────┘   └──────────┘
 *                     SPSC FIFO     Overlap FFT     Peak detect
 *
 * Demonstrates:
 *   - Ring buffer: producer/consumer circular FIFO
 *   - Frame processor: overlap-add windowed FFT analysis
 *   - Streaming spectrogram: frequency content over time
 *   - Latency measurement: µs-resolution processing budget
 *   - Real-time constraints: deterministic memory, no malloc
 *
 * Build & run:
 *   make chapters && ./build/bin/ch28
 *
 * Read alongside: chapters/28-real-time-streaming.md
 */

#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "realtime.h"
#include "signal_gen.h"
#include "fft.h"
#include "gnuplot.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define FS          8000.0   /* Sample rate */
#define FRAME_SIZE  256      /* FFT length */
#define HOP_SIZE    128      /* 50% overlap */
#define CHUNK_SIZE  64       /* Simulated "block" from audio driver */

/* ── Section 1: Ring Buffer Demonstration ─────────────────────── */

static void demo_ring_buffer(void)
{
    printf("── Section 1: Ring Buffer (SPSC FIFO) ──\n\n");

    RingBuffer *rb = ring_buffer_create(16);  /* rounds up to 16 */

    printf("  Created ring buffer: capacity=%d\n", rb->cap);
    printf("  Initial: available=%d  space=%d\n\n",
           ring_buffer_available(rb), ring_buffer_space(rb));

    /* Write some data */
    double data[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    int written = ring_buffer_write(rb, data, 8);
    printf("  Wrote %d samples: available=%d  space=%d\n",
           written, ring_buffer_available(rb), ring_buffer_space(rb));

    /* Peek without consuming */
    double peek_buf[4];
    int peeked = ring_buffer_peek(rb, peek_buf, 4);
    printf("  Peeked %d samples: [%.0f, %.0f, %.0f, %.0f]\n",
           peeked, peek_buf[0], peek_buf[1], peek_buf[2], peek_buf[3]);
    printf("  After peek: available=%d (unchanged)\n",
           ring_buffer_available(rb));

    /* Read (consume) some data */
    double read_buf[4];
    int nread = ring_buffer_read(rb, read_buf, 4);
    printf("  Read %d samples:   [%.0f, %.0f, %.0f, %.0f]\n",
           nread, read_buf[0], read_buf[1], read_buf[2], read_buf[3]);
    printf("  After read: available=%d  space=%d\n",
           ring_buffer_available(rb), ring_buffer_space(rb));

    /* Skip some data */
    int skipped = ring_buffer_skip(rb, 2);
    printf("  Skipped %d samples: available=%d\n",
           skipped, ring_buffer_available(rb));

    /* Write more (wraps around) */
    double more[] = {9.0, 10.0, 11.0, 12.0, 13.0, 14.0};
    written = ring_buffer_write(rb, more, 6);
    printf("  Wrote %d more (wrapping): available=%d  space=%d\n",
           written, ring_buffer_available(rb), ring_buffer_space(rb));

    /* Read all remaining */
    double out[16];
    nread = ring_buffer_read(rb, out, 16);
    printf("  Read all %d remaining: [", nread);
    for (int i = 0; i < nread; i++)
        printf("%.0f%s", out[i], i < nread - 1 ? " " : "");
    printf("]\n\n");

    ring_buffer_destroy(rb);
}

/* ── Section 2: Streaming FFT with Frame Processor ────────────── */

static void demo_streaming_fft(void)
{
    printf("── Section 2: Streaming FFT Analysis ──\n\n");
    printf("  Frame size: %d samples (%.1f ms)\n", FRAME_SIZE, 1000.0 * FRAME_SIZE / FS);
    printf("  Hop size:   %d samples (%.1f ms, %.0f%% overlap)\n",
           HOP_SIZE, 1000.0 * HOP_SIZE / FS,
           100.0 * (1.0 - (double)HOP_SIZE / FRAME_SIZE));
    printf("  Chunk size: %d samples (simulated audio blocks)\n\n", CHUNK_SIZE);

    RingBuffer *rb = ring_buffer_create(4096);
    FrameProcessor *fp = frame_processor_create(FRAME_SIZE, HOP_SIZE);
    LatencyStats lat;
    latency_init(&lat);

    /* Generate a test signal: chirp from 200 Hz → 2000 Hz over 1 second */
    int total_samples = (int)FS;  /* 1 second */
    double *signal = (double *)malloc((size_t)total_samples * sizeof(double));
    for (int i = 0; i < total_samples; i++) {
        double t = (double)i / FS;
        /* Linear chirp: f(t) = f0 + (f1 - f0) * t / T */
        double f0 = 200.0, f1 = 2000.0, T = 1.0;
        double freq = f0 + (f1 - f0) * t / T;
        double phase = 2.0 * M_PI * (f0 * t + 0.5 * (f1 - f0) * t * t / T);
        signal[i] = sin(phase);
        (void)freq;  /* Used conceptually */
    }

    printf("  Signal: linear chirp 200 → 2000 Hz, %.0f Hz sample rate\n\n", FS);
    printf("  Streaming spectrogram (time → frequency):\n");
    printf("  ─────────────────────────────────────────\n");

    int sample_pos = 0;
    int frame_count = 0;

    while (sample_pos < total_samples) {
        /* Simulate: producer writes a chunk to ring buffer */
        int chunk = CHUNK_SIZE;
        if (sample_pos + chunk > total_samples)
            chunk = total_samples - sample_pos;

        ring_buffer_write(rb, signal + sample_pos, chunk);
        sample_pos += chunk;

        /* Consumer: read from ring buffer and feed to frame processor */
        double block[CHUNK_SIZE];
        int got = ring_buffer_read(rb, block, chunk);

        double t0 = timer_usec();
        int frames = frame_processor_feed(fp, block, got);
        double t1 = timer_usec();

        if (frames > 0) {
            latency_record(&lat, t1 - t0);

            /* Print spectrogram row */
            double peak_hz = frame_processor_peak_freq(fp, FS);
            double time_sec = (double)sample_pos / FS;

            printf("  t=%5.3fs | peak=%7.1f Hz | ", time_sec, peak_hz);

            /* Simple ASCII spectrogram: show energy in frequency bands */
            int nbands = 16;
            int half = FRAME_SIZE / 2;
            int bins_per_band = half / nbands;
            for (int b = 0; b < nbands; b++) {
                double sum = 0;
                for (int k = b * bins_per_band; k < (b + 1) * bins_per_band; k++)
                    sum += fp->magnitude[k];
                sum /= bins_per_band;

                char ch = ' ';
                if (sum > 0.30) ch = '#';
                else if (sum > 0.15) ch = '=';
                else if (sum > 0.05) ch = '-';
                else if (sum > 0.01) ch = '.';
                printf("%c", ch);
            }
            printf("|\n");

            frame_count++;
            if (frame_count >= 20) {
                printf("  ... (showing first 20 frames)\n");
                break;
            }
        }
    }

    printf("\n  Processed %d frames total\n", fp->frames_processed);
    printf("  Latency: min=%.1f µs  avg=%.1f µs  max=%.1f µs\n",
           lat.min_us, latency_avg(&lat), lat.max_us);

    /* Budget check */
    double budget_us = 1e6 * HOP_SIZE / FS;
    printf("  Budget per frame: %.0f µs (hop=%d at %.0f Hz)\n", budget_us, HOP_SIZE, FS);
    printf("  Margin: %.1f%% of budget used\n",
           100.0 * latency_avg(&lat) / budget_us);

    free(signal);
    frame_processor_destroy(fp);
    ring_buffer_destroy(rb);
    printf("\n");
}

/* ── Section 3: Latency Measurement ──────────────────────────── */

static void demo_latency(void)
{
    printf("── Section 3: Latency Budget Analysis ──\n\n");

    int sizes[] = {64, 128, 256, 512, 1024};
    int nsizes = 5;
    int runs = 50;

    printf("  FFT Frame Processing Latency (Hann window + FFT + magnitude):\n\n");
    printf("  Frame Size    Budget(µs)    Avg Latency    Margin\n");
    printf("  ─────────────────────────────────────────────────\n");

    for (int s = 0; s < nsizes; s++) {
        int N = sizes[s];
        FrameProcessor *fp = frame_processor_create(N, N / 2);
        LatencyStats lat;
        latency_init(&lat);

        /* Generate N samples */
        double *sig = (double *)calloc((size_t)N, sizeof(double));
        for (int i = 0; i < N; i++)
            sig[i] = sin(2.0 * M_PI * 440.0 * i / FS);

        for (int r = 0; r < runs; r++) {
            double t0 = timer_usec();
            frame_processor_feed(fp, sig, N);
            double t1 = timer_usec();
            latency_record(&lat, t1 - t0);
        }

        double budget = 1e6 * (N / 2) / FS;  /* hop = N/2 */
        double avg = latency_avg(&lat);
        double pct = 100.0 * avg / budget;

        printf("  %5d         %7.0f       %7.1f µs     %5.1f%%  %s\n",
               N, budget, avg, pct,
               pct < 50 ? "[OK]" : pct < 80 ? "[WARN]" : "[TIGHT]");

        free(sig);
        frame_processor_destroy(fp);
    }
    printf("\n");
}

/* ── Section 4: Plots ────────────────────────────────────────── */

static void demo_plots(void)
{
    printf("── Section 4: Generating Plots ──\n\n");

    gp_init("ch28");

    /* Plot 1: Ring buffer fill level over time */
    {
        RingBuffer *rb = ring_buffer_create(256);
        double x[100], y[100];
        int idx = 0;

        /* Simulate: bursty writes, steady reads */
        for (int t = 0; t < 100; t++) {
            if (t % 3 == 0) {
                double data[20];
                for (int i = 0; i < 20; i++) data[i] = 0;
                ring_buffer_write(rb, data, 20);
            }
            double buf[8];
            ring_buffer_read(rb, buf, 8);

            x[idx] = (double)t;
            y[idx] = 100.0 * ring_buffer_available(rb) / (rb->cap - 1);
            idx++;
        }

        gp_plot_1("ch28", "ring_buffer_fill", "Ring Buffer Fill Level",
                  "Time (chunks)", "Fill %", x, y, idx, "lines");
        printf("  [saved] plots/ch28/ring_buffer_fill.png\n");
        ring_buffer_destroy(rb);
    }

    /* Plot 2: Streaming spectrogram (frequency vs time) */
    {
        FrameProcessor *fp = frame_processor_create(FRAME_SIZE, HOP_SIZE);
        double time_arr[200], freq_arr[200];
        int npts = 0;

        int total = (int)FS;
        double *sig = (double *)malloc((size_t)total * sizeof(double));
        for (int i = 0; i < total; i++) {
            double t = (double)i / FS;
            double phase = 2.0 * M_PI * (200.0 * t + 0.5 * 1800.0 * t * t);
            sig[i] = sin(phase);
        }

        for (int pos = 0; pos + CHUNK_SIZE <= total && npts < 200; pos += CHUNK_SIZE) {
            int frames = frame_processor_feed(fp, sig + pos, CHUNK_SIZE);
            if (frames > 0) {
                time_arr[npts] = (double)pos / FS;
                freq_arr[npts] = frame_processor_peak_freq(fp, FS);
                npts++;
            }
        }

        gp_plot_1("ch28", "streaming_spectrogram",
                  "Streaming Spectrogram (Chirp 200-2000 Hz)",
                  "Time (s)", "Peak Frequency (Hz)",
                  time_arr, freq_arr, npts, "linespoints");
        printf("  [saved] plots/ch28/streaming_spectrogram.png\n");
        free(sig);
        frame_processor_destroy(fp);
    }

    printf("\n");
}

/* ── Main ────────────────────────────────────────────────────── */

int main(void)
{
    printf("=== Chapter 28: Real-Time Streaming DSP ===\n\n");

    demo_ring_buffer();
    demo_streaming_fft();
    demo_latency();
    demo_plots();

    printf("Key takeaways:\n");
    printf("  1. Ring buffers decouple producer/consumer rates\n");
    printf("  2. Frame processors provide overlap-add spectral analysis\n");
    printf("  3. Latency budgets must be checked against real-time constraints\n");
    printf("  4. Pre-allocation eliminates malloc jitter in the processing loop\n");

    return 0;
}

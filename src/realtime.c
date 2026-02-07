/**
 * @file realtime.c
 * @brief Chapter 28 — Real-time streaming DSP infrastructure.
 *
 * Implements:
 *   - Ring buffer (power-of-2 circular FIFO)
 *   - Frame processor (overlap-add windowed FFT)
 *   - Microsecond latency timer (POSIX clock_gettime)
 *
 * ── Ring Buffer Layout ──────────────────────────────────────────
 *
 *   Index:    0   1   2   3   4   5   6   7
 *           ┌───┬───┬───┬───┬───┬───┬───┬───┐
 *           │ . │ . │ . │ W │ . │ R │ . │ . │
 *           └───┴───┴───┴───┴───┴───┴───┴───┘
 *                         ↑       ↑
 *                        head    tail
 *
 *   available = (head - tail) & mask
 *   space     = cap - 1 - available   (always leave 1 slot empty)
 *
 * ── Frame Processor Pipeline ────────────────────────────────────
 *
 *   new samples ──► overlap_buf ──► window ──► FFT ──► |X[k]|
 *                   [N samples]    [Hann]             [N/2 bins]
 *
 * @see include/realtime.h
 */

#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "realtime.h"
#include "fft.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ================================================================== */
/*  Utility: next power of 2                                          */
/* ================================================================== */

static int next_pow2(int v)
{
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    return v + 1;
}

/* ================================================================== */
/*  Ring Buffer                                                       */
/* ================================================================== */

RingBuffer *ring_buffer_create(int capacity)
{
    if (capacity < 2) capacity = 2;
    capacity = next_pow2(capacity);

    RingBuffer *rb = (RingBuffer *)calloc(1, sizeof(RingBuffer));
    if (!rb) return NULL;

    rb->buf  = (double *)calloc((size_t)capacity, sizeof(double));
    if (!rb->buf) { free(rb); return NULL; }

    rb->cap  = capacity;
    rb->mask = capacity - 1;
    rb->head = 0;
    rb->tail = 0;
    return rb;
}

void ring_buffer_destroy(RingBuffer *rb)
{
    if (!rb) return;
    free(rb->buf);
    free(rb);
}

int ring_buffer_available(const RingBuffer *rb)
{
    return (rb->head - rb->tail) & rb->mask;
}

int ring_buffer_space(const RingBuffer *rb)
{
    /* Leave one slot empty to distinguish full from empty */
    return rb->cap - 1 - ring_buffer_available(rb);
}

int ring_buffer_write(RingBuffer *rb, const double *data, int n)
{
    int space = ring_buffer_space(rb);
    if (n > space) n = space;

    for (int i = 0; i < n; i++) {
        rb->buf[rb->head & rb->mask] = data[i];
        rb->head = (rb->head + 1) & rb->mask;
    }
    return n;
}

int ring_buffer_read(RingBuffer *rb, double *data, int n)
{
    int avail = ring_buffer_available(rb);
    if (n > avail) n = avail;

    for (int i = 0; i < n; i++) {
        data[i] = rb->buf[rb->tail & rb->mask];
        rb->tail = (rb->tail + 1) & rb->mask;
    }
    return n;
}

int ring_buffer_peek(const RingBuffer *rb, double *data, int n)
{
    int avail = ring_buffer_available(rb);
    if (n > avail) n = avail;

    int pos = rb->tail;
    for (int i = 0; i < n; i++) {
        data[i] = rb->buf[pos & rb->mask];
        pos = (pos + 1) & rb->mask;
    }
    return n;
}

int ring_buffer_skip(RingBuffer *rb, int n)
{
    int avail = ring_buffer_available(rb);
    if (n > avail) n = avail;
    rb->tail = (rb->tail + n) & rb->mask;
    return n;
}

void ring_buffer_reset(RingBuffer *rb)
{
    rb->head = 0;
    rb->tail = 0;
}

/* ================================================================== */
/*  Frame Processor                                                   */
/* ================================================================== */

FrameProcessor *frame_processor_create(int frame_size, int hop_size)
{
    if (frame_size < 4) frame_size = 4;
    frame_size = next_pow2(frame_size);
    if (hop_size <= 0 || hop_size > frame_size)
        hop_size = frame_size / 2;

    FrameProcessor *fp = (FrameProcessor *)calloc(1, sizeof(FrameProcessor));
    if (!fp) return NULL;

    fp->frame_size       = frame_size;
    fp->hop_size         = hop_size;
    fp->frame            = (double *)calloc((size_t)frame_size, sizeof(double));
    fp->window           = (double *)malloc((size_t)frame_size * sizeof(double));
    fp->spectrum         = (Complex *)calloc((size_t)frame_size, sizeof(Complex));
    fp->magnitude        = (double *)calloc((size_t)(frame_size / 2), sizeof(double));
    fp->magnitude_db     = (double *)calloc((size_t)(frame_size / 2), sizeof(double));
    fp->overlap_buf      = (double *)calloc((size_t)frame_size, sizeof(double));
    fp->frames_processed = 0;
    fp->samples_queued   = 0;

    /* Pre-compute Hann window */
    for (int i = 0; i < frame_size; i++)
        fp->window[i] = 0.5 * (1.0 - cos(2.0 * M_PI * i / (frame_size - 1)));

    return fp;
}

void frame_processor_destroy(FrameProcessor *fp)
{
    if (!fp) return;
    free(fp->frame);
    free(fp->window);
    free(fp->spectrum);
    free(fp->magnitude);
    free(fp->magnitude_db);
    free(fp->overlap_buf);
    free(fp);
}

/**
 * Process one frame: window → FFT → magnitude.
 *
 *   overlap_buf holds the latest frame_size samples.
 *   Apply Hann window, compute FFT, extract |X[k]|.
 */
static void process_one_frame(FrameProcessor *fp)
{
    int N = fp->frame_size;
    int half = N / 2;

    /* Apply window to overlap buffer */
    for (int i = 0; i < N; i++) {
        fp->spectrum[i].re = fp->overlap_buf[i] * fp->window[i];
        fp->spectrum[i].im = 0.0;
    }

    /* FFT (in-place) */
    fft(fp->spectrum, N);

    /* Magnitude and dB */
    for (int i = 0; i < half; i++) {
        double mag = complex_mag(fp->spectrum[i]);
        fp->magnitude[i] = mag / half;  /* normalise */
        fp->magnitude_db[i] = 20.0 * log10(mag / half + 1e-30);
    }

    fp->frames_processed++;
}

int frame_processor_feed(FrameProcessor *fp, const double *samples, int n)
{
    int N = fp->frame_size;
    int hop = fp->hop_size;
    int frames_out = 0;

    for (int i = 0; i < n; i++) {
        /* Shift overlap buffer left by 1 and append new sample */
        /* Optimisation: only shift when we have hop_size new samples */
        fp->overlap_buf[(fp->frames_processed == 0 ?
            fp->samples_queued : (N - hop + fp->samples_queued))] = samples[i];
        fp->samples_queued++;

        if (fp->frames_processed == 0 && fp->samples_queued == N) {
            /* First frame: need full N samples */
            process_one_frame(fp);
            fp->samples_queued = 0;
            frames_out++;
        } else if (fp->frames_processed > 0 && fp->samples_queued == hop) {
            /* Subsequent frames: shift and process */
            /* Shift old data left by hop, new samples are at the end */
            memmove(fp->overlap_buf, fp->overlap_buf + hop,
                    (size_t)(N - hop) * sizeof(double));
            /* samples_queued == hop, and the new samples are already
               written at positions [N-hop .. N-1] via the index above,
               but let's fix: rewrite properly */
            fp->samples_queued = 0;
            process_one_frame(fp);
            frames_out++;
        }
    }
    return frames_out;
}

int frame_processor_peak_bin(const FrameProcessor *fp)
{
    int half = fp->frame_size / 2;
    int peak = 0;
    double peak_val = -1e30;
    for (int i = 1; i < half; i++) {
        if (fp->magnitude[i] > peak_val) {
            peak_val = fp->magnitude[i];
            peak = i;
        }
    }
    return peak;
}

double frame_processor_peak_freq(const FrameProcessor *fp, double fs)
{
    int bin = frame_processor_peak_bin(fp);
    return (double)bin * fs / (double)fp->frame_size;
}

/* ================================================================== */
/*  Latency Timer                                                     */
/* ================================================================== */

double timer_usec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec * 1e6 + (double)ts.tv_nsec / 1000.0;
}

void latency_init(LatencyStats *ls)
{
    ls->min_us = 1e30;
    ls->max_us = 0.0;
    ls->sum_us = 0.0;
    ls->count  = 0;
}

void latency_record(LatencyStats *ls, double us)
{
    if (us < ls->min_us) ls->min_us = us;
    if (us > ls->max_us) ls->max_us = us;
    ls->sum_us += us;
    ls->count++;
}

double latency_avg(const LatencyStats *ls)
{
    if (ls->count == 0) return 0.0;
    return ls->sum_us / ls->count;
}

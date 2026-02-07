/**
 * @file realtime.h
 * @brief Chapter 28 — Real-time streaming DSP infrastructure.
 *
 * Provides the building blocks for real-time DSP pipelines:
 *
 *   ┌──────────┐    ┌──────────┐    ┌────────────┐    ┌──────────┐
 *   │  Audio   │───►│  Ring    │───►│  Frame     │───►│  Result  │
 *   │  Source  │    │  Buffer  │    │ Processor  │    │  Buffer  │
 *   └──────────┘    └──────────┘    └────────────┘    └──────────┘
 *
 * Components:
 *   - **RingBuffer**: Lock-free SPSC circular FIFO (power-of-2 capacity)
 *   - **FrameProcessor**: Overlap-add frame extraction + windowed FFT
 *   - **LatencyTimer**: Microsecond-resolution timing for budget tracking
 *
 * All structures are pre-allocated — zero malloc during processing.
 *
 * @note Uses POSIX clock_gettime for timing. Define _POSIX_C_SOURCE
 *       200809L before including <time.h>.
 *
 * @see chapters/28-real-time-streaming.md
 */

#ifndef REALTIME_H
#define REALTIME_H

#include "dsp_utils.h"

/* ================================================================== */
/*  Ring Buffer — lock-free SPSC circular FIFO                        */
/* ================================================================== */

/**
 * @brief Circular buffer with power-of-2 capacity.
 *
 * Uses bitwise AND mask for O(1) index wrapping.
 * Safe for single-producer / single-consumer without locks.
 *
 *   write ──►  [5][6][7][_][_][_][1][2][3][4]  ◄── read
 *               ↑                               ↑
 *              head                             tail
 */
typedef struct {
    double *buf;     /**< Pre-allocated sample storage */
    int     cap;     /**< Capacity (must be power of 2) */
    int     mask;    /**< cap - 1 for bitwise wrap */
    int     head;    /**< Write position */
    int     tail;    /**< Read position */
} RingBuffer;

/** Create a ring buffer.  capacity is rounded up to next power of 2. */
RingBuffer *ring_buffer_create(int capacity);

/** Destroy a ring buffer and free memory. */
void ring_buffer_destroy(RingBuffer *rb);

/** Write up to n samples.  Returns number actually written. */
int ring_buffer_write(RingBuffer *rb, const double *data, int n);

/** Read up to n samples.  Returns number actually read. */
int ring_buffer_read(RingBuffer *rb, double *data, int n);

/** Number of samples available for reading. */
int ring_buffer_available(const RingBuffer *rb);

/** Free space available for writing. */
int ring_buffer_space(const RingBuffer *rb);

/** Peek at n samples without consuming them.  Returns count peeked. */
int ring_buffer_peek(const RingBuffer *rb, double *data, int n);

/** Discard up to n samples from the read side. */
int ring_buffer_skip(RingBuffer *rb, int n);

/** Reset to empty state. */
void ring_buffer_reset(RingBuffer *rb);

/* ================================================================== */
/*  Frame Processor — overlap windowed FFT analysis                   */
/* ================================================================== */

/**
 * @brief Real-time frame processor for streaming FFT analysis.
 *
 * Extracts overlapping frames from a ring buffer, applies a window,
 * computes the FFT, and stores the magnitude spectrum.
 *
 *   ┌─────────────────────────────────┐
 *   │ Frame 0:  [=========]           │  frame_size = 8
 *   │ Frame 1:      [=========]       │  hop_size   = 4  (50% overlap)
 *   │ Frame 2:          [=========]   │
 *   └─────────────────────────────────┘
 */
typedef struct {
    int      frame_size;       /**< FFT length N (power of 2) */
    int      hop_size;         /**< Samples between frames */
    double  *frame;            /**< Current windowed frame [N] */
    double  *window;           /**< Pre-computed window [N] */
    Complex *spectrum;         /**< FFT output [N] */
    double  *magnitude;        /**< |X[k]| for k=0..N/2-1 */
    double  *magnitude_db;     /**< 20·log10(|X[k]|) */
    int      frames_processed; /**< Counter */
    int      samples_queued;   /**< Samples received since last frame */
    double  *overlap_buf;      /**< Internal overlap storage [N] */
} FrameProcessor;

/** Create a frame processor.  frame_size must be power of 2. */
FrameProcessor *frame_processor_create(int frame_size, int hop_size);

/** Destroy and free all internal buffers. */
void frame_processor_destroy(FrameProcessor *fp);

/**
 * @brief Feed samples into the processor.
 *
 * Internally accumulates samples.  When enough for a new frame, it:
 *   1. Extracts a frame (with overlap from previous)
 *   2. Applies the Hann window
 *   3. Computes the FFT
 *   4. Stores magnitude + dB in fp->magnitude / fp->magnitude_db
 *
 * @return Number of frames produced (0 or 1+).
 */
int frame_processor_feed(FrameProcessor *fp, const double *samples, int n);

/** Get the dominant frequency bin from the most recent frame. */
int frame_processor_peak_bin(const FrameProcessor *fp);

/** Get the dominant frequency in Hz given sample rate fs. */
double frame_processor_peak_freq(const FrameProcessor *fp, double fs);

/* ================================================================== */
/*  Latency Timer — microsecond-resolution timing                     */
/* ================================================================== */

/** Return current time in microseconds (monotonic clock). */
double timer_usec(void);

/**
 * @brief Latency statistics accumulator.
 *
 * Call timer_usec() before and after processing, then
 * latency_record(lat, after - before) to track min/max/avg.
 */
typedef struct {
    double min_us;   /**< Minimum latency observed */
    double max_us;   /**< Maximum latency observed */
    double sum_us;   /**< Running sum for average */
    int    count;    /**< Number of measurements */
} LatencyStats;

/** Initialise latency stats (min = huge, max = 0). */
void latency_init(LatencyStats *ls);

/** Record one latency measurement in microseconds. */
void latency_record(LatencyStats *ls, double us);

/** Get average latency.  Returns 0 if no measurements. */
double latency_avg(const LatencyStats *ls);

#endif /* REALTIME_H */

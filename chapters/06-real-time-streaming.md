# Chapter 6 — Real-Time Streaming

Real-time DSP processes audio or sensor data as it arrives, with strict
latency constraints. This chapter covers the architecture and
algorithms used to do continuous FFT analysis on a live stream.

> **Status:** This chapter describes the design. The streaming
> module (`ring_buffer.c`, `streaming.c`) is planned for Phase 3 of
> the project.

---

## 6.1 The Challenge

Batch processing (load file → FFT → done) is simple. Real-time adds
three constraints:

1. **Fixed latency** — processing must complete before the next frame
   arrives
2. **No allocation** — `malloc` during processing causes unpredictable
   delays
3. **No data loss** — every sample must be processed

> **📊 Real-Time Architecture** — [View full-size diagram →](../reference/diagrams/realtime_architecture.png)

## 6.2 The Overlap-Add Method

To do continuous FFT analysis, we process the signal in overlapping
frames:

```
Frame 0:  [═══════════════]
Frame 1:        [═══════════════]       ← 50% overlap
Frame 2:              [═══════════════]
```

Each frame:
1. Copy $N$ samples from the input buffer (50% new, 50% from previous)
2. Apply window function
3. Compute FFT
4. Extract features (magnitude, phase, etc.)
5. For reconstruction: IFFT → overlap-add with previous frame

### Why Overlap?

Windowing attenuates the signal edges. With 50% overlap, every sample
falls near the centre of at least one window, so no energy is lost.
The Hann window with 50% overlap gives perfect reconstruction:

$$\sum_{m} w[n - mH] = 1 \quad \text{(constant for all } n \text{)}$$

where $H = N/2$ is the hop size.

## 6.3 Ring Buffer Design

A **ring buffer** (circular FIFO) is the standard data structure for
real-time streaming:

```
   write ──→  [5][6][7][_][_][_][1][2][3][4]  ←── read
              ↑                               ↑
              head                            tail
```

Properties:
- Fixed-size, pre-allocated memory
- **Lock-free** for single-producer / single-consumer
- O(1) read and write
- Wraps around — no copying needed

### Planned API

```c
typedef struct {
    double *buffer;
    int     capacity;    /* must be power of 2 */
    int     head;        /* write position (atomic) */
    int     tail;        /* read position (atomic) */
} RingBuffer;

RingBuffer *ring_buffer_create(int capacity);
void        ring_buffer_destroy(RingBuffer *rb);
int         ring_buffer_write(RingBuffer *rb, const double *data, int n);
int         ring_buffer_read(RingBuffer *rb, double *data, int n);
int         ring_buffer_available(const RingBuffer *rb);
```

The capacity must be a power of 2 so that modular indexing can use
bitwise AND: `index & (capacity - 1)` instead of `index % capacity`.

## 6.4 Latency Budget

For audio at 48000 Hz with 1024-point FFT:

| Stage | Time |
|-------|------|
| Fill buffer (1024 samples) | 21.3 ms |
| Window + FFT | ~0.1 ms |
| Feature extraction | ~0.01 ms |
| **Total frame latency** | **~21.4 ms** |

With 50% overlap (hop size 512):
- New data arrives every $512/48000 = 10.7$ ms
- Processing must complete in <10.7 ms to keep up

For music applications, <20 ms total latency is acceptable. For live
monitoring, <5 ms is better.

## 6.5 Memory Pre-Allocation

In real-time code, never call `malloc` in the processing loop:

```c
/* BAD — malloc can take microseconds to milliseconds */
double *frame = malloc(N * sizeof(double));

/* GOOD — allocate once at startup */
typedef struct {
    double     frame[1024];
    Complex    spectrum[1024];
    double     magnitudes[1024];
    RingBuffer input_ring;
} StreamContext;

StreamContext *ctx = create_stream_context();  /* once at init */
```

All temporary arrays should be part of a pre-allocated context struct.

## 6.6 Processing Thread Architecture

```
Audio Input Thread          Processing Thread
────────────────            ─────────────────
  ALSA/sensor read     →     ring_buffer_read
                              apply_window
                              fft
                              fft_magnitude
  ring_buffer_write    ←     (results available)
```

The two threads communicate via the lock-free ring buffer. The
processing thread:
1. Waits until $N/2$ new samples appear in the ring
2. Reads $N$ samples (50% overlap with previous frame)
3. Processes the frame
4. Writes results to an output ring buffer
5. Repeats

## 6.7 Connection to Existing Code

The streaming module will use:
- `apply_window()` from [`src/dsp_utils.c`](../src/dsp_utils.c) — Chapter 03
- `fft()` from [`src/fft.c`](../src/fft.c) — Chapter 02
- `fft_magnitude()` from [`src/fft.c`](../src/fft.c) — Chapter 02
- `fir_filter()` from [`src/filter.c`](../src/filter.c) — Chapter 04

Every function we've built so far is a building block for the real-time
pipeline.

## 6.8 Exercises

1. **Code exercise:** Implement a basic `RingBuffer` using the API
   above. Make sure `write` and `read` use `& (capacity - 1)` for
   wrap-around.

2. **Thinking question:** Why must the ring buffer capacity be a power
   of 2?

3. **Design exercise:** Sketch the memory layout for a streaming FFT
   that processes 1024-point frames with 50% overlap. How many samples
   of history must be kept?

4. **Advanced:** Implement a simple stdin-based streaming demo that
   reads raw PCM samples from a pipe and prints the dominant frequency
   every frame.

---

**Previous:** [Chapter 05 — Spectral Analysis](05-spectral-analysis.md)
| **Next:** [Chapter 07 — Optimisation →](07-optimisation.md)

/**
 * @file cepstrum.h
 * @brief Cepstral analysis — real/complex cepstrum, liftering, MFCCs.
 *
 * Chapter 26: Cepstrum Analysis and Mel-Frequency Cepstral Coefficients
 *
 * Dependencies: fft (dsp_utils.h / fft.h), math.h
 */

#ifndef CEPSTRUM_H
#define CEPSTRUM_H

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------ */
/*  Real / Complex Cepstrum                                            */
/* ------------------------------------------------------------------ */

/**
 * Real cepstrum:  c[n] = IFFT{ log|FFT{x}| }
 *
 * @param x     Input signal (length n)
 * @param n     Signal length (need NOT be power of 2; zero-padded to nfft)
 * @param c     Output real cepstrum (length nfft)
 * @param nfft  FFT size (must be power of 2 and >= n)
 */
void cepstrum_real(const double *x, int n, double *c, int nfft);

/**
 * Complex cepstrum:  c[n] = IFFT{ log(FFT{x}) }
 *
 * Uses unwrapped phase. Returns real part of complex cepstrum.
 *
 * @param x     Input signal (length n)
 * @param n     Signal length
 * @param c     Output complex cepstrum (length nfft)
 * @param nfft  FFT size (power of 2, >= n)
 */
void cepstrum_complex(const double *x, int n, double *c, int nfft);

/**
 * Low-time liftering — extract spectral envelope.
 * Keeps cepstral coefficients 0..L-1, zeros the rest, then FFT back.
 *
 * @param c         Input cepstrum (length nfft)
 * @param nfft      Cepstrum/FFT length
 * @param L         Lifter length (how many low-quefrency coeffs to keep)
 * @param envelope  Output log spectral envelope (length nfft/2, in dB)
 */
void cepstrum_lifter(const double *c, int nfft, int L, double *envelope);

/* ------------------------------------------------------------------ */
/*  Mel-Frequency Cepstral Coefficients (MFCC)                         */
/* ------------------------------------------------------------------ */

/**
 * Convert frequency in Hz to Mel scale.
 * m = 2595 · log10(1 + f/700)
 */
double hz_to_mel(double f_hz);

/**
 * Convert Mel to Hz.
 * f = 700 · (10^(m/2595) - 1)
 */
double mel_to_hz(double mel);

/**
 * Compute Mel filterbank energies.
 *
 * @param power_spec    Power spectrum |X[k]|² (length nfft/2)
 * @param nfft          FFT size
 * @param fs            Sampling rate (Hz)
 * @param n_filters     Number of Mel filters
 * @param fbank         Output filterbank energies (length n_filters)
 * @param f_low         Lowest frequency (Hz), typically 0 or 20
 * @param f_high        Highest frequency (Hz), typically fs/2
 */
void mel_filterbank(const double *power_spec, int nfft, double fs,
                    int n_filters, double *fbank,
                    double f_low, double f_high);

/**
 * Compute MFCCs from a windowed frame.
 *
 * Pipeline: Hamming window → FFT → |X|² → Mel filterbank → log → DCT-II
 *
 * @param frame     Input time-domain frame (length frame_len)
 * @param frame_len Frame size in samples
 * @param nfft      FFT size (power of 2, >= frame_len)
 * @param fs        Sampling rate (Hz)
 * @param n_filters Number of Mel filters (e.g. 26)
 * @param n_mfcc    Number of MFCCs to return (e.g. 13)
 * @param mfcc      Output MFCC coefficients (length n_mfcc)
 */
void compute_mfcc(const double *frame, int frame_len, int nfft,
                  double fs, int n_filters, int n_mfcc, double *mfcc);

/**
 * Type-II DCT (unitary normalisation).
 *
 * Y[k] = Σ_{n=0}^{N-1} x[n] · cos(π(2n+1)k / 2N)
 *
 * @param x     Input (length n)
 * @param y     Output DCT coefficients (length n)
 * @param n     Length
 */
void dct_ii(const double *x, double *y, int n);

#ifdef __cplusplus
}
#endif

#endif /* CEPSTRUM_H */

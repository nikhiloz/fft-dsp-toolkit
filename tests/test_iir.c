/**
 * @file test_iir.c
 * @brief Unit tests for IIR filter library — design, processing, analysis.
 *
 * Tests cover:
 *   1. Biquad DF1 identity (b0=1, rest=0 → passthrough)
 *   2. Biquad DF1 vs DF2T produce identical output
 *   3. Butterworth lowpass DC gain = 1 (0 dB)
 *   4. Butterworth lowpass -3 dB at cutoff frequency
 *   5. Butterworth lowpass attenuates in stopband
 *   6. Butterworth highpass blocks DC
 *   7. Chebyshev Type I steeper than Butterworth (same order)
 *   8. SOS cascade impulse response matches iir_filter expansion
 *   9. All designed filters are stable (poles inside unit circle)
 *  10. Group delay of symmetric FIR is constant
 *
 * Run with: make test
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "test_framework.h"
#include "iir.h"
#include "dsp_utils.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void)
{
    TEST_SUITE("IIR Filter Functions");

    /* ── Test 1: Biquad identity passthrough ─────────────────── */
    TEST_CASE_BEGIN("Biquad identity (passthrough)");
    {
        Biquad bq = {1.0, 0.0, 0.0, 0.0, 0.0};
        BiquadDF1State st;
        biquad_df1_init(&st);

        double input[] = {1.0, 2.0, 3.0, -1.0, 0.5};
        int ok = 1;
        for (int i = 0; i < 5; i++) {
            double y = biquad_process_df1(&bq, &st, input[i]);
            if (fabs(y - input[i]) > 1e-10) { ok = 0; break; }
        }
        if (ok) { TEST_PASS_STMT; }
        else    { TEST_FAIL_STMT("b0=1, rest=0 should pass through unchanged"); }
    }

    /* ── Test 2: DF1 matches DF2T ────────────────────────────── */
    TEST_CASE_BEGIN("DF1 matches DF2T (same biquad)");
    {
        Biquad bq = {0.2, 0.4, 0.2, -0.5, 0.1};
        BiquadDF1State  df1;
        BiquadDF2TState df2t;
        biquad_df1_init(&df1);
        biquad_df2t_init(&df2t);

        int ok = 1;
        for (int i = 0; i < 50; i++) {
            double x = (i == 0) ? 1.0 : 0.0;
            double y1 = biquad_process_df1(&bq, &df1, x);
            double y2 = biquad_process_df2t(&bq, &df2t, x);
            if (fabs(y1 - y2) > 1e-10) { ok = 0; break; }
        }
        if (ok) { TEST_PASS_STMT; }
        else    { TEST_FAIL_STMT("DF1 and DF2T should give identical output"); }
    }

    /* ── Test 3: Butterworth lowpass DC gain = 1 ─────────────── */
    TEST_CASE_BEGIN("Butterworth LP: unity DC gain");
    {
        int orders[] = {1, 2, 4, 6, 8};
        int ok = 1;
        for (int oi = 0; oi < 5; oi++) {
            SOSCascade sos;
            butterworth_lowpass(orders[oi], 0.2, &sos);

            /* DC gain = product of all section DC gains × overall gain */
            double dc = sos.gain;
            for (int s = 0; s < sos.n_sections; s++) {
                Biquad *bq = &sos.sections[s];
                double sec_dc = (bq->b0 + bq->b1 + bq->b2) /
                                (1.0 + bq->a1 + bq->a2);
                dc *= sec_dc;
            }
            if (fabs(dc - 1.0) > 0.01) {
                ok = 0;
                break;
            }
        }
        if (ok) { TEST_PASS_STMT; }
        else    { TEST_FAIL_STMT("DC gain should be 1.0 for all orders"); }
    }

    /* ── Test 4: Butterworth -3 dB at cutoff ─────────────────── */
    TEST_CASE_BEGIN("Butterworth LP: -3 dB at cutoff");
    {
        double cutoff = 0.2;
        SOSCascade sos;
        butterworth_lowpass(4, cutoff, &sos);

        /* Evaluate magnitude at the cutoff frequency */
        int np = 257;
        double mag[257], phase[257];
        sos_freq_response(&sos, mag, phase, np);

        /* cutoff = 0.2 → index = 0.2 / 0.5 * (np-1) ≈ 102.4 */
        int idx = (int)(cutoff / 0.5 * (np - 1) + 0.5);
        double mag_db = 20.0 * log10(mag[idx]);

        /* Should be approximately -3 dB */
        if (fabs(mag_db - (-3.0)) < 0.5) {
            TEST_PASS_STMT;
        } else {
            char msg[64];
            snprintf(msg, sizeof(msg),
                     "Expected ~-3dB at cutoff, got %.2f dB", mag_db);
            TEST_FAIL_STMT(msg);
        }
    }

    /* ── Test 5: Butterworth stopband attenuation ────────────── */
    TEST_CASE_BEGIN("Butterworth LP: stopband attenuation");
    {
        SOSCascade sos;
        butterworth_lowpass(4, 0.1, &sos);

        /* At 2× cutoff (0.2), should have significant attenuation */
        int np = 257;
        double mag[257], phase[257];
        sos_freq_response(&sos, mag, phase, np);

        /* 0.3 fs → index = 0.3/0.5 * 256 ≈ 153 */
        int idx = (int)(0.3 / 0.5 * (np - 1) + 0.5);
        double mag_db = 20.0 * log10(mag[idx] + 1e-30);

        /* 4th-order = 80 dB/decade, at 3× cutoff should be > 20 dB down */
        if (mag_db < -20.0) {
            TEST_PASS_STMT;
        } else {
            char msg[64];
            snprintf(msg, sizeof(msg),
                     "Expected < -20dB at 0.3, got %.1f dB", mag_db);
            TEST_FAIL_STMT(msg);
        }
    }

    /* ── Test 6: Butterworth highpass blocks DC ──────────────── */
    TEST_CASE_BEGIN("Butterworth HP: blocks DC");
    {
        SOSCascade sos;
        butterworth_highpass(4, 0.2, &sos);

        /* DC gain should be ~0 */
        double dc = sos.gain;
        for (int s = 0; s < sos.n_sections; s++) {
            Biquad *bq = &sos.sections[s];
            double sec_dc = (bq->b0 + bq->b1 + bq->b2) /
                            (1.0 + bq->a1 + bq->a2);
            dc *= sec_dc;
        }
        if (fabs(dc) < 0.01) {
            TEST_PASS_STMT;
        } else {
            char msg[64];
            snprintf(msg, sizeof(msg), "DC gain should ≈ 0, got %.4f", dc);
            TEST_FAIL_STMT(msg);
        }
    }

    /* ── Test 7: Chebyshev I steeper than Butterworth ────────── */
    TEST_CASE_BEGIN("Chebyshev I steeper than Butterworth");
    {
        double cutoff = 0.15;
        SOSCascade bw, ch;
        butterworth_lowpass(4, cutoff, &bw);
        chebyshev1_lowpass(4, 1.0, cutoff, &ch);

        /* Compare attenuation at 2× cutoff */
        int np = 257;
        double bw_mag[257], bw_ph[257], ch_mag[257], ch_ph[257];
        sos_freq_response(&bw, bw_mag, bw_ph, np);
        sos_freq_response(&ch, ch_mag, ch_ph, np);

        int idx = (int)(0.3 / 0.5 * (np - 1) + 0.5);
        double bw_db = 20.0 * log10(bw_mag[idx] + 1e-30);
        double ch_db = 20.0 * log10(ch_mag[idx] + 1e-30);

        /* Chebyshev should have more attenuation (more negative dB) */
        if (ch_db < bw_db) {
            TEST_PASS_STMT;
        } else {
            char msg[80];
            snprintf(msg, sizeof(msg),
                     "Cheby (%.1f dB) should be < Butter (%.1f dB) at 2x cutoff",
                     ch_db, bw_db);
            TEST_FAIL_STMT(msg);
        }
    }

    /* ── Test 8: SOS impulse matches iir_filter ──────────────── */
    TEST_CASE_BEGIN("SOS impulse matches iir_filter");
    {
        SOSCascade sos;
        butterworth_lowpass(2, 0.25, &sos);

        /* Feed impulse through SOS */
        double h_sos[64];
        for (int i = 0; i < sos.n_sections; i++)
            biquad_df1_init(&sos.states[i]);

        for (int i = 0; i < 64; i++) {
            double x = (i == 0) ? 1.0 : 0.0;
            h_sos[i] = sos_process_sample(&sos, x);
        }

        /* Feed impulse through iir_filter with expanded coefficients */
        double b[3] = {sos.sections[0].b0, sos.sections[0].b1,
                       sos.sections[0].b2};
        double a[3] = {1.0, sos.sections[0].a1, sos.sections[0].a2};
        double impulse[64];
        memset(impulse, 0, sizeof(impulse));
        impulse[0] = 1.0;

        double h_direct[64];
        iir_filter(b, 3, a, 3, impulse, h_direct, 64);

        /* Apply gain */
        for (int i = 0; i < 64; i++) h_direct[i] *= sos.gain;

        int ok = 1;
        for (int i = 0; i < 64; i++) {
            if (fabs(h_sos[i] - h_direct[i]) > 1e-8) { ok = 0; break; }
        }
        if (ok) { TEST_PASS_STMT; }
        else    { TEST_FAIL_STMT("SOS and direct should match"); }
    }

    /* ── Test 9: Designed filters are stable ─────────────────── */
    TEST_CASE_BEGIN("All designed filters are stable");
    {
        /* A stable IIR has impulse response that decays to 0 */
        int ok = 1;

        SOSCascade filters[4];
        butterworth_lowpass(4, 0.2, &filters[0]);
        butterworth_lowpass(8, 0.3, &filters[1]);
        butterworth_highpass(4, 0.2, &filters[2]);
        chebyshev1_lowpass(4, 1.0, 0.2, &filters[3]);

        for (int f = 0; f < 4; f++) {
            SOSCascade *sos = &filters[f];
            for (int s = 0; s < sos->n_sections; s++)
                biquad_df1_init(&sos->states[s]);

            /* Feed impulse, check last sample is near zero */
            double y = 0.0;
            for (int i = 0; i < 1000; i++) {
                double x = (i == 0) ? 1.0 : 0.0;
                y = sos_process_sample(sos, x);
            }
            if (fabs(y) > 0.01) { ok = 0; break; }
        }
        if (ok) { TEST_PASS_STMT; }
        else    { TEST_FAIL_STMT("Impulse response should decay to ~0"); }
    }

    /* ── Test 10: Group delay of symmetric FIR is constant ───── */
    TEST_CASE_BEGIN("FIR group delay is constant");
    {
        /* Symmetric 5-tap FIR → expected group delay = 2.0 */
        double b[] = {0.1, 0.25, 0.3, 0.25, 0.1};
        double a[] = {1.0};

        int ok = 1;
        for (int i = 1; i < 16; i++) {
            double omega = M_PI * (double)i / 16.0;
            double gd = group_delay_at(b, 5, a, 1, omega);
            if (fabs(gd - 2.0) > 0.1) { ok = 0; break; }
        }
        if (ok) { TEST_PASS_STMT; }
        else    { TEST_FAIL_STMT("Symmetric FIR should have τ ≈ (N-1)/2"); }
    }

    printf("\n=== Test Summary ===\n");
    printf("Total: %d, Passed: %d, Failed: %d\n",
           test_count, test_passed, test_failed);
    printf("Pass Rate: %.1f%%\n",
           test_count > 0 ? (100.0 * test_passed / test_count) : 0.0);
    return (test_failed == 0) ? 0 : 1;
}

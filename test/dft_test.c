#include "dft.h"
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>

#define ABS(a) ((a) < 0 ? (-a) : (a))
#define DFT_TEST_ASSERT(a)                                           \
    if(!(a)) { printf("\nAssertion failed: \t" #a "\n"); return false; }

#define DFT_TEST(name) static bool dft_test_##name ()
#define DFT_TEST_END   return true;

#define DFT_TEST_ASSERT_EQ_APPROX_DELTA(a, b, d) \
    do { \
        if(((a) > (b) ? (a) - (b) : (b) - (a)) > (d)) { \
            printf("\nAssert approx failed: \t" #a " !~ " #b \
                    " \t%f !~ %f\n", (float)(a), (float)(b)); \
            return false; \
        } \
    } while(0);

#define DFT_TEST_ASSERT_EQ_APPROX(a, b) \
    DFT_TEST_ASSERT_EQ_APPROX_DELTA((a), (b), 0.001)

#define DFT_TEST_MAIN int main()
#define DFT_TEST_RUN(name)        \
    do {                          \
        printf(#name ":");        \
        if(dft_test_##name ()) {  \
            printf("PASS\n");     \
        }                         \
    } while(0);


DFT_TEST(run_with_2_zero) {

    dft_cmplx_t result[2];
    float samples[2] = { 0, 0 };

    dft_transform_cmplx(result, samples, 2);

    DFT_TEST_ASSERT_EQ_APPROX(result[0].r, 0);
    DFT_TEST_ASSERT_EQ_APPROX(result[0].i, 0);
    DFT_TEST_ASSERT_EQ_APPROX(result[1].r, 0);
    DFT_TEST_ASSERT_EQ_APPROX(result[1].i, 0);

    DFT_TEST_END;
}

DFT_TEST(run_with_2_cos) {

    dft_cmplx_t result[2];
    float samples[2] = { 1, -1 };

    dft_transform_cmplx(result, samples, 2);

    DFT_TEST_ASSERT_EQ_APPROX(result[0].r, 0);
    DFT_TEST_ASSERT_EQ_APPROX(result[0].i, 0);
    DFT_TEST_ASSERT_EQ_APPROX(result[1].r, 2);
    DFT_TEST_ASSERT_EQ_APPROX(result[1].i, 0);

    DFT_TEST_END;
}

static float msin(float a) {
    return sin(a * 2 * M_PI);
}

DFT_TEST(run_with_4_sin) {

    dft_cmplx_t result[4];
    float samples[4] = { msin(0), msin(0.25), msin(0.5), msin(0.75) };

    dft_transform_cmplx(result, samples, 4);

    DFT_TEST_ASSERT_EQ_APPROX(result[0].r, 0);
    DFT_TEST_ASSERT_EQ_APPROX(result[0].i, 0);
    DFT_TEST_ASSERT_EQ_APPROX(result[1].r, 0);
    DFT_TEST_ASSERT_EQ_APPROX(result[1].i, 2);

    DFT_TEST_END;
}

DFT_TEST(run_with_many) {
    int samples_size = 128;
    float samples[128];
    dft_cmplx_t result[128];
    for(int f = 0; f < samples_size / 2; f++) {
        for(int i = 0; i < samples_size; i++) {
            samples[i] = cos(i * f * 2 * M_PI / samples_size);
        }
        dft_transform_cmplx(result, samples, samples_size);
        int mf = -1;
        int mm = 0;
        for(int tf = 0; tf < samples_size / 2; tf++) {
            if(result[tf].r > mm) {
                mm = result[tf].r;
                mf = tf;
            }
        }
        DFT_TEST_ASSERT(mf == f);


        for(int i = 0; i < samples_size; i++) {
            samples[i] = sin(i * f * 2 * M_PI / samples_size);
        }
        dft_transform_cmplx(result, samples, samples_size);
        mf = -1;
        mm = 0;
        for(int tf = 0; tf < samples_size / 2; tf++) {
            if(result[tf].r > mm) {
                mm = result[tf].r;
                mf = tf;
            }
        }
    }

    DFT_TEST_END;
}

DFT_TEST(run_with_inverse) {
    int samples_size = 128;
    float samples[128];
    float samples_resynth[128];
    dft_cmplx_t result[128];
    for(int f = 0; f < samples_size / 2; f++) {
        for(int i = 0; i < samples_size; i++) {
            samples[i] = cos(i * f * 2 * M_PI / samples_size);
        }
        dft_transform_cmplx(result, samples, samples_size);
        dft_inverse_transform_cmplx(result, samples_resynth, samples_size);
        for(int i = 0; i < samples_size; i++) {
            DFT_TEST_ASSERT_EQ_APPROX(samples[i], samples_resynth[i]);
        }

        for(int i = 0; i < samples_size; i++) {
            samples[i] = sin(i * f * 2 * M_PI / samples_size);
        }
        dft_transform_cmplx(result, samples, samples_size);
        dft_inverse_transform_cmplx(result, samples_resynth, samples_size);
        for(int i = 0; i < samples_size; i++) {
            DFT_TEST_ASSERT_EQ_APPROX(samples[i], samples_resynth[i]);
        }
    }

    DFT_TEST_END;
}

DFT_DEFINE_FFT_WITH_TYPE(float);

DFT_TEST(run_with_fft) {
    int samples_size = 128;
    float samples[128];
    DFT_CMPLX_T(float) result[128];
    DFT_CMPLX_T(float) result_fft[128];
    for(int f = 0; f < samples_size / 2; f++) {
        for(int i = 0; i < samples_size; i++) {
            samples[i] = cos(i * f * 2 * M_PI / samples_size);
        }
        dft_transform_cmplx((dft_cmplx_t*)result, samples, samples_size);
        DFT_FFT(float)(result_fft, samples, samples_size);
        for(int i = 0; i < samples_size; i++) {
            DFT_TEST_ASSERT_EQ_APPROX(result[i].r, result_fft[i].r);
            DFT_TEST_ASSERT_EQ_APPROX(result[i].i, -result_fft[i].i);
        }

        for(int i = 0; i < samples_size; i++) {
            samples[i] = sin(i * f * 2 * M_PI / samples_size);
        }
        dft_transform_cmplx((dft_cmplx_t*)result, samples, samples_size);
        DFT_FFT(float)(result_fft, samples, samples_size);
        for(int i = 0; i < samples_size; i++) {
            DFT_TEST_ASSERT_EQ_APPROX(result[i].r, result_fft[i].r);
            DFT_TEST_ASSERT_EQ_APPROX(result[i].i, -result_fft[i].i);
        }
    }

    DFT_TEST_END;
}

DFT_TEST(run_with_ifft) {
    int samples_size = 128;
    float samples[128];
    float samples_resynth[128];
    DFT_CMPLX_T(float) result[128];
    for(int f = 0; f < samples_size / 2; f++) {
        for(int i = 0; i < samples_size; i++) {
            samples[i] = cos(i * f * 2 * M_PI / samples_size);
        }
        DFT_FFT(float)(result, samples, samples_size);
        DFT_IFFT(float)(result, samples_resynth, samples_size);
        for(int i = 0; i < samples_size; i++) {
            DFT_TEST_ASSERT_EQ_APPROX(samples[i], samples_resynth[i]);
        }

        for(int i = 0; i < samples_size; i++) {
            samples[i] = sin(i * f * 2 * M_PI / samples_size);
        }
        DFT_FFT(float)(result, samples, samples_size);
        DFT_IFFT(float)(result, samples_resynth, samples_size);
        for(int i = 0; i < samples_size; i++) {
            DFT_TEST_ASSERT_EQ_APPROX(samples[i], samples_resynth[i]);
        }
    }

    DFT_TEST_END;
}

DFT_TEST(run_with_ifft_random) {
    int samples_size = 128;
    srandom(samples_size);
    float samples[128];
    float samples_resynth[128];
    DFT_CMPLX_T(float) result[128];
    for(int i = 0; i < samples_size; i++) {
        samples[i] = random() % (1 << 8);
    }
    DFT_FFT(float)(result, samples, samples_size);
    DFT_IFFT(float)(result, samples_resynth, samples_size);
    for(int i = 0; i < samples_size; i++) {
        DFT_TEST_ASSERT_EQ_APPROX(samples[i], samples_resynth[i]);
    }

    DFT_TEST_END;
}

DFT_DEFINE_FFT_WITH_TYPE(int);

DFT_TEST(run_with_ifft_random_int) {

    int samples_size = 128;
    srandom(samples_size);
    int samples[128];
    int samples_resynth[128];
    DFT_CMPLX_T(int) result[128];
    for(int i = 0; i < samples_size; i++) {
        samples[i] = random() % (1 << 8);
    }
    DFT_FFT(int)(result, samples, samples_size);
    DFT_IFFT(int)(result, samples_resynth, samples_size);
    for(int i = 0; i < samples_size; i++) {
        DFT_TEST_ASSERT_EQ_APPROX_DELTA(samples[i], samples_resynth[i], 2);
    }

    DFT_TEST_END;
}

DFT_TEST(run_with_ifft_random_int_1024) {

    int samples_size = 1024;
    srandom(samples_size);
    int samples[1024];
    int samples_resynth[1024];
    DFT_CMPLX_T(int) result[1024];
    for(int i = 0; i < samples_size; i++) {
        samples[i] = random() % (1 << 8);
    }
    DFT_FFT(int)(result, samples, samples_size);
    DFT_IFFT(int)(result, samples_resynth, samples_size);
    for(int i = 0; i < samples_size; i++) {
        DFT_TEST_ASSERT_EQ_APPROX_DELTA(samples[i], samples_resynth[i], 2);
    }

    DFT_TEST_END;
}

DFT_DEFINE_FFT_WITH_TYPE(long);

DFT_TEST(run_with_ifft_random_long_1024) {

    int samples_size = 1024;
    srandom(samples_size);
    long samples[1024];
    long samples_resynth[1024];
    DFT_CMPLX_T(long) result[1024];
    for(int i = 0; i < samples_size; i++) {
        samples[i] = random() % (1 << 16);
    }
    DFT_FFT(long)(result, samples, samples_size);
    DFT_IFFT(long)(result, samples_resynth, samples_size);
    for(int i = 0; i < samples_size; i++) {
        DFT_TEST_ASSERT_EQ_APPROX_DELTA(samples[i], samples_resynth[i], 2);
    }

    DFT_TEST_END;
}

DFT_TEST(run_with_ifft_random_long_16384) {

    int samples_size = 16384;
    srandom(samples_size);
    long samples[16384];
    long samples_resynth[16384];
    DFT_CMPLX_T(long) result[16384];
    for(int i = 0; i < samples_size; i++) {
        samples[i] = random() % (1 << 16);
    }
    DFT_FFT(long)(result, samples, samples_size);
    DFT_IFFT(long)(result, samples_resynth, samples_size);
    for(int i = 0; i < samples_size; i++) {
        DFT_TEST_ASSERT_EQ_APPROX_DELTA(samples[i], samples_resynth[i], 2);
    }

    DFT_TEST_END;
}


DFT_TEST_MAIN {
    DFT_TEST_RUN(run_with_2_zero);
    DFT_TEST_RUN(run_with_2_cos);
    DFT_TEST_RUN(run_with_4_sin);
    DFT_TEST_RUN(run_with_many);
    DFT_TEST_RUN(run_with_inverse);
    DFT_TEST_RUN(run_with_fft);
    DFT_TEST_RUN(run_with_ifft);
    DFT_TEST_RUN(run_with_ifft_random);
    DFT_TEST_RUN(run_with_ifft_random_int);
    DFT_TEST_RUN(run_with_ifft_random_int_1024);
    DFT_TEST_RUN(run_with_ifft_random_long_1024);
    DFT_TEST_RUN(run_with_ifft_random_long_16384);
}


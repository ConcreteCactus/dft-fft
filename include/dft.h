#ifndef DFT_H_
#define DFT_H_

#include <stdint.h>
#include <assert.h>
#include <math.h>

typedef struct {
    float r, i; // r: real (cos), i: imaginary (sin)
} dft_cmplx_t;

void dft_transform(const float* samples, const uint32_t samples_size, 
        const float samples_length, const float min_freq,
        const uint32_t transform_max_size,
        float* transformation);
void dft_transform_cmplx(dft_cmplx_t* result, const float* samples,
        const uint32_t samples_size);
void dft_inverse_transform_cmplx(const dft_cmplx_t* transform, float* samples,
        const uint32_t transform_size);

#define DFT_CMPLX_T(TYP) dft_##TYP##_cmplx_t
#define DFT_FFT(TYP) dft_##TYP##_fft
#define DFT_IFFT(TYP) dft_##TYP##_ifft
#define DFT_DEFINE_FFT_WITH_TYPE(TYP)                                          \
typedef struct {                                                               \
    TYP r, i;                                                                  \
} dft_##TYP##_cmplx_t;                                                         \
                                                                               \
static inline dft_##TYP##_cmplx_t                                              \
dft_##TYP##_cmplx_mul(const dft_##TYP##_cmplx_t a, const dft_cmplx_t b) {      \
    return (dft_##TYP##_cmplx_t){a.r * b.r - a.i * b.i,                        \
        a.r * b.i + a.i * b.r};                                                \
}                                                                              \
                                                                               \
static inline dft_##TYP##_cmplx_t                                              \
dft_##TYP##_cmplx_add(const dft_##TYP##_cmplx_t a,                             \
        const dft_##TYP##_cmplx_t b) {                                         \
    return (dft_##TYP##_cmplx_t){a.r + b.r, a.i + b.i};                        \
}                                                                              \
                                                                               \
static void dft_##TYP##_multiply_w_fourier_mat(dft_##TYP##_cmplx_t* result,    \
        const TYP* samples, const uint32_t samples_size,                       \
        const uint32_t samples_distance) {                                     \
    assert(samples_size >= 2);                                                 \
                                                                               \
    if(samples_size > 2) {                                                     \
                                                                               \
        const int half_size = samples_size / 2;                                \
        assert(half_size * 2 == samples_size);                                 \
                                                                               \
        dft_##TYP##_multiply_w_fourier_mat(result, samples, half_size,         \
                samples_distance * 2);                                         \
        dft_##TYP##_multiply_w_fourier_mat(result + half_size,                 \
                samples + samples_distance, half_size, samples_distance * 2);  \
                                                                               \
        for(int i = 0; i < half_size; i++) {                                   \
            dft_cmplx_t k, mk;                                                 \
            k.r = cosf((-2.0 * M_PI * i) / samples_size);                      \
            k.i = sinf((-2.0 * M_PI * i) / samples_size);                      \
            mk.r = -k.r;                                                       \
            mk.i = -k.i;                                                       \
                                                                               \
            dft_##TYP##_cmplx_t upper_half = result[i];                        \
            dft_##TYP##_cmplx_t lower_half = result[i + half_size];            \
            result[i] = dft_##TYP##_cmplx_add(                                 \
                    upper_half,                                                \
                    dft_##TYP##_cmplx_mul(lower_half, k));                     \
            result[i + half_size] = dft_##TYP##_cmplx_add(                     \
                    upper_half,                                                \
                    dft_##TYP##_cmplx_mul(lower_half, mk));                    \
        }                                                                      \
                                                                               \
    } else {                                                                   \
        result[0] = (dft_##TYP##_cmplx_t){0, 0};                               \
        result[1] = (dft_##TYP##_cmplx_t){0, 0};                               \
        result[0].r = samples[0] + samples[samples_distance];                  \
        result[1].r = samples[0] - samples[samples_distance];                  \
    }                                                                          \
}                                                                              \
                                                                               \
static void dft_##TYP##_multiply_w_fourier_mat_cmplx_inverse(                  \
        dft_##TYP##_cmplx_t* result, const dft_##TYP##_cmplx_t* samples,       \
        const uint32_t samples_size, const uint32_t samples_distance) {        \
    assert(samples_size >= 2);                                                 \
                                                                               \
    if(samples_size > 2) {                                                     \
                                                                               \
        const int half_size = samples_size / 2;                                \
        assert(half_size * 2 == samples_size);                                 \
                                                                               \
        dft_##TYP##_multiply_w_fourier_mat_cmplx_inverse(result, samples,      \
                half_size, samples_distance * 2);                              \
        dft_##TYP##_multiply_w_fourier_mat_cmplx_inverse(result + half_size,   \
                samples + samples_distance, half_size, samples_distance * 2);  \
                                                                               \
                                                                               \
        for(int i = 0; i < half_size; i++) {                                   \
            dft_cmplx_t k, mk;                                                 \
            k.r = cosf((-2.0 * M_PI * i) / samples_size);                      \
            k.i = sinf((-2.0 * M_PI * i) / samples_size);                      \
            mk.r = -k.r;                                                       \
            mk.i = -k.i;                                                       \
                                                                               \
            dft_##TYP##_cmplx_t upper_half = result[i];                        \
            dft_##TYP##_cmplx_t lower_half = result[i + half_size];            \
            result[i] = dft_##TYP##_cmplx_add(                                 \
                    upper_half,                                                \
                    dft_##TYP##_cmplx_mul(lower_half, k));                     \
            result[i + half_size] = dft_##TYP##_cmplx_add(                     \
                    upper_half,                                                \
                    dft_##TYP##_cmplx_mul(lower_half, mk));                    \
                                                                               \
        }                                                                      \
                                                                               \
    } else {                                                                   \
        result[0] = (dft_##TYP##_cmplx_t){0, 0};                               \
        result[1] = (dft_##TYP##_cmplx_t){0, 0};                               \
        result[0].r = samples[0].r + samples[samples_distance].r;              \
        result[1].r = samples[0].r - samples[samples_distance].r;              \
        result[0].i = - samples[0].i - samples[samples_distance].i;            \
        result[1].i = - samples[0].i + samples[samples_distance].i;            \
    }                                                                          \
}                                                                              \
                                                                               \
                                                                               \
void dft_##TYP##_fft(dft_##TYP##_cmplx_t* result, const TYP* samples,          \
        const uint32_t samples_size) {                                         \
    dft_##TYP##_multiply_w_fourier_mat(result, samples, samples_size, 1);      \
}                                                                              \
                                                                               \
void dft_##TYP##_ifft(const dft_##TYP##_cmplx_t* transform, TYP* samples,      \
        const uint32_t samples_size) {                                         \
    dft_##TYP##_cmplx_t* samples_temp =                                        \
        malloc(samples_size * sizeof(dft_##TYP##_cmplx_t));                    \
    dft_##TYP##_multiply_w_fourier_mat_cmplx_inverse(samples_temp, transform,  \
            samples_size, 1);                                                  \
    for(int i = 0; i < samples_size; i++) {                                    \
        samples[i] = samples_temp[i].r / samples_size;                         \
    }                                                                          \
    free(samples_temp);                                                        \
}

#endif

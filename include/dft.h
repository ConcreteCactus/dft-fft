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
void dft_fft(dft_cmplx_t* result, const float* samples,
        const uint32_t samples_size);
void dft_ifft(const dft_cmplx_t* transform, float* samples,
        const uint32_t samples_size);

#endif

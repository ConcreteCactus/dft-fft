#include <stdint.h>

typedef struct {
    float r, i; // r: real (cos), i: imaginary (sin)
} dft_cmplx_t;

void dft_hello(void);
void dft_transform(const float* samples, const uint32_t samples_size, 
        const float samples_length, const float min_freq,
        const uint32_t transform_max_size,
        float* transformation);
void dft_fft(dft_cmplx_t* result, const float* samples,
        const uint32_t samples_size);

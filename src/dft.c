#include "dft.h"
#include <assert.h>
#include <math.h>
#include <stdbool.h>

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

void
dft_transform(const float* samples, const uint32_t samples_size, 
        const float samples_length, const float min_freq, 
        const uint32_t transform_max_size, 
        float* transformation) {

    assert(samples_size   > 0);
    assert(samples_length > 0);

    int measured_freq_cnt = samples_size / 2;
    float freq_step = 1.0 / samples_length;
    float sample_rate = samples_size / samples_length;

    for(int i = 0; i < measured_freq_cnt && i < transform_max_size; i++) {
        float freq = freq_step * i + min_freq;

        float sin_sum = 0, cos_sum = 0;
        for(int j = 0; j < samples_size; j++) {
            float angle = M_PI * 2 * freq * j / sample_rate;
            sin_sum += sin(angle) * samples[j];
            cos_sum += cos(angle) * samples[j];
        }

        sin_sum /= samples_size;
        cos_sum /= samples_size;

        transformation[i] = sqrt(sin_sum * sin_sum + cos_sum * cos_sum);
    }
}

void dft_transform_cmplx(dft_cmplx_t* result, const float* samples,
        const uint32_t samples_size) {

    for(int f = 0; f < samples_size; f++) {
        result[f] = (dft_cmplx_t){0, 0};
        for(int i = 0; i < samples_size; i++) {
            float angle = i * f * 2 * M_PI / samples_size;
            result[f].r += cos(angle) * samples[i];
            result[f].i += sin(angle) * samples[i];
        }
    }
}

void dft_inverse_transform_cmplx(const dft_cmplx_t* transform, float* samples,
        const uint32_t transform_size) {

    for(int i = 0; i < transform_size; i++) {
        samples[i] = 0;
        for(int f = 0; f < transform_size; f++) {
            float angle = i * f * 2 * M_PI / transform_size;
            samples[i] += cos(angle) * transform[f].r;
            samples[i] += sin(angle) * transform[f].i;
        }
        samples[i] /= transform_size;
    }
}


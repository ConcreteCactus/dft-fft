#include "dft.h"
#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

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

static inline dft_cmplx_t 
dft_cmplx_mul(const dft_cmplx_t a, const dft_cmplx_t b) {
    return (dft_cmplx_t){a.r * b.r - a.i * b.i, a.r * b.i + a.i * b.r};
}

static inline dft_cmplx_t
dft_cmplx_add(const dft_cmplx_t a, const dft_cmplx_t b) {
    return (dft_cmplx_t){a.r + b.r, a.i + b.i};
}

static void dft_multiply_w_fourier_mat(dft_cmplx_t* result, const float* samples, 
        const uint32_t samples_size, const uint32_t samples_distance) {
    // Assume that samples_size is a power of 2
    // Freq count will be equal to the samples_size
    assert(samples_size >= 2);
    
    if(samples_size > 2) {

        const int half_size = samples_size / 2;
        assert(half_size * 2 == samples_size);

        dft_multiply_w_fourier_mat(result, samples, half_size,
                samples_distance * 2);
        dft_multiply_w_fourier_mat(result + half_size,
                samples + samples_distance, half_size, samples_distance * 2);


        for(int i = 0; i < half_size; i++) {
            dft_cmplx_t k, mk;
            k.r = cosf((-2.0 * M_PI * i) / samples_size);
            k.i = sinf((-2.0 * M_PI * i) / samples_size);
            /* printf("k_%i_%i: \t%f \t+ %fi\n", samples_size, i, k.r, k.i); */
            mk.r = -k.r;
            mk.i = -k.i;

            dft_cmplx_t upper_half = result[i];
            dft_cmplx_t lower_half = result[i + half_size];
            result[i] = dft_cmplx_add(
                    upper_half,
                    dft_cmplx_mul(lower_half, k));
            result[i + half_size] = dft_cmplx_add(
                    upper_half, 
                    dft_cmplx_mul(lower_half, mk));

        }

    } else {
        // so samples_size == 2
        result[0] = (dft_cmplx_t){0, 0}; 
        result[1] = (dft_cmplx_t){0, 0}; 
        result[0].r = samples[0] + samples[samples_distance];
        result[1].r = samples[0] - samples[samples_distance];
    }
}

static void dft_multiply_w_fourier_mat_cmplx_inverse(dft_cmplx_t* result,
        const dft_cmplx_t* samples, const uint32_t samples_size, 
        const uint32_t samples_distance) {
    // Assume that samples_size is a power of 2
    // Freq count will be equal to the samples_size
    assert(samples_size >= 2);
    
    if(samples_size > 2) {

        const int half_size = samples_size / 2;
        assert(half_size * 2 == samples_size);

        dft_multiply_w_fourier_mat_cmplx_inverse(result, samples, half_size,
                samples_distance * 2);
        dft_multiply_w_fourier_mat_cmplx_inverse(result + half_size,
                samples + samples_distance, half_size, samples_distance * 2);


        for(int i = 0; i < half_size; i++) {
            dft_cmplx_t k, mk;
            k.r = cosf((-2.0 * M_PI * i) / samples_size);
            k.i = sinf((-2.0 * M_PI * i) / samples_size);
            /* printf("k_%i_%i: \t%f \t+ %fi\n", samples_size, i, k.r, k.i); */
            mk.r = -k.r;
            mk.i = -k.i;

            dft_cmplx_t upper_half = result[i];
            dft_cmplx_t lower_half = result[i + half_size];
            result[i] = dft_cmplx_add(
                    upper_half,
                    dft_cmplx_mul(lower_half, k));
            result[i + half_size] = dft_cmplx_add(
                    upper_half, 
                    dft_cmplx_mul(lower_half, mk));

        }

    } else {
        // so samples_size == 2
        result[0] = (dft_cmplx_t){0, 0}; 
        result[1] = (dft_cmplx_t){0, 0}; 
        result[0].r = samples[0].r + samples[samples_distance].r;
        result[1].r = samples[0].r - samples[samples_distance].r;
        result[0].i = - samples[0].i - samples[samples_distance].i;
        result[1].i = - samples[0].i + samples[samples_distance].i;
    }
}


void dft_fft(dft_cmplx_t* result, const float* samples,
        const uint32_t samples_size) {
    dft_multiply_w_fourier_mat(result, samples, samples_size, 1);
}

void dft_ifft(const dft_cmplx_t* transform, float* samples,
        const uint32_t samples_size) {
    dft_cmplx_t* samples_temp = malloc(samples_size * sizeof(dft_cmplx_t));
    dft_multiply_w_fourier_mat_cmplx_inverse(samples_temp, transform,
            samples_size, 1);
    for(int i = 0; i < samples_size; i++) {
        samples[i] = samples_temp[i].r / samples_size;
    }
    free(samples_temp);
}

#include "dft.h"
#include "SDL2/SDL.h"
#include "SDL2/SDL_events.h"
#include "SDL2/SDL_image.h"
#include "SDL2/SDL_audio.h"
#include "time.h"
#include <assert.h>
#include <math.h>
#include <stdatomic.h>
#include <stdbool.h>

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

typedef struct {
    uint16_t num_channels;
    uint16_t bits_per_sample;
    uint32_t sample_rate;
    uint32_t data_size;
    char* data;
} wav_t;

typedef struct {
    wav_t* wav;
    atomic_int pos;
    atomic_long set_time;
} playing_wav_t;

static void dft_render_graph(SDL_Renderer* renderer, const float* vals, 
        const uint32_t size, const float scale) {
    int width, height;
    SDL_GetWindowSize(SDL_RenderGetWindow(renderer), &width, &height);
    int vertical_center = height / 2;
    int bar_width = MAX(1, width / size);
    for(int i = 0; i < size; i++) {
        SDL_Rect rect = {
            .x = i * bar_width, 
            .y = vertical_center,
            .w = bar_width,
            .h = -vals[i] * scale
        };
        SDL_RenderFillRect(renderer, &rect);
    }
}

static void dft_render_cmplx_graph(SDL_Renderer* renderer,
        const dft_cmplx_t* vals, const uint32_t size, const float scale) {
    int width, height;
    SDL_GetWindowSize(SDL_RenderGetWindow(renderer), &width, &height);
    int vertical_center = height / 2;
    int bar_width = MAX(1, width / size / 2);
    for(int i = 0; i < size; i++) {
        SDL_Rect rect = {
            .x = 2 * i * bar_width, 
            .y = vertical_center,
            .w = bar_width,
            .h = -sqrt(vals[i].r * vals[i].r + vals[i].i * vals[i].i) * scale
        };
        SDL_RenderFillRect(renderer, &rect);
        /* rect = (SDL_Rect){ */
        /*     .x = (2 * i + 1) * bar_width,  */
        /*     .y = vertical_center, */
        /*     .w = bar_width, */
        /*     .h = -vals[i].i * scale */
        /* }; */
        /* SDL_RenderFillRect(renderer, &rect); */
    }
}

static void dft_render_cmplx_graph_lg(SDL_Renderer* renderer,
        const dft_cmplx_t* vals, const uint32_t size, const float scale) {
    int width, height;
    SDL_GetWindowSize(SDL_RenderGetWindow(renderer), &width, &height);
    int vertical_center = height / 2;
    int bar_width = MAX(1, width / size);
    int cnt_bar = bar_width == 1 ? ceil(size / (float)width) : 1;
    float max = 0;
    int maxi = 0;
    for(int i = 0; i < size / cnt_bar; i++) {
        float avg = 0;
        for(int j = 0; j < cnt_bar; j++) {
            float val = sqrt(vals[i].r * vals[i].r + vals[i].i * vals[i].i) * scale;
            avg += val;
            if(val > max) {
                max = val;
                maxi = i;
            }
        }
        avg /= cnt_bar;
        SDL_Rect rect = {
            .x = i * bar_width, 
            .y = vertical_center,
            .w = bar_width,
            .h = -avg
        };
        SDL_RenderFillRect(renderer, &rect);
    }
    printf("\rmaxI: %i", maxi);
    fflush(stdout);
}

static wav_t* dft_load_wav(const char* path) {
    FILE* file = fopen(path, "rb");
    char ids[4];
    wav_t* wav = NULL;
    int cnt = fread(ids, sizeof(char), 4, file);
    if(cnt != 4 || strncmp("RIFF", ids, 4) != 0) {
        fprintf(stderr, "File format unrecognized. (RIFF)\n");
        goto done;
    }
    fseek(file, 4, SEEK_CUR);
    cnt = fread(ids, sizeof(char), 4, file);
    if(cnt != 4 || strncmp("WAVE", ids, 4) != 0) {
        fprintf(stderr, "File format unrecognized. (WAVE)\n");
        goto done;
    }

    cnt = fread(ids, sizeof(char), 4, file);
    if(cnt != 4 || strncmp("fmt ", ids, 4) != 0) {
        fprintf(stderr, "File format unrecognized. (fmt )\n");
        goto done;
    }

    uint32_t header_size;
    cnt = fread(&header_size, sizeof(header_size), 1, file);
    if(cnt != 1) {
        fprintf(stderr, "File ends abruptly at header size. (%i)\n", cnt);
        goto done;
    }

    uint16_t format;
    cnt = fread(&format, sizeof(format), 1, file);
    if(cnt != 1 || format != 1) {
        fprintf(stderr, "File format unrecognized. (format)\n");
        goto done;
    }

    wav = calloc(1, sizeof(wav_t));
    cnt = fread(&wav->num_channels, sizeof(wav->num_channels), 1, file);
    if(cnt != 1) {
        fprintf(stderr, "File ends abruptly at num channels.\n");
        goto done;
    }
    cnt = fread(&wav->sample_rate, sizeof(wav->sample_rate), 1, file);
    if(cnt != 1) {
        fprintf(stderr, "File ends abruptly at sample rate.\n");
        goto done;
    }
    fseek(file, 6, SEEK_CUR);
    cnt = fread(&wav->bits_per_sample, sizeof(wav->bits_per_sample), 1, file);
    if(cnt != 1) {
        fprintf(stderr, "File ends abruptly at bits per sample.\n");
        goto done;
    }

    uint32_t seek_cnt = header_size - 16;
    do {
        fseek(file, seek_cnt, SEEK_CUR);
        cnt = fread(ids, sizeof(char), 4, file);
        cnt += fread(&seek_cnt, sizeof(seek_cnt), 1, file);
        if(cnt != 5) {
            fprintf(stderr, "File ends abruptly while searching for" 
                    " data section.\n");
            goto done;
        }
    } while(strncmp("data", ids, 4) != 0);

    wav->data_size = seek_cnt;
    wav->data = malloc(wav->data_size * sizeof(char));
    if(!wav->data) {
        fprintf(stderr, "Couldn't allocate memory for data.\n");
        goto done;
    }
    cnt = fread(wav->data, sizeof(char) * wav->bits_per_sample / 8, 
            wav->data_size / (wav->bits_per_sample / 8), file);
    if(cnt != wav->data_size / (wav->bits_per_sample / 8)) {
        fprintf(stderr, "File ends abruptly at data.\n");
        goto done;
    }

done:
    fclose(file);
    return wav;
}

static void dft_destroy_wav(wav_t* wav) {
    if(wav->data) {
        free(wav->data);
    } 
    free(wav);
}

static int dft_extract_floats_wav(wav_t* wav, int channel_id, int start_at,
        int length, float* data) {
    int bytes_per_sample = wav->bits_per_sample >> 3;
    int i = 0;
    for(;
        i < length && 
            (start_at + i) * bytes_per_sample * wav->num_channels 
            < wav->data_size;
        i++) {
        if(bytes_per_sample == 1) {
            uint8_t* wav_data = (uint8_t*)wav->data;
            data[i] = ((float)wav_data[
                        (start_at + i) * wav->num_channels + channel_id] 
                    / 128.0) - 1;
        } else if(bytes_per_sample == 2) {
            int16_t* wav_data = (int16_t*)wav->data;
            data[i] = (float)wav_data[
                        (start_at + i) * wav->num_channels + channel_id] 
                    / (1 << (sizeof(int16_t) * 8 - 1));
        } else {
            break;
        }
    }
    return i;
}

static void dft_add_wave(const float freq, const float phase, 
        const float ampl, const uint32_t size, const float length, 
        float* wave) {
    const float sample_rate = size / length;
    for(int i = 0; i < size; i++) {
        const float angle = M_PI * 2.0 * (i / sample_rate) * freq;
        const float phase_angle = M_PI * 2.0 * phase * freq;
        wave[i] += cosf(phase_angle + angle) * ampl;
    }
}

static void dft_sdl_audio_callback_440hz(void* data, uint8_t* stream, int len) {
    playing_wav_t* pwav = data;
    const int len_float = len / sizeof(float);

    memset(stream, 0, len);
    dft_add_wave(
            440, 
            pwav->pos / 44100.0,
            0.4,
            len_float, 
            len_float / 44100.0, (float*)stream);

    pwav->pos = (pwav->pos + len_float) % 44100;

}

static void dft_sdl_audio_callback_wav_one_channel(void* data, uint8_t* stream,
        int len) {
    playing_wav_t* pwav = data;
    int wav_channels = pwav->wav->num_channels;
    int bytes_per_sample = pwav->wav->bits_per_sample >> 3;
    int required_len = len / bytes_per_sample;
    int data_frame_count = 
        pwav->wav->data_size / wav_channels / bytes_per_sample;

    uint16_t* frame_stream = (uint16_t*)stream;
    uint16_t* wav_data = (uint16_t*)pwav->wav->data;
    for(int i = 0;
        i < required_len || pwav->pos + i > data_frame_count;
        i += wav_channels * bytes_per_sample) {
        frame_stream[i] = wav_data[(pwav->pos + i) * wav_channels];
    }
    pwav->pos += required_len;
}

static void dft_sdl_audio_callback_wav(void* data, uint8_t* stream,
        int len) {
    playing_wav_t* pwav = data;
    int32_t pos = atomic_load_explicit(&pwav->pos, memory_order_acquire);
    float* sample_stream = (float*)stream;
    int16_t* wav_data = (int16_t*)pwav->wav->data;
    uint32_t sample_count = pwav->wav->data_size / sizeof(int16_t);
    for(int i = 0; i < len / sizeof(float); i++) {
        sample_stream[i] = 
            (float)wav_data[(pos + i) % sample_count] 
            / (1 << (sizeof(int16_t) * 8));
    }
    if(pos + len / sizeof(float) >= sample_count) {
        pos -= sample_count;
    }
    atomic_store_explicit(
            &pwav->pos, pos + len / sizeof(float), memory_order_release);
    atomic_store_explicit(&pwav->set_time, SDL_GetTicks64(), memory_order_release);
}

void dft_hello(void) {
    SDL_Init(SDL_INIT_VIDEO | SDL_INIT_AUDIO);
    SDL_Window* window = SDL_CreateWindow("Window", SDL_WINDOWPOS_CENTERED,
            SDL_WINDOWPOS_CENTERED, 1920, 1080, SDL_WINDOW_SHOWN | SDL_WINDOW_RESIZABLE);
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, 0);

    SDL_Surface* surface = IMG_Load("images/table-food.jpg");
    SDL_Texture* texture = SDL_CreateTextureFromSurface(renderer, surface);

    wav_t* wav = dft_load_wav("sounds/kalimba.wav");
    playing_wav_t pwav = {
        .wav = wav
    };
    atomic_init(&pwav.pos, 0);
    atomic_init(&pwav.set_time, SDL_GetTicks64());

    SDL_AudioSpec spec = {0};
    spec.freq = wav->sample_rate;
    spec.format = AUDIO_F32LSB;
    spec.channels = wav->num_channels;
    spec.samples = 4096;
    spec.callback = dft_sdl_audio_callback_wav;
    spec.userdata = &pwav;
    SDL_AudioDeviceID device_id = SDL_OpenAudioDevice(NULL, 0, &spec, NULL, 0);
    SDL_PauseAudioDevice(device_id, 0);

    int wave_len = 1 << 12;
    float* wave = malloc(wave_len * sizeof(float));

    dft_cmplx_t* transform = calloc(1, wave_len * sizeof(dft_cmplx_t));

    bool render_transform = true;

    SDL_Event e;
    bool should_exit = false;
    while(!should_exit) {
        while(SDL_PollEvent(&e)) {
            switch(e.type) {
                case SDL_QUIT:
                    should_exit = true;
                    break;
                case SDL_KEYUP:
                    if(e.key.keysym.sym == SDLK_p) {
                        render_transform = !render_transform;
                    }
                    break;
            }
        }

        int pos = atomic_load_explicit(&pwav.pos, memory_order_acquire) / wav->num_channels;
        int set_time = atomic_load_explicit(&pwav.set_time, memory_order_acquire);
        float time_elapsed = (SDL_GetTicks64() - set_time) / 1000.0;
        int corrected_pos = time_elapsed * wav->sample_rate + pos;
        int extracted = dft_extract_floats_wav(
                wav, 0, corrected_pos, wave_len, wave);
        fflush(stdout);
        /* memset(wave, 0, sizeof(float) * wave_len); */
        /* dft_add_wave(1, 0, 1, wave_len, 1, wave); */
        /* dft_add_wave(3, 0, 1, wave_len, 1, wave); */
        /* dft_add_wave(9, 0, 1, wave_len, 1, wave); */
        dft_fft(transform, wave, wave_len);
        

        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);
        SDL_SetRenderDrawColor(renderer, 0, 255, 255, 255);
        if(render_transform) {
            dft_render_cmplx_graph_lg(
                    renderer, transform, wave_len / 8, 5000.0 / wave_len);
        } else {
            dft_render_graph(renderer, wave, wave_len, 100);
        }
        SDL_RenderPresent(renderer);
    }

    free(wave);
    free(transform);

    SDL_DestroyTexture(texture);
    SDL_FreeSurface(surface);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    dft_destroy_wav(wav);
}

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

static inline dft_cmplx_t 
dft_cmplx_mul(const dft_cmplx_t a, const dft_cmplx_t b) {
    return (dft_cmplx_t){a.r * b.r - a.i * b.i, a.r * b.i + a.i * b.r};
}

static inline dft_cmplx_t
dft_cmplx_add(const dft_cmplx_t a, const dft_cmplx_t b) {
    return (dft_cmplx_t){a.r + b.r, a.i + b.i};
}

dft_cmplx_t* g_result;

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

void dft_fft(dft_cmplx_t* result, const float* samples,
        const uint32_t samples_size) {
    g_result = result;
    dft_multiply_w_fourier_mat(result, samples, samples_size, 1);
}

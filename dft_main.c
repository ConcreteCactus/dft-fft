#include <SDL2/SDL.h>
#include <SDL2/SDL_video.h>
#include <SDL2/SDL_render.h>
#include <stdbool.h>

bool fill_texture(SDL_Texture* texture, uint8_t r, uint8_t g, uint8_t b) {
    SDL_Surface* surface;
    if(SDL_LockTextureToSurface(texture, NULL, &surface)) {
        return false;
    }

    SDL_PixelFormat* format = surface->format;

    int w = surface->w;
    int h = surface->h;

    if(format->BytesPerPixel != 2) {
        fprintf(stderr, "Not supported pixel format");
        abort();
    }

    uint16_t* pixels = surface->pixels;

    for(int i = 0; i < h; i++) {
        for(int j = 0; j < w; j++) {
            pixels[i * w + j] = SDL_MapRGB(format, r, g, b);
        }
    }

    SDL_UnlockTexture(texture);
    return true;
}

bool transform_surface(SDL_Surface* in_surface, SDL_Surface* out_surface) {
}

int main() {
    SDL_Init(SDL_INIT_VIDEO);

    SDL_Window* window = SDL_CreateWindow("IMAGE!!!! ", SDL_WINDOWPOS_UNDEFINED, 
            SDL_WINDOWPOS_UNDEFINED, 800, 600, SDL_WINDOW_RESIZABLE | 
            SDL_WINDOW_SHOWN);
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1,
            SDL_RENDERER_PRESENTVSYNC);

    SDL_Texture* texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGB565,
            SDL_TEXTUREACCESS_STREAMING, 800, 600);
    if(!texture) {
        printf("Couldn't allocate texture: %s\n", SDL_GetError());
        exit(1);
    }

    if(!fill_texture(texture, 128, 0, 255)) {
        printf("Couldn't lock surface for filling.\n");
    }

    bool shouldExit = false;
    SDL_Event event;
    while(!shouldExit) {
        while(SDL_PollEvent(&event)) {
            switch(event.type) {
                case SDL_QUIT:
                    shouldExit = true;
                    break;
                default:
                    break;
            }
        }

        SDL_RenderClear(renderer);
        SDL_RenderCopy(renderer, texture, NULL, NULL);
        SDL_RenderPresent(renderer);
    }

    SDL_DestroyTexture(texture);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
}

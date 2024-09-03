cmake_minimum_required(VERSION 3.29)
project(discrete_fourier_transform)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(PkgConfig REQUIRED)
pkg_check_modules(sdl2 REQUIRED sdl2)
pkg_check_modules(sdl2_image REQUIRED SDL2_image)

add_library(dft_lib STATIC src/dft.c)
target_link_libraries(dft_lib 
    m
    ${sdl2_LIBRARIES}
    ${sdl2_image_LIBRARIES})
target_include_directories(dft_lib 
    PUBLIC 
    include 
    ${sdl2_INCLUDE_DIRS}
    ${sdl2_image_INCLUDE_DIRS})
target_compile_options(dft_lib 
    PUBLIC 
    ${sdl2_CFLAGS_OTHER}
    ${sdl2_image_CFLAGS_OTHER})

add_executable(dft_main dft_main.c)
target_link_libraries(dft_main
    PRIVATE dft_lib)


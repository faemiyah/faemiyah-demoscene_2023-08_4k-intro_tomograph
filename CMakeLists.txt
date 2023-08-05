cmake_minimum_required(VERSION 3.1)
set(CMAKE_CONFIGURATION_TYPES "Debug;Release")
project(tomograph)

enable_language(C)
enable_language(CXX)

include("defaults.cmake")

find_boost("filesystem" "program_options" "system" "wave")
find_opengl()
find_glew()
find_png()
find_sdl2()
find_threads()

output_flags("DEBUG" on)

add_definitions(-DUSE_LD)

include_directories("${PROJECT_SOURCE_DIR}/src")

add_executable(tomograph
    "src/bsd_rand.c"
    "src/bsd_rand.h"
    "src/dnload.h"
    "src/glsl_pipeline.cpp"
    "src/glsl_pipeline.hpp"
    "src/glsl_program.cpp"
    "src/glsl_program.hpp"
    "src/glsl_shader.cpp"
    "src/glsl_shader.hpp"
    "src/glsl_shader_program.cpp"
    "src/glsl_shader_program.hpp"
    "src/glsl_shader_source.cpp"
    "src/glsl_shader_source.hpp"
    "src/glsl_wave.cpp"
    "src/glsl_wave.hpp"
    "src/gnu_rand.c"
    "src/gnu_rand.h"
    "src/header.glsl"
    "src/image_png.cpp"
    "src/image_png.hpp"
    "src/intro.cpp"
    "src/post.frag.glsl"
    "src/post.frag.glsl.hpp"
    "src/quad.frag.glsl"
    "src/quad.frag.glsl.hpp"
    "src/synth.h")
if(MSVC)
    target_link_libraries("tomograph" "Imm32")
    target_link_libraries("tomograph" "Setupapi")
    target_link_libraries("tomograph" "version")
    target_link_libraries("tomograph" "WINMM")
    target_link_libraries("tomograph" "${OPENGL_gl_LIBRARY}")
    target_link_libraries("tomograph" debug "${GLEW_LIBRARY_DEBUG}" optimized "${GLEW_LIBRARY}")
    target_link_libraries("tomograph" debug "${PNG_LIBRARY_DEBUG}" optimized "${PNG_LIBRARY}")
    target_link_libraries("tomograph" debug "${SDL2_LIBRARY_DEBUG}" optimized "${SDL2_LIBRARY}")
    target_link_libraries("tomograph" debug "${ZLIB_LIBRARY_DEBUG}" optimized "${ZLIB_LIBRARY}")
else()
    target_link_libraries("tomograph" "${GLEW_LIBRARY}")
    target_link_libraries("tomograph" "${BOOST_FILESYSTEM_LIBRARY}")
    target_link_libraries("tomograph" "${BOOST_PROGRAM_OPTIONS_LIBRARY}")
    target_link_libraries("tomograph" "${BOOST_SYSTEM_LIBRARY}")
    target_link_libraries("tomograph" "${BOOST_THREAD_LIBRARY}")
    target_link_libraries("tomograph" "${BOOST_WAVE_LIBRARY}")
    target_link_libraries("tomograph" "${PNG_LIBRARY}")
    target_link_libraries("tomograph" "${SDL2_LIBRARY}")
endif()
target_link_libraries("tomograph" "${THREADS_LIBRARY}")

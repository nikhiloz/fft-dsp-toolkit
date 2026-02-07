cmake_minimum_required(VERSION 3.15)
project(dsp-tutorial-suite 
    VERSION 1.0.0 
    DESCRIPTION "Hands-on C tutorial for Digital Signal Processing"
    LANGUAGES C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -Werror -fPIC")
set(CMAKE_C_FLAGS_DEBUG "-g -O0 -DDEBUG")
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG")

# Default to Release build
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Enable testing
enable_testing()

# Include directories
include_directories(include)

# Core library sources
set(DSP_SOURCES
    src/fft.c
    src/filter.c
    src/dsp_utils.c
    src/signal_gen.c
    src/convolution.c
)

# Create shared library
add_library(fft_dsp_shared SHARED ${DSP_SOURCES})
set_target_properties(fft_dsp_shared PROPERTIES
    OUTPUT_NAME fft_dsp
    VERSION 1.0.0
    SOVERSION 1
)

# Create static library
add_library(fft_dsp_static STATIC ${DSP_SOURCES})
set_target_properties(fft_dsp_static PROPERTIES
    OUTPUT_NAME fft_dsp
)

# Chapter demos
add_executable(ch01s chapters/01-signals-and-sequences.c)
target_link_libraries(ch01s fft_dsp_static m)

add_executable(ch02s chapters/02-sampling-and-aliasing.c)
target_link_libraries(ch02s fft_dsp_static m)

add_executable(ch04s chapters/04-lti-systems.c)
target_link_libraries(ch04s fft_dsp_static m)

add_executable(ch05s chapters/05-z-transform.c)
target_link_libraries(ch05s fft_dsp_static m)

add_executable(ch07s chapters/07-dft-theory.c)
target_link_libraries(ch07s fft_dsp_static m)

add_executable(ch01 chapters/01-complex-numbers.c)
target_link_libraries(ch01 fft_dsp_static m)

add_executable(ch02 chapters/02-fft-fundamentals.c)
target_link_libraries(ch02 fft_dsp_static m)

add_executable(ch03 chapters/03-window-functions.c)
target_link_libraries(ch03 fft_dsp_static m)

add_executable(ch04 chapters/04-digital-filters.c)
target_link_libraries(ch04 fft_dsp_static m)

add_executable(ch05 chapters/05-spectral-analysis.c)
target_link_libraries(ch05 fft_dsp_static m)

add_executable(ch08 chapters/08-putting-it-together.c)
target_link_libraries(ch08 fft_dsp_static m)

# Unit Tests (minimal framework, no external dependencies)
add_executable(test_fft tests/test_fft.c)
target_include_directories(test_fft PRIVATE tests)
target_link_libraries(test_fft fft_dsp_static m)
add_test(NAME fft_tests COMMAND test_fft)

add_executable(test_filter tests/test_filter.c)
target_include_directories(test_filter PRIVATE tests)
target_link_libraries(test_filter fft_dsp_static m)
add_test(NAME filter_tests COMMAND test_filter)

# Build output directories
set_target_properties(
    ch01s ch02s ch04s ch05s ch07s ch01 ch02 ch03 ch04 ch05 ch08 test_fft test_filter
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

set_target_properties(
    fft_dsp_shared fft_dsp_static
    PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)

# Installation
install(TARGETS fft_dsp_shared fft_dsp_static
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
install(DIRECTORY include/ DESTINATION include)

# Formatting target (clang-format)
add_custom_target(format
    COMMAND clang-format -i src/*.c include/*.h chapters/*.c tests/test_*.c
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

# Code coverage target
if(ENABLE_COVERAGE)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --coverage")
    add_custom_target(coverage
        COMMAND ctest
        COMMAND lcov --directory . --capture --output-file coverage.info
        COMMAND lcov --remove coverage.info '/usr/*' '*/tests/*' --output-file coverage.info
        COMMAND genhtml coverage.info --output-directory coverage_html
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )
endif()

# Print configuration
message(STATUS "=== DSP Tutorial Suite Configuration ===")
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C Compiler: ${CMAKE_C_COMPILER}")
message(STATUS "C Flags: ${CMAKE_C_FLAGS}")
message(STATUS "C++ Flags: ${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}}")
message(STATUS "====================================")

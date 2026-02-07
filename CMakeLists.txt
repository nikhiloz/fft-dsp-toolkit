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
    src/iir.c
    src/dsp_utils.c
    src/signal_gen.c
    src/convolution.c
)

# Create shared library
add_library(dsp_core_shared SHARED ${DSP_SOURCES})
set_target_properties(dsp_core_shared PROPERTIES
    OUTPUT_NAME dsp_core
    VERSION 1.0.0
    SOVERSION 1
)

# Create static library
add_library(dsp_core_static STATIC ${DSP_SOURCES})
set_target_properties(dsp_core_static PROPERTIES
    OUTPUT_NAME dsp_core
)

# Chapter demos
add_executable(ch01 chapters/01-signals-and-sequences.c)
target_link_libraries(ch01 dsp_core_static m)

add_executable(ch02 chapters/02-sampling-and-aliasing.c)
target_link_libraries(ch02 dsp_core_static m)

add_executable(ch03 chapters/03-complex-numbers.c)
target_link_libraries(ch03 dsp_core_static m)

add_executable(ch04 chapters/04-lti-systems.c)
target_link_libraries(ch04 dsp_core_static m)

add_executable(ch05 chapters/05-z-transform.c)
target_link_libraries(ch05 dsp_core_static m)

add_executable(ch06 chapters/06-frequency-response.c)
target_link_libraries(ch06 dsp_core_static m)

add_executable(ch07 chapters/07-dft-theory.c)
target_link_libraries(ch07 dsp_core_static m)

add_executable(ch08 chapters/08-fft-fundamentals.c)
target_link_libraries(ch08 dsp_core_static m)

add_executable(ch09 chapters/09-window-functions.c)
target_link_libraries(ch09 dsp_core_static m)

add_executable(ch10 chapters/10-digital-filters.c)
target_link_libraries(ch10 dsp_core_static m)

add_executable(ch11 chapters/11-iir-filter-design.c)
target_link_libraries(ch11 dsp_core_static m)

add_executable(ch12 chapters/12-filter-structures.c)
target_link_libraries(ch12 dsp_core_static m)

add_executable(ch13 chapters/13-spectral-analysis.c)
target_link_libraries(ch13 dsp_core_static m)

add_executable(ch30 chapters/30-putting-it-together.c)
target_link_libraries(ch30 dsp_core_static m)

# Unit Tests (minimal framework, no external dependencies)
add_executable(test_fft tests/test_fft.c)
target_include_directories(test_fft PRIVATE tests)
target_link_libraries(test_fft dsp_core_static m)
add_test(NAME fft_tests COMMAND test_fft)

add_executable(test_filter tests/test_filter.c)
target_include_directories(test_filter PRIVATE tests)
target_link_libraries(test_filter dsp_core_static m)
add_test(NAME filter_tests COMMAND test_filter)

add_executable(test_iir tests/test_iir.c)
target_include_directories(test_iir PRIVATE tests)
target_link_libraries(test_iir dsp_core_static m)
add_test(NAME iir_tests COMMAND test_iir)

# Build output directories
set_target_properties(
    ch01 ch02 ch03 ch04 ch05 ch06 ch07 ch08 ch09 ch10 ch11 ch12 ch13 ch30 test_fft test_filter test_iir
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

set_target_properties(
    dsp_core_shared dsp_core_static
    PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)

# Installation
install(TARGETS dsp_core_shared dsp_core_static
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

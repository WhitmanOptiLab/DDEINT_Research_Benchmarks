cmake_minimum_required(VERSION 3.10)
project(BenchmarkTests)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add the Google Benchmark library
add_subdirectory(../utils/benchmark build_benchmark)

# Include directories
include_directories(../DDEINT ../utils/benchmark/include ../test/)

# Glucose Insulin Model
set(GI_SOURCES
    Glucose_Insulin_Model/gi_benchmark.cpp
    Glucose_Insulin_Model/gi_functions.cpp
)

add_executable(benchmark_gi ${GI_SOURCES})
target_link_libraries(benchmark_gi benchmark::benchmark benchmark::benchmark_main pthread)

# Cardiovascular Model
set(CV_SOURCES
    Cardiovascular_Model/cv_benchmark.cpp
    Cardiovascular_Model/cv_functions.cpp
)

add_executable(benchmark_cv ${CV_SOURCES})
target_link_libraries(benchmark_cv benchmark::benchmark benchmark::benchmark_main pthread)

# Breast Cancer Model
set(BC_SOURCES
    Breast_Cancer_Model/bc_benchmark.cpp
    Breast_Cancer_Model/bc_functions.cpp
)

add_executable(benchmark_bc ${BC_SOURCES})
target_link_libraries(benchmark_bc benchmark::benchmark benchmark::benchmark_main pthread)

# Set the build type to Release
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

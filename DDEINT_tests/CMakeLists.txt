cmake_minimum_required(VERSION 3.10)

# Use real GCC (Homebrew) instead of AppleClang — needed for std::execution support
set(CMAKE_C_COMPILER   /opt/homebrew/bin/gcc-16)
set(CMAKE_CXX_COMPILER /opt/homebrew/bin/g++-16)

project(BenchmarkTests)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add the Google Benchmark library
add_subdirectory(../utils/benchmark build_benchmark)

# Include directories
include_directories(../DDEINT ../utils/benchmark/include ../DDEINT_tests/Glucose_Insulin_Model ../DDEINT_tests/Cardiovascular_Model ../DDEINT_tests/Breast_Cancer_Model ../DDEINT_tests/Cobweb_Model)

# TBB (needed for std::execution::par/seq under libstdc++)
find_program(BREW_EXECUTABLE brew)
execute_process(COMMAND ${BREW_EXECUTABLE} --prefix tbb
                OUTPUT_VARIABLE TBB_PREFIX
                OUTPUT_STRIP_TRAILING_WHITESPACE)
include_directories(${TBB_PREFIX}/include)
link_directories(${TBB_PREFIX}/lib)

# Glucose Insulin Model
set(GI_SOURCES
   Glucose_Insulin_Model/gi_benchmark.cpp
    Glucose_Insulin_Model/gi_functions.cpp
)

add_executable(benchmark_gi ${GI_SOURCES})
target_link_libraries(benchmark_gi benchmark::benchmark benchmark::benchmark_main pthread tbb)

# Cardiovascular Model
#set(CV_SOURCES
#    Cardiovascular_Model/cv_benchmark.cpp
#    Cardiovascular_Model/cv_functions.cpp
#)

#add_executable(benchmark_cv ${CV_SOURCES})
#target_link_libraries(benchmark_cv benchmark::benchmark benchmark::benchmark_main pthread tbb)

# Breast Cancer Model
set(BC_SOURCES
    Breast_Cancer_Model/bc_benchmark.cpp
    Breast_Cancer_Model/bc_functions.cpp
)

add_executable(benchmark_bc ${BC_SOURCES})
target_link_libraries(benchmark_bc benchmark::benchmark benchmark::benchmark_main pthread tbb)


# Cobweb model
#set(CW_SOURCES
#    Cobweb_Model/cw_benchmark.cpp
#    Cobweb_Model/cw_functions.cpp
#)

#add_executable(benchmark_cw ${CW_SOURCES})
#target_link_libraries(benchmark_cw benchmark::benchmark benchmark::benchmark_main pthread tbb)


# Set the build type to Release
if(NOT CMAKE_BUILD_TYPE)
set(CMAKE_BUILD_TYPE Release)
endif()
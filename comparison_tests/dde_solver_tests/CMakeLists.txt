cmake_minimum_required(VERSION 3.10)
project(DDE_Solver_Comparison Fortran CXX)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add the Google Benchmark library
add_subdirectory(../../utils/benchmark build_benchmark)

# dde_solver Fortran library
add_library(dde_solver_lib STATIC
    ../../comparison_libraries/dde_solver/dde_solver_m.f90
)

# Breast Cancer Model
add_library(bc_fortran STATIC
    Breast_Cancer_Model/define_ddes.f90
    Breast_Cancer_Model/bc_benchmark.f90
)
target_link_libraries(bc_fortran dde_solver_lib)

set(BC_SOURCES
    Breast_Cancer_Model/bc_benchmark_wrapper.cpp
)
add_executable(benchmark_bc ${BC_SOURCES})
target_link_libraries(benchmark_bc bc_fortran dde_solver_lib benchmark::benchmark benchmark::benchmark_main gfortran pthread)








# Set the build type to Release
if(NOT CMAKE_BUILD_TYPE)
set(CMAKE_BUILD_TYPE Release)
endif()
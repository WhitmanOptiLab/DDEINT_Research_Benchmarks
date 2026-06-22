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
    Breast_Cancer_Model/bc_define_ddes.f90
    Breast_Cancer_Model/bc_benchmark.f90
)
target_link_libraries(bc_fortran dde_solver_lib)

set(BC_SOURCES
    Breast_Cancer_Model/bc_benchmark_wrapper.cpp
)
add_executable(benchmark_bc ${BC_SOURCES})
target_link_libraries(benchmark_bc bc_fortran dde_solver_lib benchmark::benchmark benchmark::benchmark_main gfortran pthread)

# Cardiovascular Model
add_library(cv_fortran STATIC
    Cardiovascular_Model/cv_define_ddes.f90
    Cardiovascular_Model/cv_benchmark.f90
)
target_link_libraries(cv_fortran dde_solver_lib)

set(CV_SOURCES
    Cardiovascular_Model/cv_benchmark_wrapper.cpp
)
add_executable(benchmark_cv ${CV_SOURCES})
target_link_libraries(benchmark_cv cv_fortran dde_solver_lib benchmark::benchmark benchmark::benchmark_main gfortran pthread)

# Cardiovascular Model
add_library(cw_fortran STATIC
    Cobweb_Model/cw_define_ddes.f90
    Cobweb_Model/cw_benchmark.f90
)
target_link_libraries(cw_fortran dde_solver_lib)

set(CW_SOURCES
    Cobweb_Model/cw_benchmark_wrapper.cpp
)
add_executable(benchmark_cw ${CW_SOURCES})
target_link_libraries(benchmark_cw cw_fortran dde_solver_lib benchmark::benchmark benchmark::benchmark_main gfortran pthread)

# Glucose Insulin Model
add_library(gi_fortran STATIC
    Glucose_Insulin_Model/gi_define_ddes.f90
    Glucose_Insulin_Model/gi_benchmark.f90
)
target_link_libraries(gi_fortran dde_solver_lib)

set(GI_SOURCES
    Glucose_Insulin_Model/gi_benchmark_wrapper.cpp
)
add_executable(benchmark_gi ${GI_SOURCES})
target_link_libraries(benchmark_gi gi_fortran dde_solver_lib benchmark::benchmark benchmark::benchmark_main gfortran pthread)





# Set the build type to Release
if(NOT CMAKE_BUILD_TYPE)
set(CMAKE_BUILD_TYPE Release)
endif()
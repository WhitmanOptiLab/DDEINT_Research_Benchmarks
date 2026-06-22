#!/bin/bash
# This script is used to run the benchmark tests 
BUILD_DIR="bench_build"
DATA_DIR="data/bench_data"
mkdir -p $BUILD_DIR $DATA_DIR
git submodule update --init --recursive
# compile and run benchmark test for each model
# Build the project using CMake
cd $BUILD_DIR
cmake .. -DCMAKE_BUILD_TYPE=Release -DBENCHMARK_ENABLE_TESTING=OFF 
make
cd ..

# Run the benchmark tests
# save data to benchmark_results.txt in the data/bench_data directory

# Breast Cancer Model 
./$BUILD_DIR/benchmark_bc --benchmark_format=console --benchmark_out=$DATA_DIR/bc_benchmark_results.txt
if [ $? -eq 0 ]; then
    echo "Breast Cancer Benchmark tests passed"
else
    echo "Breast Cancer Benchmark tests failed"
    exit 1
fi

# Cardiovascular Model
./$BUILD_DIR/benchmark_cv --benchmark_format=console --benchmark_out=$DATA_DIR/cv_benchmark_results.txt
if [ $? -eq 0 ]; then
    echo "Cardiovascular Benchmark tests passed"
else
    echo "Cardiovascular Benchmark tests failed"
    exit 1
fi

# Cobweb Model
./$BUILD_DIR/benchmark_cw --benchmark_format=console --benchmark_out=$DATA_DIR/cw_benchmark_results.txt
if [ $? -eq 0 ]; then
    echo "Cobweb Benchmark tests passed"
else
    echo "Cobweb Benchmark tests failed"
    exit 1
fi

# Glucose Insulin Model
./$BUILD_DIR/benchmark_gi --benchmark_format=console --benchmark_out=$DATA_DIR/gi_benchmark_results.txt
if [ $? -eq 0 ]; then
    echo "Glucose Insulin Benchmark tests passed"
else
    echo "Glucose Insulin Benchmark tests failed"
    exit 1
fi

# clean up
echo "Cleaning up"
rm -rf $BUILD_DIR

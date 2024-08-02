#!/bin/bash

# This script is used to run the benchmark tests 

BUILD_DIR="build"

mkdir -p $BUILD_DIR

git submodule update --init --recursive

# compile and run benchmark test for each model
# Build the project using CMake
cd $BUILD_DIR
cmake .. -DCMAKE_BUILD_TYPE=Release -DBENCHMARK_ENABLE_TESTING=OFF 
make
cd ..

# Run the benchmark tests
./$BUILD_DIR/benchmark_bc --benchmark_format=console --benchmark_out=benchmark_results.txt
if [ $? -eq 0 ]; then
    echo "Benchmark tests passed"
else
    echo "Benchmark tests failed"
    exit 1
fi
./$BUILD_DIR/benchmark_cv --benchmark_format=console --benchmark_out=benchmark_results.txt
if [ $? -eq 0 ]; then
    echo "Benchmark tests passed"
else
    echo "Benchmark tests failed"
    exit 1
fi
./$BUILD_DIR/benchmark_gi --benchmark_format=console --benchmark_out=benchmark_results.txt
if [ $? -eq 0 ]; then
    echo "Benchmark tests passed"
else
    echo "Benchmark tests failed"
    exit 1
fi

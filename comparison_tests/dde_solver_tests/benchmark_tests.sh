#!/bin/bash
# This script is used to run the benchmark tests 
BUILD_DIR="bench_build"
DATA_DIR="dde_solver_tests_data/bench_data"
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
./$BUILD_DIR/benchmark_bc --benchmark_format=console --benchmark_out=$DATA_DIR/bc_benchmark_results.txt
if [ $? -eq 0 ]; then
echo "Benchmark tests passed"
else
echo "Benchmark tests failed"
exit 1
fi



# clean up
echo "Cleaning up"
rm -rf $BUILD_DIR

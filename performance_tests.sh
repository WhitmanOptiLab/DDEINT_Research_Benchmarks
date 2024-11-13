#!/bin/bash

# This script is used to run performance tests on the application

# the compiler flags
CC="g++"
CFLAGS="-std=c++17 -g"
INCLUDES="-I./DDEINT"

# get all cpp files in the Tests directory
TESTS_SRC="Test_Cases/*.cpp"

# create necessary directories
BUILD_DIR="build"
PERF_DIR="perf_output"
PLOT_DIR="plots"
DATA_DIR="data"
mkdir -p $BUILD_DIR $PERF_DIR $PLOT_DIR $DATA_DIR

# make sure this is a linux machine
if [ "$(uname)" != "Linux" ]; then
    echo "This script is only supported on Linux due to the use of perf. Exiting."
    exit 1
fi

# Function to generate FlameGraph
generate_flamegraph() 
{
    local output_file=$1
    echo "Generating Flame Graph for $output_file"
    if [ -f $PERF_DIR/perf.data ]; then
        perf script -i $PERF_DIR/perf.data > $PERF_DIR/out.perf
        ./utils/FlameGraph/stackcollapse-perf.pl $PERF_DIR/out.perf > $PERF_DIR/out.folded
        ./utils/FlameGraph/flamegraph.pl --color=java --title="$output_file" $PERF_DIR/out.folded > $PLOT_DIR/$output_file.svg
        if [ -f $PLOT_DIR/$output_file.svg ]; then
            echo "Flame Graph saved to $PLOT_DIR/$output_file.svg"
        else
            echo "Failed to generate Flame Graph: $PLOT_DIR/$output_file.svg not found"
        fi
    else
        echo "Failed to generate Flame Graph: $PERF_DIR/perf.data not found"
    fi
}

# compile and run performance test for each model
compile_and_test() 
{
    local src_files=$1
    local output_file=$2
    echo "Compiling $output_file"
    $CC $CFLAGS $INCLUDES $src_files -o $BUILD_DIR/$output_file

    echo "Running $output_file Test"
    perf record -o $PERF_DIR/perf.data -g ./$BUILD_DIR/$output_file
    if [ $? -eq 0 ]; then
        echo "Perf record completed successfully"
        cd $PERF_DIR
        # save perf report to file in data directory
        perf report > ../$DATA_DIR/${output_file}_perf_report.txt
        echo "Perf report saved to $DATA_DIR/${output_file}_perf_report.txt"
        cd ..
    else
        echo "Perf record failed"
    fi

    generate_flamegraph $output_file
}

# Compile and run performance tests for each test file
for test_file in $TESTS_SRC
do
    test_name=$(basename $test_file .cpp)
    compile_and_test $test_file $test_name
done

# clean up
echo "Cleaning up"
rm -rf $BUILD_DIR $PERF_DIR
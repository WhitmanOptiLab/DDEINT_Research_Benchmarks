#!/bin/bash

# This script is used to run performance tests on the application

# the compiler flags
CC="g++"
CFLAGS="-std=c++17 -g"
INCLUDES="-IDDEINT"

# the source files
BC_SRC="Breast_Cancer_Model/bc_model_fi.cpp"
BC_OUTPUT="bc_model"

CV_SRC="Cardiovascular_Model/cv_model_fi.cpp"
CV_OUTPUT="cv_model"

GI_SRC="Glucose_Insulin_Model/gi_model_fi.cpp"
GI_OUTPUT="gi_model"

HE_SRC="2D_Heat_Equation/heat_equation.cpp"
HE_OUTPUT="heat_equation"

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
        ../utils/FlameGraph/stackcollapse-perf.pl $PERF_DIR/out.perf > $PERF_DIR/out.folded
        ../utils/FlameGraph/flamegraph.pl --color=java --title="$output_file" $PERF_DIR/out.folded > $PLOT_DIR/$output_file.svg
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
    time $CC $CFLAGS $INCLUDES $src_files -o $BUILD_DIR/$output_file

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

# Breast Cancer Model
compile_and_test "$BC_SRC" "$BC_OUTPUT"

# Cardiovascular Model
#compile_and_test "$CV_SRC" "$CV_OUTPUT"

# Glucose Insulin Model
#compile_and_test "$GI_SRC" "$GI_OUTPUT"

# 2D Heat Equation
# compile_and_test "$HE_SRC" "$HE_OUTPUT"

# clean up
echo "Cleaning up"
rm -rf $BUILD_DIR $PERF_DIR

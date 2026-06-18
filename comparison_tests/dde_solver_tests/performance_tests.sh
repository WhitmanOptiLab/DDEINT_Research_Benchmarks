#!/bin/bash
# This script is used to run performance tests on the application
# the compiler flags
FC="gfortran"
FFLAGS="-g"
DDE_SOLVER_SRC="../../comparison_libraries/dde_solver/dde_solver_m.f90"

# the source files
BC_SRC="Breast_Cancer_Model/define_ddes.f90 Breast_Cancer_Model/bc_model.f90"
BC_OUTPUT="bc_model"


# create necessary directories
BUILD_DIR="perf_build"
PERF_DIR="perf_output"
PLOT_DIR="dde_solver_tests_data/perf_plots"
DATA_DIR="dde_solver_tests_data/perf_data"
mkdir -p $BUILD_DIR $PERF_DIR $PLOT_DIR $DATA_DIR dde_solver_tests_data/csv_files

# make sure this is a linux machine
if [ "$(uname)" != "Linux" ]; then
    echo "This script is only supported on Linux due to the use of perf. Exiting."
    exit 1
fi

# Function to generate FlameGraph
generate_flamegraph() {
    local output_file=$1
    echo "Generating Flame Graph for $output_file"
    if [ -f $PERF_DIR/perf.data ]; then
        perf script -i $PERF_DIR/perf.data > $PERF_DIR/out.perf
        ../../utils/FlameGraph/stackcollapse-perf.pl $PERF_DIR/out.perf > $PERF_DIR/out.folded
        ../../utils/FlameGraph/flamegraph.pl --color=java --title="$output_file" $PERF_DIR/out.folded > $PLOT_DIR/$output_file.svg
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
compile_and_test() {
    local src_files=$1
    local output_file=$2
    echo "Compiling $output_file"
    time $FC $FFLAGS $DDE_SOLVER_SRC $src_files -o $BUILD_DIR/$output_file


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




# clean up
echo "Cleaning up"
rm -rf $BUILD_DIR $PERF_DIR
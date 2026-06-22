#!/bin/bash
# This script is used to run performance tests on the application
# the compiler flags
FC="gfortran"
FFLAGS="-g"
DDE_SOLVER_SRC="../../comparison_libraries/dde_solver/dde_solver_m.f90"

# the source files
BC_SRC="Breast_Cancer_Model/bc_define_ddes.f90 Breast_Cancer_Model/bc_model.f90"
BC_OUTPUT="bc_model"

CV_SRC="Cardiovascular_Model/cv_define_ddes.f90 Cardiovascular_Model/cv_model.f90"
CV_OUTPUT="cv_model"

CW_SRC="Cobweb_Model/cw_define_ddes.f90 Cobweb_Model/cw_model.f90"
CW_OUTPUT="cw_model"

GI_SRC="Glucose_Insulin_Model/gi_define_ddes.f90 Glucose_Insulin_Model/gi_model.f90"
GI_OUTPUT="gi_model"

# create necessary directories
BUILD_DIR="perf_build"
PERF_DIR="perf_output"
PERF_PLOT_DIR="data/perf_plots"
PERF_DATA_DIR="data/perf_data"
DATA_DIR="data/csv_files"

mkdir -p $BUILD_DIR $PERF_DIR $PERF_PLOT_DIR $PERF_DATA_DIR $DATA_DIR



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
        ../../utils/FlameGraph/flamegraph.pl --color=java --title="$output_file" $PERF_DIR/out.folded > $PERF_PLOT_DIR/$output_file.svg
        if [ -f $PERF_PLOT_DIR/$output_file.svg ]; then
            echo "Flame Graph saved to $PERF_PLOT_DIR/$output_file.svg"
        else
            echo "Failed to generate Flame Graph: $PERF_PLOT_DIR/$output_file.svg not found"
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
        perf report > ../$PERF_DATA_DIR/${output_file}_perf_report.txt
        echo "Perf report saved to $PERF_DATA_DIR/${output_file}_perf_report.txt"
        cd ..
    else
        echo "Perf record failed"
    fi

    generate_flamegraph $output_file
}

# Breast Cancer Model
compile_and_test "$BC_SRC" "$BC_OUTPUT"
# Cardiovascular Model
compile_and_test "$CV_SRC" "$CV_OUTPUT"
# Cobweb Model
compile_and_test "$CW_SRC" "$CW_OUTPUT"
# Glucose Insulin Model
compile_and_test "$GI_SRC" "$GI_OUTPUT"



# clean up for files that are needed !!
echo "Cleaning up"
rm -rf $BUILD_DIR $PERF_DIR
rm -f bc_define_ddes.mod cv_define_ddes.mod cw_define_ddes.mod gi_define_ddes.mod dde_solver_m.mod
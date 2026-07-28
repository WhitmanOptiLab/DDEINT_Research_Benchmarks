#include "../../utils/benchmark/include/benchmark/benchmark.h"

#include "cw_functions.hpp"
#include "../../DDEINT/Methods/Dormand_Prince/DoPri_5.hpp"

// Define the absolute and relative tolerances
#define ABS_TOL 1e-12
#define REL_TOL 1e-12

// Benchmark function specifically for the run function
static void DDEINT_BM_CW(benchmark::State& state) {
    // --- Struct params ---
    // No cells here: everything lives in "pub" (pri stays empty). num_cells
    // = 1 since this is one flat system, not a repeated cell pattern.
    SimStructParams<double> sp;
    sp.num_cells = 1;
    sp.pub_vars_per_cell = 2;
    sp.pri_vars_per_cell = 0;
    sp.cells_per_block = 1;
    sp.pub_max_delays = {cw_p.tau, cw_p.tau};
    sp.pri_max_delays = {};
    sp.pub_prehistory = {history_cw, history_cw};
    sp.pri_prehistory = {};
    sp.coefs_per_cell = 0;
    sp.coef_vals = {};

    // --- Solver params ---
    // Defines step sizing, initial conditions, and error boundaries for the adaptive solver.
    SimSolverParams<double> solp;
    solp.t_initial = 0.0;
    solp.initial_h = 0.1;
    solp.min_h = 1e-5;
    solp.pub_init = {0.4, 0.4};
    solp.pri_init = {};
    solp.abs_tol = ABS_TOL;
    solp.rel_tol = REL_TOL;

    // --- Runtime params ---
    // High-level constraints on execution and output.
    SimRuntimeParams<double> rp;
    rp.t_initial = 0.0;
    rp.t_final = 10000.0;
    rp.sample_points = 500.0;
    rp.max_steps_per_call = 10000.0;


    // Initial conditions and time span
    std::vector<double> u0 = {0.4, 0.4};
    double t_initial = 0.0;
    double t_final = 10000.0;


    for (auto _ : state) {
        
        Simulation<DDESim<cw_dde, double>, double> sim(sp, solp, rp, false);
        sim.initialize();
        Results results = sim.simulator.simulate_collect(rp);
        benchmark::DoNotOptimize(results);
    }
}

// Register the function as a benchmark
BENCHMARK(DDEINT_BM_CW)
    ->Unit(benchmark::kMillisecond)  // Measure time in milliseconds
    ->ReportAggregatesOnly(true)  // Only report the aggregated results
    ->UseRealTime();  // Use wall-clock time

// Run the benchmark
BENCHMARK_MAIN();
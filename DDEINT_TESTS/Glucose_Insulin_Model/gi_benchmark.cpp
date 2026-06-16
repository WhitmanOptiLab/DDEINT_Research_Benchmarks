#include "../../utils/benchmark/include/benchmark/benchmark.h"
#include "gi_functions.hpp"
#include "../../DDEINT/Methods/Dormand_Prince/DoPri_5.hpp"

// Define the absolute and relative tolerances
#define ABS_TOL 1e-9
#define REL_TOL 1e-9

// Benchmark function specifically for the run function
static void DDEINT_BM_GI(benchmark::State& state) {
    // Initial conditions and time span
    std::vector<double> u0 = {gi_p.G_b, gi_p.I_b};
    double t_initial = 0.0;
    double t_final = 200.0;

    // Create the DDE problem and solve it
    std::vector<std::function<double(double)>> prehistory = {history_gi, history_gi};
    std::vector<double> max_delays = {gi_p.tau, gi_p.tau}; // Ensure the size matches the number of equations


    for (auto _ : state) {
        DoPri_5<gi_dde> dde_solver(2, max_delays, prehistory);
        dde_solver.initialize(0, 0.1, 1e-5, u0, ABS_TOL, REL_TOL, false, false);
        Results results = dde_solver.solve(t_initial, t_final, 500, 10000);
        benchmark::DoNotOptimize(results);
    }
}

// Register the function as a benchmark
BENCHMARK(DDEINT_BM_GI)
    ->Unit(benchmark::kMillisecond)  // Measure time in milliseconds
    ->ReportAggregatesOnly(true)  // Only report the aggregated results
    ->UseRealTime();  // Use wall-clock time

// Run the benchmark
BENCHMARK_MAIN();


#include "../../utils/benchmark/include/benchmark/benchmark.h"
#include "bc_functions.hpp"
#include "../../DDEINT/Methods/Dormand_Prince/DoPri_5.hpp"

// Define the absolute and relative tolerances
#define ABS_TOL 1e-9
#define REL_TOL 1e-9

static void DDEINT_BM_BC(benchmark::State& state)
{
    // Initial conditions and time span
    std::vector<double> u0 = {1.0, 1.0, 1.0};
    double t_initial = 0.0;
    double t_final = 10.0;

    // Create the DDE problem and solve it
    std::vector<std::function<double(double)>> prehistory = {history_bc, history_bc, history_bc};
    std::vector<double> max_delays = {1.0, 1.0, 1.0}; // Ensure the size matches the number of equations

    for (auto _ : state)
    {

        DoPri_5<bc_dde> dde_solver(3, max_delays, prehistory);
        dde_solver.initialize(0, 0.1, 1e-5, u0, ABS_TOL, REL_TOL, false, false);
        Results results = dde_solver.solve(t_initial, t_final, 500, 20000);
        benchmark::DoNotOptimize(results);
    }
}

BENCHMARK(DDEINT_BM_BC)
    ->Unit(benchmark::kMillisecond)  // Measure time in milliseconds
    ->ReportAggregatesOnly(true)  // Only report the aggregated results
    ->UseRealTime();  // Use wall-clock time

BENCHMARK_MAIN();

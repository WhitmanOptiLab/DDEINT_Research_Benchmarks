#include "../../utils/benchmark/include/benchmark/benchmark.h"
#include "bc_functions.hpp"
#include "../../DDEINT/dopri/ddeint_dopri_5.hpp"

// Define the absolute and relative tolerances
#define ABS_TOL 1e-9
#define REL_TOL 1e-9

static void BM_DDEINT_BC(benchmark::State& state)
{
    // Initial conditions and time span
    std::vector<double> u0 = {1.0, 1.0, 1.0};
    double t_i = 0.0;
    double t_f = 10.0;

    // Create the DDE problem and solve it
    std::vector<std::function<double(double)>> prehistory = {history_bc, history_bc, history_bc};
    std::vector<double> max_delays = {1.0, 1.0, 1.0}; // Ensure the size matches the number of equations

    DDEint_dopri_5<bc_dde> dde_solver(3, max_delays, prehistory);

    for (auto _ : state)
    {
        std::vector<std::vector<double>> solution = dde_solver.run(t_i, t_f, u0, 0.1, 1e-5, 20000, ABS_TOL, REL_TOL);
        benchmark::DoNotOptimize(solution);
    }
}

BENCHMARK(BM_DDEINT_BC)
    ->Unit(benchmark::kMillisecond)  // Measure time in milliseconds
    ->ReportAggregatesOnly(true)  // Only report the aggregated results
    ->UseRealTime();  // Use wall-clock time

BENCHMARK_MAIN();

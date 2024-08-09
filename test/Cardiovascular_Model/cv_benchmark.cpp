#include "../../utils/benchmark/include/benchmark/benchmark.h"
#include "cv_functions.hpp"
#include "../../DDEINT/dopri/ddeint_dopri_5.hpp"

// Define the absolute and relative tolerances
#define ABS_TOL 1e-9
#define REL_TOL 1e-9

static void BM_DDEINT_CV(benchmark::State& state)
{
    // Initial conditions and time span
    std::vector<double> u0 = {93, (1 / (1 + cv_p.R / cv_p.r)) * 93, (1 / (cv_p.R * cv_p.Vstr)) * (1 / (1 + cv_p.r / cv_p.R)) * 93};
    double t_initial = 0.0;
    double t_final = 1000.0;

    // Create the DDE problem and solve it
    std::vector<std::function<double(double)>> prehistory = {history_Pa, history_Pv, history_H};
    std::vector<double> max_delays = {4.0, 4.0, 4.0}; // Ensure the size matches the number of equations

    DDEint_dopri_5<cv_dde> dde_solver(3, max_delays, prehistory);

    for (auto _ : state)
    {
        std::vector<std::vector<double>> solution = dde_solver.run(t_initial, t_final, u0, 0.1, 1e-5, 20000, ABS_TOL, REL_TOL);
        benchmark::DoNotOptimize(solution);
    }
}

BENCHMARK(BM_DDEINT_CV)
    ->Unit(benchmark::kMillisecond)  // Measure time in milliseconds
    ->ReportAggregatesOnly(true)  // Only report the aggregated results
    ->UseRealTime();  // Use wall-clock time

BENCHMARK_MAIN();

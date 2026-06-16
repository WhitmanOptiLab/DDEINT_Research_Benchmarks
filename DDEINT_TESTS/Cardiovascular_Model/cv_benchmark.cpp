#include "../../utils/benchmark/include/benchmark/benchmark.h"
#include "cv_functions.hpp"
#include "../../DDEINT/Methods/Dormand_Prince/DoPri_5.hpp"

// Define the absolute and relative tolerances
#define ABS_TOL 1e-9
#define REL_TOL 1e-9

static void DDEINT_BM_CV(benchmark::State& state)
{
    // Initial conditions and time span
    std::vector<double> u0 = {93, (1 / (1 + cv_p.R / cv_p.r)) * 93, (1 / (cv_p.R * cv_p.Vstr)) * (1 / (1 + cv_p.r / cv_p.R)) * 93};
    double t_initial = 0.0;
    double t_final = 1000.0;

    // Create the DDE problem and solve it
    std::vector<std::function<double(double)>> prehistory = {history_Pa, history_Pv, history_H};
    std::vector<double> max_delays = {4.0, 4.0, 4.0}; // Ensure the size matches the number of equations


    for (auto _ : state)
    {
        DoPri_5<cv_dde> dde_solver(3, max_delays, prehistory);
        dde_solver.initialize(0, 0.1, 1e-5, u0, ABS_TOL, REL_TOL, false, false);
        Results results = dde_solver.solve(t_initial, t_final, 500, 20000);
        benchmark::DoNotOptimize(results);
    }
}

BENCHMARK(DDEINT_BM_CV)
    ->Unit(benchmark::kMillisecond)  // Measure time in milliseconds
    ->ReportAggregatesOnly(true)  // Only report the aggregated results
    ->UseRealTime();  // Use wall-clock time

BENCHMARK_MAIN();

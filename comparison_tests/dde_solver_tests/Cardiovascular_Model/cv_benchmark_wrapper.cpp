#include "../../../utils/benchmark/include/benchmark/benchmark.h"
 
extern "C" void run_cv_solver_();
 
static void DDE_SOLVER_BM_CV(benchmark::State& state)
{
    for (auto _ : state)
    {
        run_cv_solver_();
    }
}
 
BENCHMARK(DDE_SOLVER_BM_CV)
    ->Unit(benchmark::kMillisecond)
    ->ReportAggregatesOnly(true)
    ->UseRealTime();
 
BENCHMARK_MAIN();

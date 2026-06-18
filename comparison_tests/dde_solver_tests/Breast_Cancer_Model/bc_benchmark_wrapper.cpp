#include "../../../utils/benchmark/include/benchmark/benchmark.h"
 
extern "C" void run_bc_solver_();
 
static void DDE_SOLVER_BM_BC(benchmark::State& state)
{
    for (auto _ : state)
    {
        run_bc_solver_();
    }
}
 
BENCHMARK(DDE_SOLVER_BM_BC)
    ->Unit(benchmark::kMillisecond)
    ->ReportAggregatesOnly(true)
    ->UseRealTime();
 
BENCHMARK_MAIN();

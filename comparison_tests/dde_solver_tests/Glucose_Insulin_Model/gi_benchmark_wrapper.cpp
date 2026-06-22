#include "../../../utils/benchmark/include/benchmark/benchmark.h"

extern "C" void run_gi_solver_();

static void DDE_SOLVER_BM_GI(benchmark::State& state)
{
    for (auto _ : state)
    {
        run_gi_solver_();
    }
}

BENCHMARK(DDE_SOLVER_BM_GI)
    ->Unit(benchmark::kMillisecond)
    ->ReportAggregatesOnly(true)
    ->UseRealTime();

BENCHMARK_MAIN();
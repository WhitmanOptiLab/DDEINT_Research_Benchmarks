#include "../../../utils/benchmark/include/benchmark/benchmark.h"

extern "C" void run_cw_solver_();

static void DDE_SOLVER_BM_CW(benchmark::State& state)
{
    for (auto _ : state)
    {
        run_cw_solver_();
    }
}

BENCHMARK(DDE_SOLVER_BM_CW)
    ->Unit(benchmark::kMillisecond)
    ->ReportAggregatesOnly(true)
    ->UseRealTime();

BENCHMARK_MAIN();
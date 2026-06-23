#!/bin/bash

# Run all benchmarks across all solvers

echo "Running DDEINT benchmarks..."
cd DDEINT_tests
bash benchmark_tests.sh
cd ..

echo "Running dde_solver benchmarks..."
cd comparison_tests/dde_solver_tests
bash benchmark_tests.sh
cd ../..

echo "Running dde (R) benchmarks..."
cd comparison_tests/dde_tests
bash benchmark_tests.sh
cd ../..

echo "All benchmarks complete. Results saved to:"
echo "  DDEINT_tests/data/bench_data"
echo "  comparison_tests/dde_solver_tests/data/bench_data"
echo "  comparison_tests/dde_tests/data/bench_data"

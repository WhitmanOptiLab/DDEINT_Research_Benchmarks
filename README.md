# DDEINT_Research_Benchmarks

This repository contains experimental cases to test the DDEINT library against other Delay Differential Equation (DDE) and Ordinary Differential Equation (ODE) solvers. The goal is to benchmark the performance and accuracy of DDEINT in various scenarios.

## Contributors

Luke Samuels, Terence Mahlatini, Uli Raudales, Sebastian Wiedenhoeft, Julio De Jesus

## Folder Structure
- `DDEINT_tests/` — benchmarks and performance tests for our library (DDEINT) 
- `comparison_tests/dde_solver_tests/` — equivalent tests for [dde_solver](https://github.com/WarrenWeckesser/dde_solver) (Fortran) 
- `comparison_tests/dde_tests/` — equivalent tests for [dde](https://github.com/mrc-ide/dde) (R-based DDE solver)
- `comparison_tests/dde_tests_compiled/` - equivalent tests using the same `dde` solver, but with model equations implemented as compiled C functions. This directory exists to isolate and measure the performance impact of R callback overhead while keeping the underlying solver unchanged.
- `comparison_tests/copasi_tests/` — validation tests comparing Copasi output against DDEINT output
- `comparison_libraries/` — git submodules for external solvers

Note: that [DDEINT](https://github.com/WhitmanOptiLab/DDEINT/tree/0265012fa5e706271ba44148650359d4c2869693) and [dde_solver](https://github.com/WarrenWeckesser/dde_solver) must be initialized as git submodules before running their respective tests. See Set Up step 3.

## Set Up

1. Clone Repo
```
git clone https://github.com/WhitmanOptiLab/DDEINT_Research_Benchmarks.git
```
2. Change directory to the repository
```
cd DDEINT_RESEARCH_BENCHMARKS
```
3. Update the submodules
```
git submodule update --init --recursive
```
4. (optional) Switch DDEINT to a different branch

If you want to use a specific branch of DDEINT instead of `main`:
```
cd DDEINT
git checkout <branch-name>
cd ..
```

## DDEINT

The primary library being tested. Written in C++ and benchmarked using Google Benchmark

### Running Benchmark Tests

Uses Google Benchmark to measure solver speed across different models. To run test on DDEINT:

```
cd DDEINT_tests
bash benchmark_tests.sh
```

The results will be saved to `DDEINT_tests/data/bench_data`

### Running Performance Tests

Uses Linux `perf` and [FlameGraph](https://github.com/brendangregg/FlameGraph) to profile solver CPU performance and generate flame graphs for each model.

```
cd DDEINT_tests
bash performance_tests.sh
```
The results will be saved to `DDEINT_tests/data/perf_data` & `DDEINT_tests/data/perf_plots` & `DDEINT_tests/data/csv_files`
<!-- Comparison tests -->

## Running Comparison Libraries Tests

### dde_solver (Fortran)
 
Uses Google Benchmark and Linux `perf` with FlameGraph. Requires `gfortran` and [dde_solver](https://github.com/WarrenWeckesser/dde_solver) initialized as a git submodule (See Set Up step 3).
 
```
cd comparison_tests/dde_solver_tests
bash benchmark_tests.sh   # results → data/bench_data
bash performance_tests.sh # results → data/perf_data, data/perf_plots, data/csv_files
```
---
### dde (R)
 
Uses [microbenchmark](https://cran.r-project.org/package=microbenchmark) for benchmarking only — no performance profiling. Requires R with `dde` and `microbenchmark` packages installed.
 
```
cd comparison_tests/dde_tests
bash benchmark_tests.sh   # results → data/bench_data, data/csv_files
```
---
### dde (Compiled C)

Uses the same `dde` solver and `microbenchmark` package, but model equations are implemented as compiled C functions.

These tests are meant to measure the R callback overhead while keeping the solver unchanged.

```
cd comparison_tests/dde_tests_compiled
bash benchmark_tests.sh   # results → data/bench_data, data/csv_files
```

## Validation Against COPASI

These tests validate DDEINT's numerical accuracy by comparing its output against
[COPASI](http://copasi.org/), an established biochemical simulation tool.
Models validated: breast cancer (`bc`), cardiovascular (`cv`), and cobweb (`cw`). The reference COPASI CSVs are already committed — no need to re-run `.cps` files.

Running the script generates comparison plots into `copasi_tests/validation/`.

```
cd comparison_tests/copasi_tests
python3 -m venv .venv
source .venv/bin/activate
pip install matplotlib
python3 test_copasi_accuracy.py
```

## Validation COPASI (ddeint method and lsoda method)

These tests validate DDEINT Method accuracy by comparing its output directly against LSODA METHOD. Models validated: breast cancer (`bc`), cardiovascular (`cv`), and cobweb (`cw`). The reference COPASI CSVs are not present; however, their output is stored in `copasi_methods_tests/data` as CSV files. To view the COPASI files `.cps,`, they are found in COPASI_DELAYS.

Running the script generates comparison plots into `copasi_methods_tests/validation/`.

```
cd comparison_tests/copasi_tests
python3 -m venv .venv
source .venv/bin/activate
pip install matplotlib
python3 test_methods_accuracy.py
```

<!-- Running all benchmarks -->
## Running All Benchmarks

To run benchmarks across all solvers at once from repo root:
```
bash benchmark_all_tests.sh
```
Results will be saved:
- `DDEINT_tests/data/bench_data`
- `comparison_tests/dde_solver_tests/data/bench_data`
- `comparison_tests/dde_tests/data/bench_data`
- `comparison_tests/dde_tests_compiled/data/bench_data`

Note: There is no top-level performance test runner. Performance profiling must be run individually per solver.

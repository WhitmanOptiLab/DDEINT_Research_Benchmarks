# DDEINT_Research_Benchmarks

This repository contains experimental cases to test the DDEINT library against other Delay Differential Equation (DDE) and Ordinary Differential Equation (ODE) solvers. The goal is to benchmark the performance and accuracy of DDEINT in various scenarios.

## Contributors

Luke Samuels, Terence Mahlatini, Uli Raudales, Sebastian Wiedenhoeft, Julio De Jesus

## Folder Structure
- `DDEINT_tests/` — benchmarks and performance tests for our library (DDEINT) **REQUIRED FOLDER**
- `comparison_tests/dde_solver_tests/` — equivalent tests for [dde_solver](https://github.com/WarrenWeckesser/dde_solver) (Fortran) **REQUIRED FOLDER**
- `comparison_tests/dde_tests/` — equivalent tests for [dde](https://github.com/mrc-ide/dde) (R-based DDE solver)
- `comparison_libraries/` — git submodules for external solvers

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

## Comparison Libraries Tests

Same structure as DDEINT - will be interacting with `comparison_tests/dde_solver_tests` & `comparison_tests/dde_tests`

All results will be saved inside each respective folder under `/data`.

### Benchmarking
```
cd comparison_tests/dde_solver_tests
bash benchmark_tests.sh
```
or 
```
cd comparison_tests/dde_tests
bash benchmark_tests.sh
```
### Performance (Linux)
```
cd comparison_tests/dde_solver_tests
bash performance_tests.sh
```
or 
```
cd comparison_tests/dde_tests
bash performance_tests.sh
```

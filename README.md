# DDEINT_Research_Benchmarks

This repository contains experimental cases to test the DDEINT library against other Delay Differential Equation (DDE) and Ordinary Differential Equation (ODE) solvers. The goal is to benchmark the performance and accuracy of DDEINT in various scenarios.

## Contributors

Luke Samuels, Terence Mahlatini, Uli Raudales, Sebastian Wiedenhoeft, Julio De Jesus

## Set up

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

## Running Benchmark Tests

Uses Google Benchmark to measure solver speed across different models.

```
cd test
bash benchmark_tests.sh
```

The results will be saved to `test/data/`

## Running Performance Tests

Uses Linux `perf` and [FlameGraph](https://github.com/brendangregg/FlameGraph) to profile solver CPU performance and generate flame graphs for each model.

```
cd test
bash performance_tests.sh
```
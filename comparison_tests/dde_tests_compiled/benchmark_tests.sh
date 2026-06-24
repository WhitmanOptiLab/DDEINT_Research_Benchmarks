#!/bin/bash
set -e

DDE_INCLUDE="../../comparison_libraries/dde/inst/include"
DATA_DIR="data/bench_data"
CSV_DIR="data/csv_files"
mkdir -p "$DATA_DIR" "$CSV_DIR"

compile_model() {
  local model=$1
  echo "Compiling ${model}..."
  PKG_CFLAGS="-I${DDE_INCLUDE}" R CMD SHLIB "${model}.c" -o "${model}.so" 2>&1 \
    | grep -v "^using\|^clang\|^make" || true
}

cleanup_model() {
  local model=$1
  rm -f "${model}.o" "${model}.so"
}

# ─── Breast Cancer ─────────────────────────────────────────────────────────────
compile_model "Breast_Cancer_Model/bc_model"
Rscript - <<'EOF'
library(dde); library(microbenchmark)
dyn.load("Breast_Cancer_Model/bc_model.so")
pars <- c(p0=0.2, q0=0.3, v0=1.0, d0=5.0,
          p1=0.2, q1=0.3, v1=1.0, d1=1.0,
          d2=1.0, beta0=1.0, beta1=1.0, tau=1.0)
y0 <- c(1.0, 1.0, 1.0)
tt <- seq(0.0, 10000.0, length.out = 500)

# ── Save trajectory ──
res <- dopri(y0, tt, "bc_model", pars, dllname="bc_model",
             n_history=10000L, atol=1e-12, rtol=1e-12, return_history=FALSE)
dir.create("data/csv_files", recursive=TRUE, showWarnings=FALSE)
df <- data.frame(t=res[,1], u1=res[,2], u2=res[,3], u3=res[,4])
write.csv(df, "data/csv_files/bc_model_dde_c.csv", row.names=FALSE)
cat("Done. Results saved to data/csv_files/bc_model_dde_c.csv\n")

# ── Benchmark ──
mb <- microbenchmark(
  dopri(y0, tt, "bc_model", pars, dllname="bc_model",
        n_history=10000L, atol=1e-12, rtol=1e-12, return_history=FALSE),
  times = 100)
print(mb)
dir.create("data/bench_data", recursive=TRUE, showWarnings=FALSE)
writeLines(capture.output(print(mb)), "data/bench_data/bc_benchmark_results_c.txt")
cat("Benchmark saved to data/bench_data/bc_benchmark_results_c.txt\n")
dyn.unload("Breast_Cancer_Model/bc_model.so")
EOF
cleanup_model "Breast_Cancer_Model/bc_model"
echo "Breast Cancer benchmark passed"

# ─── Cardiovascular ────────────────────────────────────────────────────────────
compile_model "Cardiovascular_Model/cv_model"
Rscript - <<'EOF'
library(dde); library(microbenchmark)
dyn.load("Cardiovascular_Model/cv_model.so")
pars <- c(ca=1.55, cv=519, R=1.05, r=0.068, Vstr=67.9,
          alpha0=93, alphas=93, alphap=93, alphaH=0.84,
          beta0=7, betas=7, betap=7, betaH=1.17, gammaH=0, tau=4.0)
R0 <- 1.05; r0 <- 0.068; Pa0 <- 93
Pv0 <- (1 / (1 + R0 / r0)) * 93
H0  <- (1 / (R0 * 67.9)) * (1 / (1 + r0 / R0)) * 93
y0 <- c(Pa0, Pv0, H0)
tt <- seq(0.0, 10000.0, length.out = 500)

# ── Save trajectory ──
res <- dopri(y0, tt, "cv_model", pars, dllname="cv_model",
             n_history=10000L, atol=1e-12, rtol=1e-12, return_history=FALSE)
dir.create("data/csv_files", recursive=TRUE, showWarnings=FALSE)
df <- data.frame(t=res[,1], Pa=res[,2], Pv=res[,3], H=res[,4])
write.csv(df, "data/csv_files/cv_model_dde_c.csv", row.names=FALSE)
cat("Done. Results saved to data/csv_files/cv_model_dde_c.csv\n")

# ── Benchmark ──
mb <- microbenchmark(
  dopri(y0, tt, "cv_model", pars, dllname="cv_model",
        n_history=10000L, atol=1e-12, rtol=1e-12, return_history=FALSE),
  times = 100)
print(mb)
dir.create("data/bench_data", recursive=TRUE, showWarnings=FALSE)
writeLines(capture.output(print(mb)), "data/bench_data/cv_benchmark_results_c.txt")
cat("Benchmark saved to data/bench_data/cv_benchmark_results_c.txt\n")
dyn.unload("Cardiovascular_Model/cv_model.so")
EOF
cleanup_model "Cardiovascular_Model/cv_model"
echo "Cardiovascular benchmark passed"

# ─── Cobweb ────────────────────────────────────────────────────────────────────
compile_model "Cobweb_Model/cw_model"
Rscript - <<'EOF'
library(dde); library(microbenchmark)
dyn.load("Cobweb_Model/cw_model.so")
pars <- c(a=-0.5, b=1.5, c=0.8, d=0.0, tau=1.5, beta=0.6, speed=1.2)
y0 <- c(0.4, 0.4)
tt <- seq(0.0, 10000.0, length.out = 500)

# ── Save trajectory ──
res <- dopri(y0, tt, "cw_model", pars, dllname="cw_model",
             n_history=1000L, atol=1e-12, rtol=1e-12, return_history=FALSE)
dir.create("data/csv_files", recursive=TRUE, showWarnings=FALSE)
df <- data.frame(t=res[,1], price=res[,2], expected_price=res[,3])
write.csv(df, "data/csv_files/cw_model_dde_c.csv", row.names=FALSE)
cat("Done. Results saved to data/csv_files/cw_model_dde_c.csv\n")

# ── Benchmark ──
mb <- microbenchmark(
  dopri(y0, tt, "cw_model", pars, dllname="cw_model",
        n_history=1000L, atol=1e-12, rtol=1e-12, return_history=FALSE),
  times = 100)
print(mb)
dir.create("data/bench_data", recursive=TRUE, showWarnings=FALSE)
writeLines(capture.output(print(mb)), "data/bench_data/cw_benchmark_results_c.txt")
cat("Benchmark saved to data/bench_data/cw_benchmark_results_c.txt\n")
dyn.unload("Cobweb_Model/cw_model.so")
EOF
cleanup_model "Cobweb_Model/cw_model"
echo "Cobweb benchmark passed"

# ─── Glucose-Insulin ───────────────────────────────────────────────────────────
compile_model "Glucose_Insulin_Model/gi_model"
Rscript - <<'EOF'
library(dde); library(microbenchmark)
dyn.load("Glucose_Insulin_Model/gi_model.so")
pars <- c(V_G=10.0, V_I=5.0, S_I=1.0, G_b=90.0,
          I_b=10.0, n=0.1, gamma=0.05, tau=5.0)
y0 <- c(90.0, 10.0)
tt <- seq(0.0, 10000.0, length.out = 500)

# ── Save trajectory ──
res <- dopri(y0, tt, "gi_model", pars, dllname="gi_model",
             n_history=1000L, atol=1e-12, rtol=1e-12, return_history=FALSE)
dir.create("data/csv_files", recursive=TRUE, showWarnings=FALSE)
df <- data.frame(t=res[,1], G=res[,2], I=res[,3])
write.csv(df, "data/csv_files/gi_model_dde_c.csv", row.names=FALSE)
cat("Done. Results saved to data/csv_files/gi_model_dde_c.csv\n")

# ── Benchmark ──
mb <- microbenchmark(
  dopri(y0, tt, "gi_model", pars, dllname="gi_model",
        n_history=1000L, atol=1e-12, rtol=1e-12, return_history=FALSE),
  times = 100)
print(mb)
dir.create("data/bench_data", recursive=TRUE, showWarnings=FALSE)
writeLines(capture.output(print(mb)), "data/bench_data/gi_benchmark_results_c.txt")
cat("Benchmark saved to data/bench_data/gi_benchmark_results_c.txt\n")
dyn.unload("Glucose_Insulin_Model/gi_model.so")
EOF
cleanup_model "Glucose_Insulin_Model/gi_model"
echo "Glucose-Insulin benchmark passed"

echo ""
echo "All dde C benchmarks complete. Results saved to:"
echo "  ${CSV_DIR}"
echo "  ${DATA_DIR}"
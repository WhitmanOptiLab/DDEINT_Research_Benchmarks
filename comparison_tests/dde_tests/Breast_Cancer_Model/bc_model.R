library(dde)
library(microbenchmark)

# ─── Parameters ───────────────────────────────────────────────────────────────
pars <- c(
  p0 = 0.2, q0 = 0.3, v0 = 1.0, d0 = 5.0,
  p1 = 0.2, q1 = 0.3, v1 = 1.0, d1 = 1.0,
  d2 = 1.0, beta0 = 1.0, beta1 = 1.0, tau = 1.0
)

# ─── Model ────────────────────────────────────────────────────────────────────
bc_model <- function(t, y, pars) {
  hist3 <- ylag(t - pars["tau"], 3L)
  du1 <- (pars["v0"] / (1 + pars["beta0"] * hist3^2)) * (pars["p0"] - pars["q0"]) * y[1] - pars["d0"] * y[1]
  du2 <- (pars["v0"] / (1 + pars["beta0"] * hist3^2)) * (1 - pars["p0"] + pars["q0"]) * y[1] +
         (pars["v1"] / (1 + pars["beta1"] * hist3^2)) * (pars["p1"] - pars["q1"]) * y[2] - pars["d1"] * y[2]
  du3 <- (pars["v1"] / (1 + pars["beta1"] * hist3^2)) * (1 - pars["p1"] + pars["q1"]) * y[2] - pars["d2"] * y[3]
  c(du1, du2, du3)
}

# ─── Initial conditions ───────────────────────────────────────────────────────
y0 <- c(1.0, 1.0, 1.0)
tt <- seq(0.0, 10.0, length.out = 500)

# ─── Solve ────────────────────────────────────────────────────────────────────
res <- dopri(y0, tt, bc_model, pars,
             n_history = 1000L,
             atol = 1e-9, rtol = 1e-9,
             return_history = FALSE)

# ─── Save to CSV ──────────────────────────────────────────────────────────────
dir.create("data/csv_files", recursive = TRUE, showWarnings = FALSE)
df <- data.frame(t = res[,1], u1 = res[,2], u2 = res[,3], u3 = res[,4])
write.csv(df, "data/csv_files/bc_model_dde.csv", row.names = FALSE)
cat("Done. Results saved to data/csv_files/bc_model_dde.csv\n")

# ─── Benchmark ────────────────────────────────────────────────────────────────
mb <- microbenchmark(
  dopri(y0, tt, bc_model, pars,
        n_history = 1000L,
        atol = 1e-9, rtol = 1e-9,
        return_history = FALSE),
  times = 100
)

print(mb)
dir.create("data/bench_data", recursive = TRUE, showWarnings = FALSE)
write.csv(summary(mb), "data/bench_data/bc_benchmark_results.csv", row.names = FALSE)
cat("Benchmark saved to data/bench_data/bc_benchmark_results.csv\n")
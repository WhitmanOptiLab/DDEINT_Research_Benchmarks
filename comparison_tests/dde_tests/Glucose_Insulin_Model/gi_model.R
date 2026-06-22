library(dde)
library(microbenchmark)

# ─── Parameters ───────────────────────────────────────────────────────────────
pars <- c(
  V_G   = 10.0,
  V_I   =  5.0,
  S_I   =  1.0,
  G_b   = 90.0,
  I_b   = 10.0,
  n     =  0.1,
  gamma =  0.05,
  tau   =  5.0
)

# ─── Model ────────────────────────────────────────────────────────────────────
gi_model <- function(t, y, pars) {
  G_tau <- ylag(t - pars["tau"], 1L)
  I_tau <- ylag(t - pars["tau"], 2L)

  dG <- (y[1] - pars["G_b"]) / pars["V_G"] -
        (pars["S_I"] * I_tau * y[1]) / pars["V_G"] + 10.0

  dI <- -pars["n"] * (y[2] - pars["I_b"]) +
         (pars["gamma"] * G_tau) / pars["V_I"]

  c(dG, dI)
}

# ─── Initial conditions ───────────────────────────────────────────────────────
y0 <- c(90.0, 10.0)
tt <- seq(0.0, 200.0, length.out = 500)

# ─── Solve ────────────────────────────────────────────────────────────────────
res <- dopri(y0, tt, gi_model, pars,
             n_history = 1000L,
             atol = 1e-9, rtol = 1e-9,
             return_history = FALSE)

# ─── Save to CSV ──────────────────────────────────────────────────────────────
dir.create("data/csv_files", recursive = TRUE, showWarnings = FALSE)
df <- data.frame(
  t = res[, 1],
  G = res[, 2],
  I = res[, 3]
)
write.csv(df, "data/csv_files/gi_model_dde.csv", row.names = FALSE)
cat("Done. Results saved to data/csv_files/gi_model_dde.csv\n")

# ─── Benchmark ────────────────────────────────────────────────────────────────
mb <- microbenchmark(
  dopri(y0, tt, gi_model, pars,
        n_history = 1000L,
        atol = 1e-9, rtol = 1e-9,
        return_history = FALSE),
  times = 100
)

print(mb)
dir.create("data/bench_data", recursive = TRUE, showWarnings = FALSE)
write.csv(summary(mb), "data/bench_data/gi_benchmark_results.csv", row.names = FALSE)
cat("Benchmark saved to data/bench_data/gi_benchmark_results.csv\n")
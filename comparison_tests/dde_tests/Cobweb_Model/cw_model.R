library(dde)
library(microbenchmark)

# ─── Parameters ───────────────────────────────────────────────────────────────
pars <- c(
  a     = -0.5,
  b     =  1.5,
  c     =  0.8,
  d     =  0.0,
  tau   =  1.5,
  beta  =  0.6,
  speed =  1.2
)

# ─── Model ────────────────────────────────────────────────────────────────────
cw_model <- function(t, y, pars) {
  p_now              <- y[1]
  p_expected_delayed <- ylag(t - pars["tau"], 2L)

  demand <- pars["a"] * p_now + pars["b"]
  supply <- pars["c"] * p_expected_delayed + pars["d"]

  dp  <- pars["speed"] * (demand - supply)
  dpe <- pars["beta"]  * (p_now - y[2])

  c(dp, dpe)
}

# ─── Initial conditions ───────────────────────────────────────────────────────
y0 <- c(0.4, 0.4)
tt <- seq(0.0, 1000.0, length.out = 500)

# ─── Solve ────────────────────────────────────────────────────────────────────
res <- dopri(y0, tt, cw_model, pars,
             n_history = 1000L,
             atol = 1e-9, rtol = 1e-9,
             return_history = FALSE)

# ─── Save to CSV ──────────────────────────────────────────────────────────────
dir.create("data/csv_files", recursive = TRUE, showWarnings = FALSE)
df <- data.frame(
  t                    = res[, 1],
  price_values         = res[, 2],
  expected_price_values = res[, 3]
)
write.csv(df, "data/csv_files/cw_model_dde.csv", row.names = FALSE)
cat("Done. Results saved to data/csv_files/cw_model_dde.csv\n")

# ─── Benchmark ────────────────────────────────────────────────────────────────
mb <- microbenchmark(
  dopri(y0, tt, cw_model, pars,
        n_history = 1000L,
        atol = 1e-9, rtol = 1e-9,
        return_history = FALSE),
  times = 100
)

print(mb)
dir.create("data/bench_data", recursive = TRUE, showWarnings = FALSE)
write.csv(summary(mb), "data/bench_data/cw_benchmark_results.csv", row.names = FALSE)
cat("Benchmark saved to data/bench_data/cw_benchmark_results.csv\n")
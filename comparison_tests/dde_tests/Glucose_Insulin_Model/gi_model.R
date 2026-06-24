library(microbenchmark)
source("Glucose_Insulin_Model/gi_functions.R")

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
        atol = 1e-12, rtol = 1e-12,
        return_history = FALSE),
  times = 100
)
print(mb)
dir.create("data/bench_data", recursive = TRUE, showWarnings = FALSE)
write.csv(summary(mb), "data/bench_data/gi_benchmark_results.txt", row.names = FALSE)
cat("Benchmark saved to data/bench_data/gi_benchmark_results.txt\n")
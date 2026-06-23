library(microbenchmark)
source("Cardiovascular_Model/cv_functions.R")

# ─── Solve ────────────────────────────────────────────────────────────────────
res <- dopri(y0, tt, cv_model, pars,
             n_history = 1000L,
             atol = 1e-9, rtol = 1e-9,
             return_history = FALSE)

# ─── Save to CSV ──────────────────────────────────────────────────────────────
dir.create("data/csv_files", recursive = TRUE, showWarnings = FALSE)
df <- data.frame(t = res[,1], Pa = res[,2], Pv = res[,3], H = res[,4])
write.csv(df, "data/csv_files/cv_model_dde.csv", row.names = FALSE)
cat("Done. Results saved to data/csv_files/cv_model_dde.csv\n")

# ─── Benchmark ────────────────────────────────────────────────────────────────
mb <- microbenchmark(
  dopri(y0, tt, cv_model, pars,
        n_history = 1000L,
        atol = 1e-9, rtol = 1e-9,
        return_history = FALSE),
  times = 100
)
print(mb)
dir.create("data/bench_data", recursive = TRUE, showWarnings = FALSE)
write.csv(summary(mb), "data/bench_data/cv_benchmark_results.csv", row.names = FALSE)
cat("Benchmark saved to data/bench_data/cv_benchmark_results.csv\n")
library(microbenchmark)
source("Breast_Cancer_Model/bc_functions.R")

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
        n_history = 10000L,
        atol = 1e-12, rtol = 1e-12,
        return_history = FALSE),
  times = 100
)
print(mb)
dir.create("data/bench_data", recursive = TRUE, showWarnings = FALSE)
write.csv(summary(mb), "data/bench_data/bc_benchmark_results.txt", row.names = FALSE)
cat("Benchmark saved to data/bench_data/bc_benchmark_results.txt\n")
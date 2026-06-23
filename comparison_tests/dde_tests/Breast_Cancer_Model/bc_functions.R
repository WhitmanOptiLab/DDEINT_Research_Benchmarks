library(dde)

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

# ─── Initial Conditions ───────────────────────────────────────────────────────
y0 <- c(1.0, 1.0, 1.0)
tt <- seq(0.0, 10000.0, length.out = 500)
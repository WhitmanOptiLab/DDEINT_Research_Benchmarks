library(dde)

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

# ─── Initial Conditions ───────────────────────────────────────────────────────
y0 <- c(90.0, 10.0)
tt <- seq(0.0, 200.0, length.out = 500)
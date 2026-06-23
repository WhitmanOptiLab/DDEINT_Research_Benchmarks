library(dde)

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

# ─── Initial Conditions ───────────────────────────────────────────────────────
y0 <- c(0.4, 0.4)
tt <- seq(0.0, 10000.0, length.out = 500)
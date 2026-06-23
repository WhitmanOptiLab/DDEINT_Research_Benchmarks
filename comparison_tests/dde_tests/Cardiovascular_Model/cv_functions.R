library(dde)

# ─── Parameters ───────────────────────────────────────────────────────────────
pars <- c(
  ca     = 1.55,
  cv     = 519,
  R      = 1.05,
  r      = 0.068,
  Vstr   = 67.9,
  alpha0 = 93,
  alphas = 93,
  alphap = 93,
  alphaH = 0.84,
  beta0  = 7,
  betas  = 7,
  betap  = 7,
  betaH  = 1.17,
  gammaH = 0,
  tau    = 4.0
)

# ─── Model ────────────────────────────────────────────────────────────────────
cv_model <- function(t, y, pars) {
  R <- if (t <= 600) 1.05 else 0.21 * exp(600 - t) + 0.84
  Patau <- ylag(t - pars["tau"], 1L)
  Pa    <- y[1]
  Pv    <- y[2]
  H     <- y[3]
  dPa <- -(1 / (pars["ca"] * R)) * Pa +
           (1 / (pars["ca"] * R)) * Pv +
           (1 / pars["ca"]) * pars["Vstr"] * H
  dPv <- (1 / (pars["cv"] * R)) * Pa -
         (1 / (pars["cv"] * R) + 1 / (pars["cv"] * pars["r"])) * Pv
  Ts  <- 1 / (1 + (Patau / pars["alphas"])^pars["betas"])
  Tp  <- 1 / (1 + (pars["alphap"] / Pa)^pars["betap"])
  dH  <- (pars["alphaH"] * Ts) / (1 + pars["gammaH"] * Tp) - pars["betaH"] * Tp
  c(dPa, dPv, dH)
}

# ─── Initial Conditions ───────────────────────────────────────────────────────
R0  <- 1.05
r0  <- 0.068
Pa0 <- 93
Pv0 <- (1 / (1 + R0 / r0)) * 93
H0  <- (1 / (R0 * 67.9)) * (1 / (1 + r0 / R0)) * 93
y0 <- c(Pa0, Pv0, H0)
tt <- seq(0.0, 10000.0, length.out = 500)
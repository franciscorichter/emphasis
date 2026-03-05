# =============================================================================
# Replication — Richter & Wit: GAM surrogate for species diversification
#
# Two experiments from the paper:
#
#  Example 1 — LDD (linear diversity-dependent speciation, model = "dd")
#    Paper:  lambda_t = lambda_0 - beta_N * N_t   (beta_N > 0 in paper)
#    Our fw: lambda   = beta_0   + beta_N  * N_t   (beta_N < 0 in our model)
#    Grid: mu_0 in [0.1,0.5], lambda_0 in [0.55,1.5], beta_N in [2e-4,5e-3]
#    Ground truth: analytical survival from Kolmogorov backward equations
#    Figures: r1–r6 (6 resolution scatter panels), cond.png (surface)
#    Table:   computational cost and mean absolute error vs M
#
#  Example 2 — LPD (N + cumulative PD dependence, model = c(1,1,0))
#    Grid: lambda_0 in [0.2,2], mu_0 in [0.05,0.8],
#          beta_N in [-0.3,0], beta_P in [-0.05,0.05]
#    Ground truth: empirical (100 reps per validation point)
#    Figure: study2 (scatter GAM vs empirical)
# =============================================================================

# =============================================================================
# Install / load emphasis from GitHub (works on any machine)
# =============================================================================
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
if (!requireNamespace("emphasis", quietly = TRUE) ||
    packageVersion("emphasis") < "0.0.0.9000") {
  remotes::install_github("franciscorichter/emphasis", quiet = TRUE)
}

for (pkg in c("mgcv", "ggplot2", "deSolve", "patchwork")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

library(emphasis)
library(mgcv)
library(ggplot2)
library(deSolve)
library(patchwork)

OUTDIR <- "/tmp"
MAX_T  <- 10     # crown age (both examples)
N_PRED <- 50     # heatmap grid resolution
set.seed(42)


# =============================================================================
# Analytical survival: Kolmogorov backward equations for LDD
#
# The C++ simulation tracks two crown sub-clades (N1, N2 independently labelled)
# and declares "extinct" when EITHER goes to 0 (line 232 of general_tree.hpp:
#   if (N1 < 1 || N2 < 1) { break_type = extinction; break; }).
# So the simulation computes P(both crown lineages survive to T).
#
# Under independence (valid for weak DD where beta_N is small):
#   P(both survive) = P(single lineage survives)^2
#
# We compute P(single lineage survives from N_0=1 over T time) via the 1D
# Kolmogorov backward equation, then square it.
#
# Uses OUR sign convention: beta_N <= 0 means lambda decreases with N.
# per-lineage speciation = max(0, lambda0 + beta_N * n)   (beta_N <= 0)
# per-lineage extinction = mu0                             (constant)
# a(n) = max(0, lambda0 + beta_N * n) * n   total speciation from n species
# b(n) = mu0 * n                            total extinction from n species
#
# Backward ODE (s = T - t = time remaining):
#   dP_n/ds = a(n)*(P_{n+1}-P_n) + b(n)*(P_{n-1}-P_n)
#   P_0 = 0 (absorbing),  P_{max_N} = 1 (boundary: certain survival)
#   Initial P_n(s=0) = I(n >= 1).  Solve to s=T, read P[n=1].
# =============================================================================

dd_survival_prob <- function(lambda0, mu0, beta_N, max_t, max_N = 100L) {
  # P(both crown lineages survive) = P(single lineage survives from N_0=1)^2
  nn    <- 0L:max_N
  a_vec <- pmax(0, lambda0 + beta_N * nn) * nn   # beta_N <= 0 for DD
  b_vec <- mu0 * nn

  rhs <- function(s, P, parms) {
    dP    <- numeric(max_N + 1L)
    n_seq <- seq_len(max_N - 1L)
    i_seq <- n_seq + 1L
    dP[i_seq] <- a_vec[i_seq] * (P[i_seq + 1L] - P[i_seq]) +
                 b_vec[i_seq] * (P[i_seq - 1L] - P[i_seq])
    list(dP)
  }

  P0  <- c(0, rep(1, max_N))
  sol <- deSolve::ode(y = P0, times = c(0, max_t), func = rhs, parms = NULL,
                       method = "lsoda", rtol = 1e-5, atol = 1e-7)
  # P(n=1 survives from N_0=1) = sol[2, 3]  (col 1 = time, col 2 = P(n=0), col 3 = P(n=1))
  p_single <- as.numeric(sol[2L, 3L])
  p_single^2   # crown conditioning: both crown lineages must survive
}

dd_survival_grid <- function(lambda0, mu0, beta_N, max_t, max_N = 100L) {
  n <- length(lambda0)
  cat("  Computing", n, "analytical survival probabilities...\n")
  pb  <- utils::txtProgressBar(min = 0, max = n, style = 3)
  out <- vapply(seq_len(n), function(i) {
    utils::setTxtProgressBar(pb, i)
    dd_survival_prob(lambda0[i], mu0[i], beta_N[i], max_t, max_N)
  }, numeric(1L))
  close(pb); cat("\n")
  out
}


# =============================================================================
# EXAMPLE 1 — LDD
# =============================================================================
cat("\n========== EXAMPLE 1: LDD model ==========\n\n")

LAM_MIN <- 0.55;  LAM_MAX <- 1.5
MU_MIN  <- 0.10;  MU_MAX  <- 0.5
# Paper uses beta_N in [2e-4, 5e-3] with POSITIVE sign (lambda decreases).
# We store as NEGATIVE (our convention).
BN_MAG_MIN <- 2e-4;  BN_MAG_MAX <- 5e-3

# ── Validation set: 100 points, analytical ground truth
set.seed(7)
n_val <- 100
val1 <- data.frame(
  lambda0 = runif(n_val, LAM_MIN,  LAM_MAX),
  betaN   = -runif(n_val, BN_MAG_MIN, BN_MAG_MAX),  # negative: our convention
  mu0     = runif(n_val, MU_MIN,   MU_MAX)
)

cat("Analytical ground truth for validation set:\n")
val1$p_analytical <- dd_survival_grid(
  val1$lambda0, val1$mu0, val1$betaN, MAX_T)

cat("  Analytical range: [", round(min(val1$p_analytical), 3),
    ",", round(max(val1$p_analytical), 3), "]\n")

# ── Six resolution levels
M_vals  <- c(300, 800, 2500, 8000, 25000, 80000)
M_labels <- formatC(M_vals, format = "d", big.mark = ",")

plots_r   <- vector("list", 6L)
cost_rows <- vector("list", 6L)
last_gam1 <- NULL   # keep the largest-M GAM for cond.png

for (j in seq_along(M_vals)) {
  M <- M_vals[j]
  cat("\n── M =", M_labels[j], "──\n")

  set.seed(j * 13L)
  train_mat <- cbind(
    beta_0  = runif(M, LAM_MIN,  LAM_MAX),
    beta_N  = -runif(M, BN_MAG_MIN, BN_MAG_MAX),  # negative!
    gamma_0 = runif(M, MU_MIN,   MU_MAX),
    gamma_N = 0
  )

  t_sim <- system.time(
    sims <- simulate_tree(pars = train_mat, max_t = MAX_T,
                          model = "dd", max_tries = 0, useDDD = FALSE,
                          max_lin = 1e4)
  )["elapsed"]

  # survived = "done" (completed) or "too_large" (clade survived, just huge)
  status   <- sapply(sims, `[[`, "status")
  survived <- as.integer(status %in% c("done", "too_large"))
  cat("  Survival rate:", round(mean(survived), 3),
      " [done:", sum(status == "done"), "| too_large:",
      sum(status == "too_large"), "| extinct:", sum(status == "extinct"), "]\n")

  df_tr <- data.frame(lambda0  = train_mat[, "beta_0"],
                      betaN    = train_mat[, "beta_N"],
                      mu0      = train_mat[, "gamma_0"],
                      survived = survived)

  t_fit <- system.time({
    # Univariate GAM
    gam_uni <- mgcv::gam(survived ~ s(lambda0) + s(mu0) + s(betaN),
                          data = df_tr, family = binomial())
    # Bivariate GAM (tensor-product smooths, k=5 per dim to control memory)
    gam_biv <- mgcv::gam(
      survived ~ te(lambda0, mu0, k = 5) + te(lambda0, betaN, k = 5) +
                 te(mu0, betaN, k = 5),
      data = df_tr, family = binomial())
  })["elapsed"]

  p_uni <- predict(gam_uni, val1, type = "response")
  p_biv <- predict(gam_biv, val1, type = "response")
  err_u <- mean(abs(p_uni - val1$p_analytical))
  err_b <- mean(abs(p_biv - val1$p_analytical))

  cat("  Error uni:", round(err_u, 4), "| Error biv:", round(err_b, 4),
      " | Fit time:", round(t_fit, 1), "s\n")

  cost_rows[[j]] <- data.frame(
    M          = M,
    sim_s      = round(t_sim, 1),
    fit_s      = round(t_fit, 1),
    error_uni  = round(err_u, 4),
    error_biv  = round(err_b, 4)
  )

  plt_df <- data.frame(
    analytical = rep(val1$p_analytical, 2L),
    predicted  = c(p_biv, p_uni),
    Splines    = rep(c("Bivariate", "Univariate"), each = n_val)
  )
  plots_r[[j]] <- ggplot(plt_df,
      aes(x = analytical, y = predicted, colour = Splines)) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    geom_point(alpha = 0.55, size = 1.0) +
    scale_colour_manual(
      values = c(Bivariate = "#1b7837", Univariate = "#2171b5"), name = NULL) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "Analytical", y = "GAM",
         title = paste0("M = ", M_labels[j])) +
    theme_minimal(base_size = 9) +
    theme(legend.position = "bottom",
          legend.text     = element_text(size = 7),
          legend.key.size = unit(0.4, "cm"),
          plot.title      = element_text(hjust = 0.5, size = 9))

  ggsave(file.path(OUTDIR, paste0("r", j, ".png")),
         plots_r[[j]], width = 4.5, height = 4.2, dpi = 150)
  last_gam1 <- gam_uni   # keep for cond.png
}

# ── 6-panel combined figure
fig_r6 <- patchwork::wrap_plots(plots_r, ncol = 2) +
  patchwork::plot_annotation(
    title   = "Example 1 — LDD: GAM vs analytical survival probability",
    caption = "Green = bivariate splines [te(x,y)]; Blue = univariate [s(x)]"
  )
ggsave(file.path(OUTDIR, "r1_to_r6.png"), fig_r6,
       width = 9, height = 13, dpi = 150)
cat("\nSaved r1_to_r6.png\n")

# ── Cost table
cost_df <- do.call(rbind, cost_rows)
cat("\nComputational cost (Example 1):\n")
print(cost_df, row.names = FALSE)

# ── Figure cond.png: survival surface (lambda0 × mu0), fixed |beta_N| = midpoint
BN_MID <- -exp(mean(log(c(BN_MAG_MIN, BN_MAG_MAX))))   # neg geometric midpoint
surf1  <- expand.grid(
  lambda0 = seq(LAM_MIN, LAM_MAX, length.out = N_PRED),
  betaN   = BN_MID,
  mu0     = seq(MU_MIN,  MU_MAX,  length.out = N_PRED)
)
surf1$p_survival <- predict(last_gam1, surf1, type = "response")

p_cond <- ggplot(surf1, aes(x = lambda0, y = mu0, fill = p_survival)) +
  geom_tile() +
  scale_fill_viridis_c(name = "P(survival)", limits = c(0, 1), option = "plasma") +
  labs(x = expression(lambda[0]), y = expression(mu[0]),
       title = "LDD: GAM survival probability surface",
       subtitle = bquote(beta[N] == .(signif(BN_MID, 2)) ~
                         "," ~ t[max] == .(MAX_T))) +
  theme_minimal(base_size = 13)
ggsave(file.path(OUTDIR, "cond.png"), p_cond, width = 7, height = 5.5, dpi = 150)
cat("Saved cond.png\n")


# =============================================================================
# EXAMPLE 2 — LPD  (model c(1,1,0): lambda = beta_0 + beta_N*N + beta_P*P)
# beta_N < 0 (DD), beta_P varies (PD dependence)
# =============================================================================
cat("\n========== EXAMPLE 2: LPD model ==========\n\n")

N_TRAIN2 <- 4500
N_VAL2   <- 100
N_REPS2  <- 100

set.seed(123)
train2 <- cbind(
  beta_0  = runif(N_TRAIN2, 0.2, 2.0),
  beta_N  = runif(N_TRAIN2, -0.3, -1e-4),   # negative: DD
  beta_P  = runif(N_TRAIN2, -0.05, 0.05),
  gamma_0 = runif(N_TRAIN2, 0.05, 0.8),
  gamma_N = 0,
  gamma_P = 0
)

cat("Simulating", N_TRAIN2, "LPD trees for training...\n")
t_sim2 <- system.time(
  sims2 <- simulate_tree(pars = train2, max_t = MAX_T,
                         model = c(1L, 1L, 0L), max_tries = 0,
                         useDDD = FALSE, max_lin = 1e4)
)["elapsed"]

status2   <- sapply(sims2, `[[`, "status")
survived2 <- as.integer(status2 %in% c("done", "too_large"))
cat("Status breakdown:\n"); print(table(status2))
cat("Survival rate:", round(mean(survived2), 3), "| sim time:", round(t_sim2, 1), "s\n")

df_tr2 <- data.frame(lambda0  = train2[, "beta_0"],
                     betaN    = train2[, "beta_N"],
                     betaP    = train2[, "beta_P"],
                     mu0      = train2[, "gamma_0"],
                     survived = survived2)
cat("Fitting univariate GAM (LPD)...\n")
t_fit2 <- system.time(
  gam2 <- mgcv::gam(survived ~ s(lambda0) + s(mu0) + s(betaN) + s(betaP),
                    data = df_tr2, family = binomial())
)["elapsed"]
cat("\nGAM summary (Example 2):\n"); print(summary(gam2))
cat("Fit time:", round(t_fit2, 1), "s\n")

# ── Empirical validation
set.seed(42)
val2 <- data.frame(lambda0 = runif(N_VAL2, 0.2, 2.0),
                   betaN   = runif(N_VAL2, -0.3, -1e-4),
                   betaP   = runif(N_VAL2, -0.05, 0.05),
                   mu0     = runif(N_VAL2, 0.05, 0.8))

cat("\nEmpirical validation at", N_VAL2, "points (", N_REPS2, "reps each)...\n")
pb <- utils::txtProgressBar(min = 0, max = N_VAL2, style = 3)
val2$p_empirical <- vapply(seq_len(N_VAL2), function(i) {
  utils::setTxtProgressBar(pb, i)
  pm <- matrix(c(val2$lambda0[i], val2$betaN[i], val2$betaP[i],
                 val2$mu0[i], 0, 0),
               nrow = N_REPS2, ncol = 6L, byrow = TRUE)
  reps <- simulate_tree(pars = pm, max_t = MAX_T,
                        model = c(1L, 1L, 0L), max_tries = 0,
                        useDDD = FALSE, max_lin = 1e4)
  mean(sapply(reps, `[[`, "status") %in% c("done", "too_large"))
}, numeric(1L))
close(pb); cat("\n")

val2$p_gam <- predict(gam2, val2, type = "response")
rmse2 <- sqrt(mean((val2$p_gam - val2$p_empirical)^2))
cat("RMSE GAM vs empirical:", round(rmse2, 4), "\n")

p_study2 <- ggplot(val2, aes(x = p_empirical, y = p_gam, colour = lambda0)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_colour_viridis_c(name = expression(lambda[0])) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Empirical P(survival)",
       y = "GAM P(survival)",
       title = "Example 2 — LPD: GAM vs empirical survival",
       subtitle = paste0(N_VAL2, " validation points,  RMSE = ",
                         round(rmse2, 3))) +
  theme_minimal(base_size = 13)
ggsave(file.path(OUTDIR, "study2.png"), p_study2,
       width = 6.5, height = 5.5, dpi = 150)
cat("Saved study2.png\n")

# =============================================================================
cat("\n── Summary ──────────────────────────────────────────────────────────\n")
cat("  r1.png .. r6.png    LDD panels (individual)\n")
cat("  r1_to_r6.png        LDD 6-panel combined figure\n")
cat("  cond.png            LDD survival-probability surface\n")
cat("  study2.png          LPD GAM vs empirical scatter\n")
cat("\n── Cost table (Example 1) ───────────────────────────────────────────\n")
print(cost_df, row.names = FALSE)

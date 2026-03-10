# emphasis — Evolutionary Modeling for Phylogenetic Inference

[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](./LICENSE)
![Lifecycle: Experimental](https://img.shields.io/badge/lifecycle-experimental-orange)

`emphasis` is an R package for studying species diversification and its ecological
drivers. It provides:

- **Simulation** — stochastic phylogenies under flexible diversification models
- **Inference** — parameter estimation via MCEM, Cross-Entropy Method (CEM), and GAM-based MLE
- **Model selection** — AIC-based comparison across competing models
- **Diagnostics** — convergence and importance-sampling quality checks

## Installation

```r
# If updating from a previous version, uninstall first to avoid stale cached code:
# remove.packages("emphasis")

devtools::install_github("franciscorichter/emphasis")
```

---

## The rate model

All models are special cases of a general per-lineage rate:

```
lambda(s,t) = f(beta_0  + beta_N*N(t) + beta_P*P(t) + beta_E*E(s,t))
mu(s,t)     = f(gamma_0 + gamma_N*N(t) + gamma_P*P(t) + gamma_E*E(s,t))
```

| Covariate | Meaning |
|-----------|---------|
| `N(t)` | Current species richness |
| `P(t)` | Total pendant branch-length (phylogenetic diversity) |
| `E(s,t)` | Time since lineage *s* last diverged (lineage-specific) |

The `link` argument controls `f`:

| `link` | `f(eta)` |
|--------|----------|
| `"linear"` (default) | `max(0, eta)` |
| `"exponential"` | `exp(eta)` |

> **Tip:** For DD, PD, and EP models, use `link = "exponential"` to avoid
> numerical instability.  With the linear link, covariate slopes can drive
> `lambda` to zero, causing IS collapse (most augmented trees get zero weight).

### Model shortcuts

The `model` argument selects active covariates. Use a string, formula, or binary vector:

| String | Formula | Active covariates | Parameters |
|--------|---------|-------------------|------------|
| `"cr"` | `~ 1` | none | `c(beta_0, gamma_0)` |
| `"dd"` | `~ N` | N | `c(beta_0, beta_N, gamma_0, gamma_N)` |
| `"pd"` | `~ PD` | P | `c(beta_0, beta_P, gamma_0, gamma_P)` |
| `"ep"` | `~ EP` | E | `c(beta_0, beta_E, gamma_0, gamma_E)` |
| — | `~ N + PD` | N, P | `c(beta_0, beta_N, beta_P, gamma_0, gamma_N, gamma_P)` |
| — | `~ N + PD + EP` | N, P, E | all 8 parameters |

The binary vector `c(use_N, use_P, use_E)` is also accepted. Parameter length is
always `2 + 2 * (number of active covariates)`.

> **EP note:** the evolutionary pendant model is supported in both simulation and
> inference (linear link only). The exponential link with EP has no closed-form
> integral and is not supported.

---

## Simulation

`simulate_tree()` handles both forward simulation and conditional simulation
(augmenting an observed tree with extinct lineages).

### Forward simulation

```r
library(emphasis)
set.seed(123)

cr <- simulate_tree(pars = c(0.5, 0.1),                         model = "cr", max_t = 5)
dd <- simulate_tree(pars = c(0.8, -0.02, 0.2, -0.005),          model = "dd", max_t = 5,
                    link = "exponential")
pd <- simulate_tree(pars = c(0.6, 0.005, 0.15, -0.003),         model = "pd", max_t = 5,
                    link = "exponential")
ep <- simulate_tree(pars = c(0.5, 0.05, 0.1, 0.01),             model = "ep", max_t = 5)
np <- simulate_tree(pars = c(0.5, -0.01, 0.005, 0.15, -0.002, 0.001),
                    model = ~ N + PD, max_t = 5, link = "exponential")

dd$status         # "done", "extinct", or "too_large"
dd$tes            # extant-only phylo
dd$survival_prob  # empirical survival probability (1 / n_attempts)
```

Each result contains: `tes` (extant phylo), `tas` (full phylo with extinct lineages),
`L` (L-table), `status`, `survival_prob`.

### Batch simulation

Pass a matrix to `pars` — one simulation per row:

```r
set.seed(1)
pars_mat <- cbind(beta_0  = runif(200, 0.3, 1.2),
                  gamma_0 = runif(200, 0.05, 0.5))
sims <- simulate_tree(pars = pars_mat, max_t = 8, model = "cr", max_tries = 0)
sims$survival_prob                          # fraction of rows that produced a tree
sims$simulations[[1]]$status               # "done", "extinct", or "too_large"
```

### Conditional simulation (augmentation)

Provide `tree =` to augment an observed extant tree with stochastically drawn
extinct lineages.

```r
set.seed(42)
obs <- simulate_tree(pars = c(0.6, 0.005, 0.15, -0.003), model = "pd", max_t = 5,
                     link = "exponential")

# Single augmented tree
aug <- simulate_tree(tree = obs, pars = c(0.6, 0.005, 0.15, -0.003),
                     model = "pd", n_trees = 1L, link = "exponential")
aug$tas    # augmented phylo (extant + stochastic extinct lineages)
aug$log_q  # log q(augmented tree | extant tree, pars)

# Many augmented trees
augN <- simulate_tree(tree = obs, pars = c(0.6, 0.005, 0.15, -0.003),
                      model = "pd", n_trees = 100L, link = "exponential")
augN$log_q       # log q values (length 100)
augN$trees[[1]]  # first augmented phylo
```

---

## Inference

`estimate_rates()` estimates parameters from a phylogenetic tree. It accepts a
`simulate_tree()` result, a `phylo` object, or a numeric branching-time vector.
All tuning parameters — including the required `lower_bound` and `upper_bound` —
are passed via `control = list(...)`. Inspect defaults with
`estimate_rates_control("cem")` or `estimate_rates_control("mcem")`.

### CEM

The Cross-Entropy Method performs a global search. No starting values are
needed; it initialises a population of parameter vectors randomly within the
bounds and concentrates them around the optimum via iterative selection and
perturbation. Stops when the population has collapsed (annealing exhausted),
the best fhat has plateaued, or `max_iter` is reached.

```r
set.seed(42)
sim <- simulate_tree(pars = c(0.6, -0.05, 0.15, -0.03), model = "pd", max_t = 10,
                     link = "exponential")

lb <- c(0.01, -0.5, 0, -0.5)
ub <- c(2,     0.5, 1,  0.5)

fit_cem <- estimate_rates(sim, method = "cem", model = "pd", link = "exponential",
  control = list(
    lower_bound   = lb,
    upper_bound   = ub,
    max_iter      = 30,
    num_particles = 80,     # parameter vectors per iteration
    num_trees     = 5,      # IS trees per particle; >= 5 enables bootstrap variance
    verbose       = TRUE
  ))

fit_cem$pars       # named estimates: beta_0, beta_P, gamma_0, gamma_P
fit_cem$loglik     # IS log-likelihood
fit_cem$loglik_var # MC variance of loglik (bootstrap, B=200); NA if num_trees = 1
fit_cem$AIC

# Convergence diagnostics: fhat trace, parameter traces, population spread, IS weights
diag_cem <- diagnose_cem(fit_cem, lower_bound = lb, upper_bound = ub)
diag_cem$convergence       # per-iteration: best_loglik, n_valid, rej_total
diag_cem$IS_quality        # ESS, ESS fraction, mean/sd of IS log-weights
```

| Control parameter | Default | Description |
|-------------------|---------|-------------|
| `lower_bound` | — | **Required.** Lower bounds vector |
| `upper_bound` | — | **Required.** Upper bounds vector |
| `verbose` | FALSE | Print iteration summaries |
| `max_iter` | 20 | Hard iteration cap |
| `num_particles` | 50 | Particles per iteration. Alias: `num_points` |
| `num_trees` | 1 | IS trees per particle per iteration. `>= 5` recommended for variance. Alias: `sample_size` |
| `shared_trees` | FALSE | Pool trees across particles (Mode 2; lower variance, slower) |
| `disc_prop` | 0.5 | Elite fraction of particles kept each iteration |
| `sd_vec` | NULL | Initial perturbation SDs; auto = `(upper - lower) / 4` |
| `tol` | 1e-4 | Minimum fhat improvement to reset the plateau counter |
| `patience` | 5 | Consecutive plateau iterations before stopping |
| `max_missing` | 1e4 | Max extinct lineages per augmented tree; trees exceeding this are discarded |
| `max_time` | 3600 | Max wall-clock seconds; NULL disables |
| `num_threads` | 1 | Parallel threads |

### MCEM

The Monte Carlo EM algorithm performs local refinement from a given starting
point. It iterates E-steps (augmenting trees at the current estimate) and
M-steps (maximising the IS-weighted Q-function via SBPLX). Convergence is
declared when all parameters stabilise: the maximum relative change (scaled
by the search range) stays below `tol` for `patience` consecutive iterations.

```r
fit_mcem <- estimate_rates(sim, method = "mcem", model = "pd", link = "exponential",
  init_pars = fit_cem$pars,   # warm-start from CEM (recommended)
  control = list(
    lower_bound = lb,
    upper_bound = ub,
    num_trees   = 200,   # IS trees per EM iteration; bootstrap variance auto-computed
    tol         = 1e-3,  # convergence: max relative parameter change < tol
    verbose     = TRUE
  ))

fit_mcem$pars
fit_mcem$loglik
fit_mcem$loglik_var  # bootstrap variance (B=200) from final E-step IS weights
fit_mcem$AIC

# Convergence diagnostics: fhat trace, parameter traces, parameter stability, IS weights
diag_mcem <- diagnose_mcem(fit_mcem, lower_bound = lb, upper_bound = ub)
diag_mcem$convergence  # per-iteration: fhat, delta_max, num_trees, time
diag_mcem$IS_quality   # ESS, mean/sd of log-weights, bootstrap variance
```

| Control parameter | Default | Description |
|-------------------|---------|-------------|
| `lower_bound` | — | **Required.** Lower bounds vector |
| `upper_bound` | — | **Required.** Upper bounds vector |
| `verbose` | FALSE | Print iteration summaries |
| `num_trees` | 200 | IS trees per EM iteration. Bootstrap variance auto-computed. Alias: `sample_size` |
| `max_iter` | 200 | Maximum EM iterations |
| `tol` | 1e-3 | Parameter-stability threshold: max relative change (fraction of search range) |
| `patience` | 3 | Consecutive stable iterations required to declare convergence |
| `xtol` | 1e-3 | Relative tolerance for the M-step optimiser |
| `maxN` | 10 | Max augmentation attempts per tree before rejection |
| `max_missing` | 1e4 | Max extinct lineages per augmented tree; trees exceeding this are discarded |
| `max_time` | 3600 | Max wall-clock seconds; NULL disables |
| `num_threads` | 1 | Parallel threads |

### GAM-based MLE

The GAM method evaluates the IS log-likelihood on a parameter grid, fits a
smooth GAM surface (via `mgcv`), and optimises it with L-BFGS-B. It is a
one-shot (non-iterative) approach — no starting values or convergence tuning
needed.

For models with 2 parameters (CR), use `grid_points` for a factorial grid.
For models with 3+ parameters (DD, PD, EP), use `n_grid` for a Latin
Hypercube Sampling grid to avoid the curse of dimensionality.

```r
# CR model (2 params) — factorial grid
fit_cr <- estimate_rates(sim, method = "gam", model = "cr",
  control = list(
    lower_bound  = c(0.1, 0.01),
    upper_bound  = c(2.0, 0.5),
    grid_points  = 15,       # 15^2 = 225 grid points
    sample_size  = 200,
    verbose      = TRUE
  ))

# DD model (4 params, gamma_N fixed) — LHS grid with survival conditioning
gam_surv <- train_GAM(sims, pars_mat, model = "dd")  # see GAM-based methods
fit_dd <- estimate_rates(sim, method = "gam", model = "dd", cond = gam_surv,
  control = list(
    lower_bound  = c(0.5, -0.15, 0.01, 0),
    upper_bound  = c(5.0, -0.001, 1.0, 0),  # gamma_N fixed at 0 (lb = ub)
    n_grid       = 200,      # 200 LHS points (replaces factorial grid)
    sample_size  = 50,
    verbose      = TRUE
  ))

fit_dd$pars
fit_dd$loglik    # conditioned log-likelihood
fit_dd$AIC

# Diagnostics: surface coverage, GAM fit quality, MLE location
diag_gam <- diagnose_gam(fit_dd)
```

| Control parameter | Default | Description |
|-------------------|---------|-------------|
| `lower_bound` | — | **Required.** Lower bounds vector |
| `upper_bound` | — | **Required.** Upper bounds vector |
| `verbose` | FALSE | Print progress |
| `grid_points` | 20 | Points per parameter axis (factorial grid). Use for <= 2 free params |
| `n_grid` | NULL | Total LHS grid points (overrides `grid_points`). Use for 3+ free params |
| `sample_size` | 200 | IS trees per grid point |
| `max_time` | 3600 | Max wall-clock seconds; NULL disables |
| `maxN` | NULL | Max augmentation attempts (NULL = auto) |
| `max_lambda` | 500 | Max speciation rate for thinning |
| `max_missing` | 1e4 | Max extinct lineages per augmented tree |
| `num_threads` | 1 | Parallel threads |

---

## Model selection

`compare_models()` takes two or more fitted models and returns a table sorted by
AIC. The model with the lowest AIC is preferred; `delta_AIC` shows the gap.

```r
set.seed(42)
sim <- simulate_tree(pars = c(0.8, -0.02, 0.2, 0), model = "dd", max_t = 8,
                     link = "exponential")

fit_cr <- estimate_rates(sim, method = "cem", model = "cr",
  control = list(lower_bound = c(0,    0),              upper_bound = c(2,    1),
                 max_iter = 15, num_particles = 40))
fit_dd <- estimate_rates(sim, method = "cem", model = "dd", link = "exponential",
  control = list(lower_bound = c(0.1, -0.1, 0, -0.01), upper_bound = c(2, 0.01, 0.5, 0.01),
                 max_iter = 15, num_particles = 40))
fit_pd <- estimate_rates(sim, method = "cem", model = "pd", link = "exponential",
  control = list(lower_bound = c(0.01, -1,  0, -1),    upper_bound = c(2,   1,   1,  1),
                 max_iter = 15, num_particles = 40))

compare_models(CR = fit_cr, DD = fit_dd, PD = fit_pd)
#   model n_pars  loglik    AIC delta_AIC
# 1    DD      4  -12.3   32.6      0.00
# 2    PD      4  -15.1   38.2      5.60
# 3    CR      2  -18.9   41.8      9.20
```

---

## Simulation study 1 — estimation power

Simulate 50 trees under a PD model, estimate parameters from each, and
evaluate accuracy via bias, RMSE, and tree-size effects.

```r
library(emphasis)

# True parameters
true_pars <- c(beta_0 = 0.6, beta_P = -0.05, gamma_0 = 0.15, gamma_P = -0.03)
model     <- "pd"
lnk       <- "exponential"
max_t     <- 10
n_rep     <- 50

ctrl <- list(lower_bound   = c(0.01, -0.5, 0,  -0.5),
             upper_bound   = c(2,     0.5, 1,   0.5),
             max_iter      = 100,
             num_particles = 100,
             num_trees     = 1)

# ── 1. Simulate trees ─────────────────────────────────────────────────────
set.seed(1)
trees <- vector("list", n_rep)
for (i in seq_len(n_rep)) {
  repeat {
    sim <- simulate_tree(pars = true_pars, model = model, max_t = max_t, link = lnk)
    if (sim$status == "done") { trees[[i]] <- sim; break }
  }
}

# ── 2. Estimate on each tree ───────────────────────────────────────────────
results  <- vector("list", n_rep)
n_tips   <- integer(n_rep)
for (i in seq_len(n_rep)) {
  cat(sprintf("Rep %d/%d\n", i, n_rep))
  n_tips[i] <- length(trees[[i]]$tes$tip.label)
  fit <- tryCatch(
    estimate_rates(trees[[i]], method = "cem", model = model, link = lnk, control = ctrl),
    warning = function(w) NULL,
    error   = function(e) NULL
  )
  if (!is.null(fit) && all(is.finite(fit$pars)))
    results[[i]] <- as.list(fit$pars)
}

# ── 3. Accuracy summary ────────────────────────────────────────────────────
ok      <- !sapply(results, is.null)
cat(sprintf("\nConverged: %d / %d\n\n", sum(ok), n_rep))

est_mat <- do.call(rbind, lapply(results[ok], as.numeric))
colnames(est_mat) <- names(true_pars)
tips_ok <- n_tips[ok]

bias <- colMeans(est_mat) - true_pars
rmse <- sqrt(colMeans((est_mat - rep(true_pars, each = nrow(est_mat)))^2))

summary_df <- data.frame(
  true     = true_pars,
  mean_est = colMeans(est_mat),
  bias     = bias,
  rmse     = rmse,
  row.names = names(true_pars)
)
print(round(summary_df, 5))

# ── 4. Violin plots: estimate distribution per parameter ──────────────────
library(ggplot2)

est_long <- as.data.frame(est_mat)
est_long$rep <- seq_len(nrow(est_long))
est_long <- tidyr::pivot_longer(est_long, -rep,
                                names_to = "parameter", values_to = "estimate")
est_long$parameter <- factor(est_long$parameter, levels = names(true_pars))

ref_df <- data.frame(
  parameter = factor(names(true_pars), levels = names(true_pars)),
  true      = true_pars,
  lb        = lb,
  ub        = ub
)

ggplot(est_long, aes(x = "", y = estimate)) +
  geom_violin(fill = "steelblue", alpha = 0.4, colour = "navy") +
  geom_jitter(width = 0.1, alpha = 0.5, colour = "navy", size = 1) +
  geom_hline(data = ref_df, aes(yintercept = true),
             colour = "red", linetype = "dashed", linewidth = 0.8) +
  geom_hline(data = ref_df, aes(yintercept = lb),
             colour = "darkorange", linetype = "dotted", linewidth = 0.8) +
  geom_hline(data = ref_df, aes(yintercept = ub),
             colour = "forestgreen", linetype = "dotted", linewidth = 0.8) +
  facet_wrap(~parameter, scales = "free_y") +
  labs(x = NULL, y = "Estimate",
       title = "Parameter estimate distributions",
       caption = paste("Red dashed = true value",
                       "| Orange dotted = lower bound",
                       "| Green dotted = upper bound")) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# ── 5. Tree size vs absolute error per parameter ──────────────────────────
err_mat  <- abs(est_mat - rep(true_pars, each = nrow(est_mat)))
err_long <- as.data.frame(err_mat)
err_long$tips <- tips_ok
err_long <- tidyr::pivot_longer(err_long, -tips,
                                names_to = "parameter", values_to = "abs_error")
err_long$parameter <- factor(err_long$parameter, levels = names(true_pars))

ggplot(err_long, aes(x = tips, y = abs_error)) +
  geom_point(alpha = 0.5, colour = "steelblue") +
  geom_smooth(method = "loess", se = TRUE, colour = "red") +
  facet_wrap(~parameter, scales = "free_y") +
  labs(x = "Extant tips (tree size)", y = "|Error|",
       title = "Estimation error vs tree size") +
  theme_bw()
```

The output `summary_df` contains, for each parameter:

| Column | Meaning |
|--------|---------|
| `true` | True value used to simulate |
| `mean_est` | Mean estimate across successful replicates |
| `bias` | `mean_est - true` |
| `rmse` | Root mean squared error |

The tree-size plot shows whether larger trees (more extant tips) yield
more accurate estimates — expected if the IS estimator is consistent.

---

## Simulation study 2 — model selection power

Simulate 25 trees from each of the four pure models (CR, DD, PD, EP), fit all four
models to every tree via CEM, select the best by AIC, and evaluate how often the
correct model is recovered.

```r
library(emphasis)

set.seed(42)
n_per_model <- 25
max_t       <- 8
ctrl        <- list(max_iter = 20, num_particles = 80, num_trees = 1)

# True parameters and bounds for each model
configs <- list(
  CR = list(
    pars = c(0.6, 0.15),
    lb   = c(0.01,  0),
    ub   = c(2,     1)
  ),
  DD = list(
    pars = c(0.8, -0.02, 0.2, 0.0),
    lb   = c(0.01, -0.1,  0,  -0.1),
    ub   = c(2,     0.0,  1,   0.1)
  ),
  PD = list(
    pars = c(0.6,  0.005, 0.15, -0.003),
    lb   = c(0.01, -0.1,  0,    -0.1),
    ub   = c(2,     0.1,  1,     0.1)
  ),
  EP = list(
    pars = c(0.5, 0.05, 0.1, 0.01),
    lb   = c(0.01, -0.1, 0,  -0.1),
    ub   = c(2,     0.2, 1,   0.2)
  )
)
model_names <- names(configs)

# ── 1. Simulate trees ──────────────────────────────────────────────────────
all_trees   <- list()
true_models <- character()
for (m in model_names) {
  cat("Simulating", m, "trees...\n")
  count <- 0L
  while (count < n_per_model) {
    lnk <- if (m == "CR") "linear" else "exponential"
    sim <- simulate_tree(pars = configs[[m]]$pars, model = tolower(m), max_t = max_t,
                         link = lnk)
    if (sim$status == "done") {
      count <- count + 1L
      all_trees   <- c(all_trees,   list(sim))
      true_models <- c(true_models, m)
    }
  }
}

# ── 2. Fit all four models to every tree ───────────────────────────────────
n_trees    <- length(all_trees)
selected   <- character(n_trees)

for (i in seq_len(n_trees)) {
  cat(sprintf("Tree %d/%d (true: %s)\n", i, n_trees, true_models[i]))
  fits <- list()
  for (m in model_names) {
    fits[[m]] <- tryCatch(
      estimate_rates(all_trees[[i]], method = "cem", model = tolower(m),
                     link = if (m == "CR") "linear" else "exponential",
                     control = c(ctrl,
                                 list(lower_bound = configs[[m]]$lb,
                                      upper_bound = configs[[m]]$ub))),
      warning = function(w) NULL,
      error   = function(e) NULL
    )
  }
  # Select by lowest AIC among models that converged
  aics <- sapply(fits, function(f) if (!is.null(f)) f$AIC else NA_real_)
  if (any(is.finite(aics))) {
    selected[i] <- model_names[which.min(aics)]
  } else {
    selected[i] <- NA_character_
  }
}

# ── 3. Confusion matrix ────────────────────────────────────────────────────
ok      <- !is.na(selected)
conf    <- table(True = true_models[ok], Selected = selected[ok])
# reorder rows/cols to canonical order
conf <- conf[model_names, model_names]
cat("\nModel selection confusion matrix (rows=true, cols=selected):\n")
print(conf)

# Overall accuracy
accuracy <- sum(diag(conf)) / sum(conf)
cat(sprintf("Overall accuracy: %.1f%%\n", 100 * accuracy))

# Per-model recall (% correctly identified)
recall <- diag(conf) / rowSums(conf)
cat("\nPer-model recall:\n")
print(round(recall, 2))

# ── 4. Confusion matrix heatmap ───────────────────────────────────────────
conf_pct <- sweep(conf, 1, rowSums(conf), "/")  # row-normalise to recall
image(t(conf_pct[nrow(conf_pct):1, ]),
      axes = FALSE, col = grDevices::colorRampPalette(c("white", "steelblue"))(20),
      main = sprintf("Model selection confusion\n(row-normalised, accuracy=%.0f%%)",
                     100 * accuracy))
axis(1, at = seq(0, 1, length.out = 4), labels = model_names)
axis(2, at = seq(1, 0, length.out = 4), labels = model_names, las = 1)
mtext("Selected", side = 1, line = 3)
mtext("True", side = 2, line = 3)
for (i in seq_len(nrow(conf_pct)))
  for (j in seq_len(ncol(conf_pct)))
    text((j - 1) / 3, 1 - (i - 1) / 3,
         sprintf("%.0f%%", 100 * conf_pct[i, j]), cex = 0.9)

# ── 5. Bar chart: selection frequency per true model ─────────────────────
op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
cols <- c(CR = "grey70", DD = "steelblue", PD = "darkorange", EP = "forestgreen")
for (m in model_names) {
  idx   <- true_models[ok] == m
  votes <- table(factor(selected[ok][idx], levels = model_names))
  barplot(as.numeric(votes), names.arg = model_names, col = cols,
          main = sprintf("True model: %s", m),
          ylab = "# trees selected", ylim = c(0, n_per_model))
  abline(h = sum(idx) * 0.5, lty = 2, col = "red")
}
par(op)
```

The confusion matrix rows are true models, columns are selected models. A perfect
classifier has all mass on the diagonal. The bar charts show, for each true model,
how often each candidate was chosen — useful for diagnosing systematic confusions
(e.g. CR often mistaken for DD when the density-dependence signal is weak).

---

## GAM-based methods

Two GAM-based workflows complement the iterative EM/CE methods. Both
require the `mgcv` package.

### Use case 1: Survival conditioning

Estimate the probability that a forward simulation survives (doesn't go
extinct), as a function of parameters. Pass the trained GAM to
`estimate_rates(cond = ...)` to condition the likelihood on survival.

```r
library(emphasis)
set.seed(1)

# 1. Batch simulation: one tree per parameter row, no retries
pars_mat <- cbind(beta_0  = runif(1000, 0.3, 1.2),
                  gamma_0 = runif(1000, 0.05, 0.5))
sims <- simulate_tree(pars = pars_mat, max_t = 8, model = "cr", max_tries = 0)

# 2. Fit survival GAM (once — reuse across all trees)
gam_surv <- train_GAM(sims$simulations, pars_mat, model = "cr")

# 3. Predict survival probability at new parameter values
predict_survival(gam_surv, data.frame(beta_0 = 0.8, gamma_0 = 0.2))

# 4. Use in estimation — conditioned log-likelihood
fit <- estimate_rates(tree, method = "gam", model = "cr", cond = gam_surv,
                      control = list(lower_bound = c(0.3, 0.05),
                                     upper_bound = c(1.2, 0.5)))
fit$loglik  # conditioned: ell_IS - log(P_tree)
```

The `cond` parameter works with all three methods (`"mcem"`, `"cem"`, `"gam"`).

### Use case 2: GAM-based MLE

Estimate the log-likelihood surface over a parameter grid via IS, fit a
GAM to it, then optimize the GAM directly — a one-shot alternative to
iterative MCEM/CEM.

```r
set.seed(42)
sim <- simulate_tree(pars = c(0.5, 0.1), model = "cr", max_t = 8, useDDD = TRUE)

# Create parameter grid
pars_grid <- expand.grid(
  beta_0  = seq(0.2, 1.0, length.out = 10),
  gamma_0 = seq(0.01, 0.5, length.out = 10)
)

# Evaluate IS log-likelihood at each grid point
surface <- estimate_likelihood_surface(
  sim, as.matrix(pars_grid), model = "cr",
  sample_size = 100, verbose = TRUE
)

# Fit GAM and find MLE
gam_fit <- train_likelihood_GAM(surface)
mle     <- find_MLE(gam_fit,
                    lower_bound = c(0.2, 0.01),
                    upper_bound = c(1.0, 0.5))
mle$par   # MLE estimates
mle$value # log-likelihood at MLE
```

---

## Diagnostics — inference quality

Both `diagnose_cem()` and `diagnose_mcem()` produce four diagnostic plots
and return the underlying data. See the CEM and MCEM subsections above for
usage examples integrated with the estimation workflow.

**CEM diagnostic fields:**

| Field | Content |
|-------|---------|
| `convergence` | Per-iteration: `best_loglik`, `n_valid`, `rej_total` |
| `population_spread` | Per-iteration: median/q25/q75 of all particle `fhat` values |
| `IS_quality` | ESS, ESS fraction, mean/sd of IS log-weights |
| `final_pop` | Parameter values + `fhat` for each particle at the final iteration |
| `best_IS` | Raw IS data (`logf`, `log_q`, `lw`) at the best particle |

**MCEM diagnostic fields:**

| Field | Content |
|-------|---------|
| `convergence` | Per-iteration: `fhat`, `delta_max`, `num_trees`, `time` |
| `IS_quality` | ESS, mean/sd of log-weights, bootstrap variance |
| `final_IS` | Raw IS data (`logf`, `logg`, `lw`) from the final E-step |

---

## Authors

**Francisco Richter** (author & maintainer) — <richtf@usi.ch>

### Former contributors

Thijs Janzen — <t.janzen@rug.nl>
Hanno Hildenbrandt — <h.hildenbrandt@rug.nl>

## License

MIT. See [LICENSE](LICENSE) for details.

# emphasis — Evolutionary Modeling for Phylogenetic Inference

[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](./LICENSE)
![Lifecycle: Experimental](https://img.shields.io/badge/lifecycle-experimental-orange)

`emphasis` is an R package for studying species diversification and its ecological
drivers. It provides:

- **Simulation** — stochastic phylogenies under flexible diversification models
- **Inference** — parameter estimation via MCEM and Cross-Entropy Method (CEM)
- **Model selection** — AIC-based comparison across competing models
- **Diagnostics** — survival-probability estimation via GAM

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
dd <- simulate_tree(pars = c(0.8, -0.02, 0.2, -0.005),          model = "dd", max_t = 5)
pd <- simulate_tree(pars = c(0.6, 0.005, 0.15, -0.003),         model = "pd", max_t = 5)
ep <- simulate_tree(pars = c(0.5, 0.05, 0.1, 0.01),             model = "ep", max_t = 5)
np <- simulate_tree(pars = c(0.5, -0.01, 0.005, 0.15, -0.002, 0.001),
                    model = ~ N + PD, max_t = 5)

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
obs <- simulate_tree(pars = c(0.6, 0.005, 0.15, -0.003), model = "pd", max_t = 5)

# Single augmented tree
aug <- simulate_tree(tree = obs, pars = c(0.6, 0.005, 0.15, -0.003),
                     model = "pd", n_trees = 1L)
aug$tas    # augmented phylo (extant + stochastic extinct lineages)
aug$log_q  # log q(augmented tree | extant tree, pars)

# Many augmented trees
augN <- simulate_tree(tree = obs, pars = c(0.6, 0.005, 0.15, -0.003),
                      model = "pd", n_trees = 100L)
augN$log_q       # log q values (length 100)
augN$trees[[1]]  # first augmented phylo
```

---

## Inference

`estimate_rates()` estimates parameters from a phylogenetic tree. It accepts a
`simulate_tree()` result, a `phylo` object, or a numeric branching-time vector.

Two methods are available:

| Method | Description |
|--------|-------------|
| `"cem"` | Monte Carlo Cross-Entropy Method — broad search, no starting values needed |
| `"mcem"` | Monte Carlo EM — local refinement, converges to a (local) optimum |

```r
set.seed(42)
sim <- simulate_tree(pars = c(0.6, -0.05, 0.15, -0.03), model = "pd", max_t = 10)

fit <- estimate_rates(sim, method = "cem", model = "pd",
  control = list(
    lower_bound   = c(0.01, -0.5, 0, -0.5),
    upper_bound   = c(2,     0.5, 1,  0.5),
    verbose       = TRUE,
    max_iter      = 50,
    num_particles = 80,
    num_trees     = 10   # > 1 enables automatic bootstrap variance (B=200)
  ))

fit              # prints pars, loglik (MC se), AIC, convergence reason
fit$pars         # named estimates: beta_0, beta_P, gamma_0, gamma_P
fit$loglik
fit$AIC
fit$converged    # stopping reason: "annealing", "plateau", or "max_iter"
fit$loglik_var   # MC sampling variance of log-likelihood
```

### Diagnostics

Both methods have dedicated diagnostic functions:

```r
# CEM diagnostics
diag <- diagnose_cem(fit_cem, lower_bound = lb, upper_bound = ub)

# MCEM diagnostics
diag <- diagnose_mcem(fit_mcem, lower_bound = lb, upper_bound = ub)
```

In a simulation study where the true parameters are known, pass them via
`true_pars` to overlay red reference lines:

```r
true_pars <- c(beta_0 = 0.6, beta_P = -0.05, gamma_0 = 0.15, gamma_P = -0.03)
diag <- diagnose_cem(fit, lower_bound = lb, upper_bound = ub,
                     true_pars = true_pars)
diag <- diagnose_mcem(fit, lower_bound = lb, upper_bound = ub,
                      true_pars = true_pars)
```

### CEM stopping rules

The CEM stops at the first of:

1. **Annealing exhausted** — perturbation SD reaches zero; population has collapsed
2. **Plateau** — best fhat has not improved by `> tol` for `patience` consecutive iterations
3. **Hard cap** — `max_iter` reached

### CEM control parameters

All parameters are passed via `control = list(...)`. Inspect defaults with
`estimate_rates_control("cem")`.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `lower_bound` | — | **Required.** Lower parameter bounds vector |
| `upper_bound` | — | **Required.** Upper parameter bounds vector |
| `verbose` | FALSE | Print iteration summaries |
| `max_iter` | 20 | Hard iteration cap |
| `num_particles` | 50 | Particles (parameter vectors) per iteration. Alias: `num_points` |
| `num_trees` | 1 | Augmented trees per particle. When `> 1`, bootstrap variance of loglik is computed automatically (B=200). Recommend `>= 5`. Alias: `sample_size` |
| `shared_trees` | FALSE | Pool trees across all particles (`TRUE` = lower variance, slower) |
| `disc_prop` | 0.5 | Elite fraction of particles kept for resampling |
| `sd_vec` | NULL | Initial perturbation SDs; auto = `(upper - lower) / 4` |
| `maxN` | 10 | Max augmented-tree attempts before a particle is rejected |
| `max_missing` | 1e4 | Max missing lineages tolerated during augmentation |
| `max_lambda` | 500 | Max speciation rate during augmentation |
| `tol` | 1e-4 | Minimum fhat improvement to reset the plateau counter |
| `patience` | 5 | Consecutive plateau iterations before early stop |
| `bias_correct` | FALSE | Moment-based IS bias correction (experimental) |
| `num_threads` | 1 | Parallel threads for particle evaluation |

**Total augmented trees per iteration** = `num_particles × num_trees`.

**How `fhat` is computed per particle:**

```
fhat(theta) = log( mean_i[ exp(logf_i - log_q_i) ] )

  logf_i  = log p(obs, z_i | theta)   # joint likelihood of observed + augmented tree
  log_q_i = log q(z_i | obs, theta)   # proposal density of the augmented tree
  mean over i = 1 ... num_trees augmented trees drawn for this particle
```

This is the importance-sampling (IS) log-likelihood estimate. A higher `fhat` means
the parameter vector `theta` is better supported by the data. Each iteration the
best particle (`fhat*`) is recorded, and the population is resampled around it.

In Mode 1 (default, `shared_trees = FALSE`): each particle draws its own `num_trees`
trees. In Mode 2 (`shared_trees = TRUE`): all particles share a common pool of trees,
reducing variance but requiring more trees per iteration.

```r
# Higher-quality run: 100 particles, 5 trees each => 500 trees/iter
# num_trees >= 5 also enables automatic bootstrap variance of log-likelihood
estimate_rates(sim, method = "cem", model = "pd",
  control = list(
    lower_bound   = c(0.01, -0.5, 0, -0.5),
    upper_bound   = c(2,     0.5, 1,  0.5),
    verbose       = TRUE,
    max_iter      = 30,
    num_particles = 100,
    num_trees     = 5,
    patience      = 8
  ))
```

### MCEM control parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `lower_bound` | — | **Required.** Lower parameter bounds vector |
| `upper_bound` | — | **Required.** Upper parameter bounds vector |
| `verbose` | FALSE | Print progress messages |
| `num_trees` | 200 | Augmented trees per EM iteration. Alias: `sample_size`. Bootstrap variance of `loglik` computed automatically (B = 200) |
| `tol` | 0.1 | Convergence: SE of fhat across post-burnin iterations |
| `burnin` | 20 | Burn-in iterations before convergence is checked |
| `xtol` | 1e-3 | Relative tolerance for the M-step optimizer |
| `max_missing` | 1e4 | Max missing lineages during augmentation |
| `max_lambda` | 500 | Max speciation rate during augmentation |
| `num_threads` | 1 | Parallel threads |

### Empirical trees

```r
data(bird.orders, package = "ape")
extant <- prune_to_extant(bird.orders)

fit <- estimate_rates(extant, method = "cem", model = "dd",
  control = list(lower_bound = c(0.01, -0.5, 0, -0.5),
                 upper_bound = c(2,     0.5, 1,  0.5)))
fit$pars
```

---

## Model selection

`compare_models()` takes two or more fitted models and returns a table sorted by
AIC. The model with the lowest AIC is preferred; `delta_AIC` shows the gap.

```r
set.seed(42)
sim <- simulate_tree(pars = c(0.8, -0.02, 0.2, 0), model = "dd", max_t = 8)

fit_cr <- estimate_rates(sim, method = "cem", model = "cr",
  control = list(lower_bound = c(0,    0),              upper_bound = c(2,    1),
                 max_iter = 15, num_particles = 40))
fit_dd <- estimate_rates(sim, method = "cem", model = "dd",
  control = list(lower_bound = c(0.1, -0.1, 0, -0.01), upper_bound = c(2, 0.01, 0.5, 0.01),
                 max_iter = 15, num_particles = 40))
fit_pd <- estimate_rates(sim, method = "cem", model = "pd",
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
    sim <- simulate_tree(pars = true_pars, model = model, max_t = max_t)
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
    estimate_rates(trees[[i]], method = "cem", model = model, control = ctrl),
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
    sim <- simulate_tree(pars = configs[[m]]$pars, model = tolower(m), max_t = max_t)
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

## Diagnostics — inference quality

### CEM diagnostics

`diagnose_cem()` takes a fitted `"cem"` model and returns convergence
statistics, IS quality metrics, and the full particle cloud at the final
iteration. It produces diagnostic plots by default.

```r
diag <- diagnose_cem(fit, lower_bound = lb, upper_bound = ub)

diag$convergence       # per-iteration: best_loglik, n_valid, rej_lambda, rej_overruns, rej_total
diag$population_spread # per-iteration: median/q25/q75/range of all fhat values
diag$IS_quality        # ESS, ESS fraction, mean/sd of IS log-weights
diag$final_pop         # data frame: parameter columns + fhat for each particle
diag$best_IS           # raw IS data: logf, log_q, lw, trees at best pars
```

| Plot | What it shows |
|------|--------------|
| Log-likelihood | Best `fhat` per iteration with stop reason |
| Parameter traces | Best particle's value per iteration; orange = lower bound, purple = upper bound, red = true value |
| Population spread | Median with q25/q75 band of all valid `fhat` values per iteration |
| IS weight distribution | Histogram of `logf - log_q` at the best particle; ESS in title |
| Rejection rates | Per-iteration count of rejected trees (shown only when rejections occurred) |

### MCEM diagnostics

`diagnose_mcem()` provides analogous diagnostics for MCEM fits:
log-likelihood trace, parameter traces, SE convergence trace, and IS
weight histogram at the final E-step.

```r
diag <- diagnose_mcem(fit, lower_bound = lb, upper_bound = ub)

diag$convergence  # per-iteration: fhat, num_trees, time
diag$se_trace     # running SE of fhat (convergence criterion)
diag$IS_quality   # ESS, mean/sd of log-weights, bootstrap variance
diag$final_IS     # raw IS data: logf, logg, lw from last E-step
```

| Plot | What it shows |
|------|--------------|
| Log-likelihood trace | `fhat` vs. EM iteration |
| Parameter traces | M-step estimates per iteration; orange/purple = bounds, red = true value |
| SE trace | Running SE(fhat) — convergence declared when SE < tol |
| IS weight distribution | Histogram of IS log-weights at the final E-step; ESS in title |

---

## Diagnostics — survival probability via GAM

`train_GAM()` estimates P(clade survives to present | parameters) by fitting a
binomial GAM to batch forward simulations. Run with `max_tries = 0` to keep
extinct realisations alongside successful ones.

```r
set.seed(1)
pars_mat <- cbind(beta_0  = runif(300, 0.3, 1.5),
                  beta_N  = -0.02,
                  gamma_0 = runif(300, 0.05, 0.6),
                  gamma_N =  0.0)
sims <- simulate_tree(pars = pars_mat, max_t = 10,
                      model = "dd", max_tries = 0, useDDD = FALSE)

gam_fit <- train_GAM(sims, pars_mat, model = "dd")

newpars <- data.frame(beta_0 = 0.8, beta_N = -0.02,
                      gamma_0 = 0.2, gamma_N = 0.0)
predict_survival(gam_fit, newpars)  # estimated survival probability
```

---

## Authors

**Francisco Richter** (author & maintainer) — <richtf@usi.ch>
Thijs Janzen — <t.janzen@rug.nl>
Hanno Hildenbrandt — <h.hildenbrandt@rug.nl>

## License

GPL-3. See [LICENSE](LICENSE) for details.

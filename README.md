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

The recommended workflow runs CEM first to find a good region, then MCEM to
refine:

```r
set.seed(42)
sim <- simulate_tree(pars = c(0.6, 0.005, 0.15, -0.003), model = "pd", max_t = 8)

lb <- c(0.01, -1, 0, -1)
ub <- c(2,     1, 1,  1)

# Stage 1: broad search with bootstrap UQ
fit_cem <- estimate_rates(sim, method = "cem", model = "pd",
  lower_bound = lb, upper_bound = ub,
  control = list(max_iter = 20, num_points = 50, n_boot = 200))

fit_cem$converged    # stopping reason: "annealing", "plateau", or "max_iter"
fit_cem$loglik_var   # MC sampling variance of log-likelihood

# Stage 2: refine
fit_em <- estimate_rates(sim, method = "mcem", model = "pd",
  lower_bound = lb, upper_bound = ub,
  init_pars   = fit_cem$pars,
  control     = list(sample_size = 200, tol = 0.05))

fit_em          # prints pars, loglik (MC se), AIC
fit_em$pars     # named estimates: beta_0, beta_P, gamma_0, gamma_P
fit_em$loglik
fit_em$AIC
```

### CEM stopping rules

The CEM stops at the first of:

1. **Annealing exhausted** — perturbation SD reaches zero; population has collapsed
2. **Plateau** — best fhat has not improved by `> tol` for `patience` consecutive iterations
3. **Hard cap** — `max_iter` reached

### Control parameters

```r
estimate_rates_control("mcem")  # sample_size=200, tol=0.1, burnin=20, ...
estimate_rates_control("cem")   # max_iter=20, num_points=50, tol=1e-4, patience=5, ...
```

Override any default by passing a named list to `control`:

```r
estimate_rates(sim, method = "cem", model = "pd",
  lower_bound = lb, upper_bound = ub,
  control = list(max_iter = 30, num_points = 100, patience = 10, n_boot = 200))
```

### Empirical trees

```r
data(bird.orders, package = "ape")
extant <- prune_to_extant(bird.orders)

fit <- estimate_rates(extant, method = "cem", model = "dd",
  lower_bound = c(0.01, -0.5, 0, -0.5),
  upper_bound = c(2,     0.5, 1,  0.5))
fit$pars
```

---

## Model selection

`compare_models()` takes two or more fitted models and returns a table sorted by
AIC. The model with the lowest AIC is preferred; `delta_AIC` shows the gap.

```r
set.seed(42)
sim <- simulate_tree(pars = c(0.8, -0.02, 0.2, 0), model = "dd", max_t = 8)

ctrl <- list(max_iter = 15, num_points = 40)

fit_cr <- estimate_rates(sim, method = "cem", model = "cr",
  lower_bound = c(0,    0),              upper_bound = c(2,    1),    control = ctrl)
fit_dd <- estimate_rates(sim, method = "cem", model = "dd",
  lower_bound = c(0.1, -0.1, 0, -0.01), upper_bound = c(2, 0.01, 0.5, 0.01), control = ctrl)
fit_pd <- estimate_rates(sim, method = "cem", model = "pd",
  lower_bound = c(0.01, -1,  0, -1),    upper_bound = c(2,   1,   1,  1),    control = ctrl)

compare_models(CR = fit_cr, DD = fit_dd, PD = fit_pd)
#   model n_pars  loglik    AIC delta_AIC
# 1    DD      4  -12.3   32.6      0.00
# 2    PD      4  -15.1   38.2      5.60
# 3    CR      2  -18.9   41.8      9.20
```

---

## Diagnostics — CEM inference quality

`diagnose_cem()` takes a fitted `"cem"` model and returns convergence
statistics, IS quality metrics, and the full particle cloud at the final
iteration. Optionally produces three diagnostic plots.

```r
fit <- estimate_rates(sim, method = "cem", model = "pd",
  lower_bound = lb, upper_bound = ub,
  control = list(max_iter = 20, n_boot = 200))

diag <- diagnose_cem(fit)          # produces plots by default

diag$convergence       # per-iteration: best_loglik, n_valid, rejections
diag$population_spread # per-iteration: median/IQR/range of all fhat values
diag$IS_quality        # ESS, ESS fraction, mean/sd of IS log-weights
diag$final_pop         # data frame: parameter columns + fhat for each particle
diag$best_IS           # raw IS data: logf, log_q, lw, trees at best pars

# IS effective sample size (fraction of 1 = perfect; close to 0 = one dominant tree)
diag$IS_quality$ESS_fraction

# Access raw augmented trees from the best particle's final evaluation
diag$best_IS$trees[[1]]   # first augmented tree (C++ data frame format)
diag$best_IS$lw            # IS log-weights for all trees
```

### Plots produced

| Plot | What it shows |
|------|--------------|
| Convergence | Best `fhat` per iteration with stop reason |
| Population spread | Median ± IQR of all valid `fhat` values per iteration |
| IS weight distribution | Histogram of `logf - log_q` at the best particle; red dashed line = mean; title shows ESS |

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

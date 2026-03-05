# emphasis — Evolutionary Modeling for Phylogenetic Inference

[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](./LICENSE)
![Lifecycle: Experimental](https://img.shields.io/badge/lifecycle-experimental-orange)

## Overview

`emphasis` is an R package for studying species diversification and its
ecological drivers. It combines fast C++ back-ends with R-friendly interfaces
for:

- **Simulation** — stochastic phylogenies under general diversification models
- **Inference** — parameter estimation via MCEM and differential evolution
- **Model selection** — AIC-based comparison across competing models
- **Diagnostics** — survival-probability estimation via GAM

## Installation

```r
devtools::install_github("franciscorichter/emphasis")
```

---

## Quick start

```r
library(emphasis)

# 1. Simulate a tree under diversity dependence
set.seed(42)
sim <- simulate_tree(pars = c(0.8, -0.02, 0.2, 0), max_t = 8, model = "dd")

# 2. Estimate parameters (two-stage workflow)
lb <- c(0.1, -0.1, 0, -0.01)
ub <- c(2,    0.01, 0.5, 0.01)

fit_de <- estimate_rates(sim, method = "mcde", model = "dd",
  lower_bound = lb, upper_bound = ub,
  control = list(num_iterations = 10, num_points = 30))

fit_em <- estimate_rates(sim, method = "mcem", model = "dd",
  lower_bound = lb, upper_bound = ub,
  init_pars = fit_de$pars,
  control = list(sample_size = 200, tol = 0.1, burnin = 10))

fit_em
# emphasis fit (mcem, model = N)
#
# Parameters:
#   beta_0   beta_N  gamma_0  gamma_N
#   0.7890  -0.0195   0.1980   0.0005
#
# Log-likelihood: -12.34
# AIC:            32.68
# n_pars:         4

# 3. Compare models
fit_cr <- estimate_rates(sim, method = "mcem", model = "cr",
  lower_bound = c(0, 0), upper_bound = c(2, 1))

compare_models(CR = fit_cr, DD = fit_em)
#   model n_pars   loglik    AIC delta_AIC
# 1    N      4  -12.34  32.68      0.00
# 2   CR      2  -18.90  41.80      9.12
```

---

## Rate model

All models are special cases of a general per-lineage rate formulation:

```
lambda(s,t) = f(beta_0  + beta_N * N(t) + beta_P * P(t) + beta_E * E(s,t))
mu(s,t)     = f(gamma_0 + gamma_N * N(t) + gamma_P * P(t) + gamma_E * E(s,t))
```

| Covariate | Meaning |
|-----------|---------|
| `N(t)` | Current species richness (diversity) |
| `P(t)` | Total pendant branch-length of alive lineages (phylogenetic diversity) |
| `E(s,t)` | Pendant edge of lineage *s* — time since its last divergence (lineage-specific) |

Two **link functions** are available via the `link` argument:

| Link | `f(eta)` | When to use |
|------|----------|-------------|
| `"linear"` (default) | `max(0, eta)` | Rates linear in covariates; clipped to non-negative |
| `"exponential"` | `exp(eta)` | Always positive; log-linear relationship |

## Model specification

The `model` argument selects active covariates. Three equivalent formats:

**R formula** (recommended for mixed models):

```r
model = ~ 1           # constant rate (intercept only)
model = ~ N           # diversity dependence
model = ~ PD          # phylogenetic diversity dependence
model = ~ EP          # evolutionary pendant
model = ~ N + PD      # mixed: diversity + phylogenetic diversity
model = ~ N + PD + EP # full model
```

**String shortcut:**

| String | Formula | Active covariates |
|--------|---------|-------------------|
| `"cr"` | `~ 1` | none — constant rate |
| `"dd"` | `~ N` | N only — diversity dependence |
| `"pd"` | `~ PD` | P only — phylogenetic diversity |
| `"ep"` | `~ EP` | E only — evolutionary pendant |

**Binary vector:** `c(use_N, use_P, use_E)`, e.g. `c(1, 1, 0)` for N + P.

## Parameter vector

Parameters follow **lambda first, then mu**, each in N–P–E order. Only active
covariates are included; length = `2 + 2 * sum(model)`:

| Model | Length | Layout |
|-------|--------|--------|
| `"cr"` | 2 | `c(beta_0, gamma_0)` |
| `"dd"` | 4 | `c(beta_0, beta_N, gamma_0, gamma_N)` |
| `"pd"` | 4 | `c(beta_0, beta_P, gamma_0, gamma_P)` |
| `"ep"` | 4 | `c(beta_0, beta_E, gamma_0, gamma_E)` |
| `c(1,1,0)` | 6 | `c(beta_0, beta_N, beta_P, gamma_0, gamma_N, gamma_P)` |
| `c(1,1,1)` | 8 | full 8-parameter model |

---

## Simulation

`simulate_tree()` is the single entry point — both forward simulation (generating new
trees) and conditional simulation (augmenting observed trees with extinct lineages).

### Forward simulation

```r
set.seed(123)

# Constant rate
cr <- simulate_tree(pars = c(0.5, 0.1), max_t = 5, model = "cr")
cr$status   # "done" / "extinct" / "too_large"

# Diversity dependence
dd <- simulate_tree(pars = c(0.8, -0.02, 0.2, -0.005), max_t = 5, model = "dd")

# Phylogenetic diversity dependence
pd <- simulate_tree(pars = c(0.6, 0.005, 0.15, -0.003), max_t = 5, model = "pd")

# Evolutionary pendant
ep <- simulate_tree(pars = c(0.5, 0.05, 0.1, 0.01), max_t = 5, model = "ep")

# Mixed N + PD model via formula
np <- simulate_tree(pars = c(0.5, -0.01, 0.005, 0.15, -0.002, 0.001),
                    max_t = 5, model = ~ N + PD)

# Exponential link
cr_exp <- simulate_tree(pars = c(-0.7, -2.3), max_t = 5,
                         model = "cr", link = "exponential")
```

Each call returns:

| Element | Content |
|---------|---------|
| `tes`   | Extant-only phylogeny (`phylo`) |
| `tas`   | Full phylogeny with extinct lineages (`phylo`) |
| `L`     | L-table matrix (DDD format) |
| `brts`  | Branching times of the extant tree |
| `status` | `"done"`, `"extinct"`, or `"too_large"` |
| `model` | Resolved binary vector `c(use_N, use_P, use_E)` |
| `pars`  | Parameter vector as supplied |

### Batch simulation

When `pars` is a matrix, each row is a separate simulation:

```r
pars_mat <- cbind(beta_0  = runif(100, 0.2, 1.0),
                  gamma_0 = runif(100, 0.0, 0.3))
sims <- simulate_tree(pars = pars_mat, max_t = 5, model = "cr")
length(sims)          # 100
sims[[1]]$status
```

### Conditional simulation (augmentation)

Pass `tree =` to augment an observed extant tree with stochastically drawn
extinct lineages. Returns importance-sampling statistics needed for inference.

```r
obs <- simulate_tree(pars = c(0.6, 0.005, 0.15, -0.003), max_t = 5, model = "pd")

# Single augmented tree
aug1 <- simulate_tree(tree = obs, pars = c(0.6, 0.005, 0.15, -0.003),
                      model = "pd", n_trees = 1L)
aug1$tas      # augmented phylo (extant + stochastic extinct lineages)
aug1$logf     # log p(augmented tree | pars)
aug1$logg     # log q(augmented tree | extant tree, pars)
aug1$log_w    # logf - logg  (importance weight)
aug1$fhat     # MC estimate of log p(extant tree | pars)

# Many augmented trees (for E-step)
augN <- simulate_tree(tree = obs, pars = c(0.6, 0.005, 0.15, -0.003),
                      model = "pd", n_trees = 200L)
augN$fhat             # MC log-likelihood (scalar)
augN$log_w            # importance weights, length 200
augN$trees[[1]]$tas   # first augmented phylo
```

---

## Inference

`estimate_rates()` is the unified entry point for parameter estimation.

### Methods

| Method | Algorithm | When to use |
|--------|-----------|-------------|
| `"mcde"` | Monte Carlo Differential Evolution | Broad exploration; find good starting values |
| `"mcem"` | Monte Carlo Expectation-Maximisation | Refine from a starting point; converges to local optimum |

### Recommended workflow: MCDE then MCEM

```r
set.seed(42)
sim <- simulate_tree(pars = c(0.6, 0.005, 0.15, -0.003), max_t = 8, model = "pd")

lb <- c(0.01, -1, 0, -1)
ub <- c(2,     1, 1,  1)

# Stage 1 — MCDE: broad exploration
fit_de <- estimate_rates(sim, method = "mcde", model = "pd",
  lower_bound = lb, upper_bound = ub,
  control = list(num_iterations = 15, num_points = 50))

# Stage 2 — MCEM: refine from MCDE estimate
fit_em <- estimate_rates(sim, method = "mcem", model = "pd",
  lower_bound = lb, upper_bound = ub,
  init_pars = fit_de$pars,
  control = list(sample_size = 200, tol = 0.05))

fit_em$pars     # named parameter estimates
fit_em$loglik   # final log-likelihood
fit_em$AIC      # Akaike Information Criterion
```

### Supported models

All models work with both simulation and inference:

```r
# Constant rate
fit_cr <- estimate_rates(tree, model = "cr",
  lower_bound = c(0, 0), upper_bound = c(2, 1))

# Diversity dependence
fit_dd <- estimate_rates(tree, model = "dd",
  lower_bound = c(0.1, -0.1, 0, -0.01),
  upper_bound = c(2, 0.01, 0.5, 0.01))

# Phylogenetic diversity dependence
fit_pd <- estimate_rates(tree, model = "pd",
  lower_bound = c(0.01, -1, 0, -1),
  upper_bound = c(2, 1, 1, 1))

# Evolutionary pendant (linear link only)
fit_ep <- estimate_rates(tree, model = "ep",
  lower_bound = c(0.1, -0.5, 0, -0.5),
  upper_bound = c(2, 0.5, 0.5, 0.5))

# Mixed models via formula
fit_np <- estimate_rates(tree, model = ~ N + PD,
  lower_bound = c(0.01, -0.5, -0.5, 0, -0.5, -0.5),
  upper_bound = c(2, 0.5, 0.5, 1, 0.5, 0.5))

# Exponential link
fit_exp <- estimate_rates(tree, model = "cr", link = "exponential",
  lower_bound = c(-5, -5), upper_bound = c(2, 2))
```

### Empirical trees

```r
data(bird.orders, package = "ape")
extant <- prune_to_extant(bird.orders)
fit <- estimate_rates(extant, method = "mcde", model = "dd",
  lower_bound = c(0.01, -0.5, 0, -0.5),
  upper_bound = c(2, 0.5, 0.5, 0.5))
fit$pars
```

### Control parameters

Use `estimate_rates_control()` to inspect defaults:

```r
estimate_rates_control("mcem")
# $max_missing  1e4
# $max_lambda   500
# $num_threads  1
# $sampling     "dynamic_fresh"
# $sample_size  200
# $xtol         1e-3
# $tol          0.1
# $burnin       20

estimate_rates_control("mcde")
# $max_missing  1e4
# $max_lambda   500
# $num_threads  1
# $num_iterations  20
# $num_points      50
# $disc_prop       0.5
```

---

## Model selection

`compare_models()` ranks fitted models by AIC. The model with the lowest AIC
is preferred; `delta_AIC` shows the gap from the best model.

```r
set.seed(42)
sim <- simulate_tree(pars = c(0.8, -0.02, 0.2, 0), max_t = 8, model = "dd")
lb_cr <- c(0, 0);        ub_cr <- c(2, 1)
lb_dd <- c(0.1, -0.1, 0, -0.01); ub_dd <- c(2, 0.01, 0.5, 0.01)
lb_pd <- c(0.01, -1, 0, -1);     ub_pd <- c(2, 1, 1, 1)

fit_cr <- estimate_rates(sim, model = "cr",
  lower_bound = lb_cr, upper_bound = ub_cr)
fit_dd <- estimate_rates(sim, model = "dd",
  lower_bound = lb_dd, upper_bound = ub_dd)
fit_pd <- estimate_rates(sim, model = "pd",
  lower_bound = lb_pd, upper_bound = ub_pd)

compare_models(CR = fit_cr, DD = fit_dd, PD = fit_pd)
#   model n_pars  loglik    AIC delta_AIC
# 1    DD      4  -12.3   32.6      0.00
# 2    PD      4  -15.1   38.2      5.60
# 3    CR      2  -18.9   41.8      9.20
```

---

## Monte Carlo log-likelihood

`mc_loglik()` augments an observed extant tree and returns a Monte Carlo
estimate of the log-likelihood via importance sampling. This is the core
engine used by `estimate_rates()` internally.

```r
set.seed(1)
tr <- simulate_tree(pars = c(0.5, 0.1), max_t = 5, model = "cr")

ll <- mc_loglik(
  brts        = tr$brts,
  pars        = c(0.5, 0, 0, 0, 0.1, 0, 0, 0),   # full 8-element layout
  sample_size = 25,
  maxN        = 250,
  max_missing = 1e4,
  max_lambda  = 500,
  lower_bound = rep(-1e6, 8),
  upper_bound = rep( 1e6, 8),
  xtol_rel    = 1e-3,
  num_threads = 1,
  model       = c(0L, 0L, 0L)
)
ll$fhat      # Monte Carlo log-likelihood estimate
ll$weights   # log importance weights (logf - logg)
```

---

## Diagnostics — survival probability via GAM

`train_GAM()` estimates P(survival | parameters) by fitting a binomial GAM to
forward simulations run with `max_tries = 0` (extinct outcomes retained).

```r
set.seed(1)
pars_mat <- cbind(beta_0  = runif(500, 0.3, 1.5),
                  beta_N  = -0.02,
                  gamma_0 = runif(500, 0.05, 0.6),
                  gamma_N =  0.0)
sims <- simulate_tree(pars = pars_mat, max_t = 10,
                      model = "dd", max_tries = 0, useDDD = FALSE)

gam_fit <- train_GAM(sims, pars_mat, model = "dd")
summary(gam_fit)

newpars <- data.frame(beta_0 = 0.8, beta_N = -0.02,
                      gamma_0 = 0.2, gamma_N = 0.0)
predict_survival(gam_fit, newpars)
```

---

## EP model notes

The evolutionary pendant (EP) model makes rates lineage-specific: each
lineage's rate depends on how long it has been since it last speciated.
In inference, a **mean-field approximation** is used: speciation event terms
and the integral use `E_avg = P/N` (average pendant length per lineage),
while extinction event terms use the exact `E_s` for each augmented lineage.

EP inference currently supports the **linear link only**. The exponential link
with EP is guarded (the integral has no closed-form solution).

---

## Authors & Contact

**Author & Maintainer**
Francisco Richter — <richtf@usi.ch>

**Contributors**
Thijs Janzen — <t.janzen@rug.nl>
Hanno Hildenbrandt — <h.hildenbrandt@rug.nl>

## License

GPL-3. See [LICENSE](LICENSE) for details.

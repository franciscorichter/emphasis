# emphasis — Evolutionary Modeling for Phylogenetic Inference

[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](./LICENSE)
![Lifecycle: Experimental](https://img.shields.io/badge/lifecycle-experimental-orange)

## Overview

`emphasis` is an R package for studying species diversification and its ecological drivers. It combines fast C++ back-ends with R-friendly interfaces for:

- stochastic simulation of phylogenies under general diversification models,
- Monte Carlo log-likelihood estimation via importance sampling with augmented trees,
- Monte Carlo expectation-maximisation (MCEM) parameter estimation, and
- differential-evolution (DE) searches for good starting values.

The package is organised around three modules: **Simulation**, **Inference**, and **Diagnostics**.

## Installation

```r
# Development version from GitHub
devtools::install_github("franciscorichter/emphasis")
```

---

## Simulation module

`simulate_tree()` is the single entry point for all tree simulation — both forward (generating new trees) and conditional (augmenting an observed extant tree with extinct lineages).

### Rate model

All models are special cases of a single general per-lineage rate formulation:

```
lambda(s,t) = max(0, beta_0  + beta_N * N(t) + beta_P * P(t) + beta_E * E(s,t))
mu(s,t)     = max(0, gamma_0 + gamma_N * N(t) + gamma_P * P(t) + gamma_E * E(s,t))
```

| Covariate | Meaning |
|-----------|---------|
| `N(t)` | current species richness |
| `P(t)` | total pendant branch-length of alive lineages (phylogenetic diversity) |
| `E(s,t)` | pendant edge of lineage *s* (time since its last divergence) — makes rates **lineage-specific** |

### Model specification

The `model` argument selects active covariates via a 3-element binary vector `c(use_N, use_P, use_E)`. Named shortcuts are also accepted:

| String | Binary vector | Active covariates |
|--------|--------------|-------------------|
| `"cr"` | `c(0, 0, 0)` | none — constant rate |
| `"dd"` | `c(1, 0, 0)` | N only — diversity dependence |
| `"pd"` | `c(0, 1, 0)` | P only — phylogenetic diversity dependence |
| `"ep"` | `c(0, 0, 1)` | E only — evolutionary pendant |

Mixed models such as `c(1, 1, 0)` (N + P) or `c(1, 1, 1)` (full) are fully supported.

### Parameter vector

`pars` always follows **lambda parameters first, then mu parameters**, each in N–P–E order. Only active covariates are included; length = `2 + 2 * sum(model)`:

```
pars = c(beta_0,  [beta_N], [beta_P], [beta_E],
         gamma_0, [gamma_N], [gamma_P], [gamma_E])
```

| Model | Length | Layout |
|-------|--------|--------|
| `"cr"` | 2 | `c(beta_0, gamma_0)` |
| `"dd"` | 4 | `c(beta_0, beta_N, gamma_0, gamma_N)` |
| `"pd"` | 4 | `c(beta_0, beta_P, gamma_0, gamma_P)` |
| `"ep"` | 4 | `c(beta_0, beta_E, gamma_0, gamma_E)` |
| `c(1,1,0)` | 6 | `c(beta_0, beta_N, beta_P, gamma_0, gamma_N, gamma_P)` |
| `c(1,1,1)` | 8 | `c(beta_0, beta_N, beta_P, beta_E, gamma_0, gamma_N, gamma_P, gamma_E)` |

### Single-tree forward simulation

```r
library(emphasis)
set.seed(123)

# CR — constant rate
cr <- simulate_tree(pars = c(0.5, 0.1), max_t = 5, model = "cr")
cr$status   # "done" / "extinct" / "too_large"

# DD — diversity dependence
dd <- simulate_tree(pars = c(0.8, -0.02, 0.2, -0.005), max_t = 5, model = "dd")

# PD — phylogenetic diversity dependence
pd <- simulate_tree(pars = c(0.6, 0.005, 0.15, -0.003), max_t = 5, model = "pd")

# Mixed N+P model  c(1,1,0), 6 parameters
np <- simulate_tree(pars = c(0.5, -0.01, 0.005, 0.15, -0.002, 0.001),
                    max_t = 5, model = c(1, 1, 0))
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

### Batch forward simulation

When `pars` is a matrix, each row is a separate simulation. Returns a list of
per-row outputs with the same structure as above:

```r
set.seed(42)
pars_mat <- cbind(
  beta_0  = runif(100, 0.2, 1.0),
  gamma_0 = runif(100, 0.0, 0.3)
)
sims <- simulate_tree(pars = pars_mat, max_t = 5, model = "cr")
length(sims)          # 100
sims[[1]]$brts
sims[[1]]$status
```

### Conditional simulation — augmentation with importance-sampling statistics

Pass `tree =` to augment an observed extant tree with stochastically drawn
extinct lineages. The output always includes importance-sampling statistics
(`logf`, `logg`, `log_w`, `fhat`) needed for downstream inference.

```r
set.seed(42)
obs <- simulate_tree(pars = c(0.6, 0.005, 0.15, -0.003),
                     max_t = 5, model = "pd")

# Single augmented tree
aug1 <- simulate_tree(tree = obs, pars = c(0.6, 0.005, 0.15, -0.003),
                      model = "pd", n_trees = 1L)
aug1$tas      # augmented phylo (extant + stochastic extinct lineages)
aug1$logf     # log p(augmented tree | pars)
aug1$logg     # log q(augmented tree | extant tree, pars)
aug1$log_w    # logf - logg  (importance weight)
aug1$fhat     # MC estimate of log p(extant tree | pars)
```

For the E-step, draw many augmented trees at once with `n_trees > 1`:

```r
augN <- simulate_tree(tree = obs, pars = c(0.6, 0.005, 0.15, -0.003),
                      model = "pd", n_trees = 200L)
augN$fhat             # MC log-likelihood (scalar)
augN$log_w            # importance weights, length 200
augN$logf             # per-tree log-likelihoods, length 200
augN$trees[[1]]$tas   # first augmented phylo
```

Batch augmentation across many parameter sets (one augmented tree per row):

```r
pars_mat <- cbind(beta_0  = runif(50, 0.3, 0.9),
                  beta_P  = runif(50, -0.01, 0.01),
                  gamma_0 = runif(50, 0.05, 0.3),
                  gamma_P = runif(50, -0.005, 0.005))
results <- simulate_tree(tree = obs, pars = pars_mat, model = "pd")
sapply(results, `[[`, "fhat")   # fhat for each parameter set
```

#### Output fields — conditional simulation

| `n_trees` | Structure |
|-----------|-----------|
| `1` | Same flat list as forward sim, plus `logf`, `logg`, `log_w`, `fhat` |
| `> 1` | `trees` (list of `n_trees` results), `logf`/`logg`/`log_w` (vectors), `fhat` (scalar), `tes`, `brts`, `model`, `pars` |

### Architecture

```
simulate_tree()                            R/simulate.R
  ├─ [pars is matrix]   lapply over rows
  ├─ [forward]          .expand_pars()
  │                     simulate_div_tree_cpp()    src/div_tree.cpp
  │                       └─ general_div::simulate_tree_ltable()
  │                            Gillespie           inst/include/general_tree.hpp
  │                              → L-table → DDD::L2phylo()
  └─ [conditional]      .sim_tree_conditional()
                          .augment_tree_internal()
                            └─ mc_loglik()         src/rcpp_mce.cpp
                                 logf, logg, log_w, fhat returned
```

---

## Monte Carlo log-likelihood

`mc_loglik()` augments an observed extant tree with stochastically drawn extinct lineages and returns a Monte Carlo estimate of the log-likelihood via importance sampling.

```r
set.seed(1)
tr <- simulate_tree(pars = c(0.5, 0.1), max_t = 5, model = "cr")

ll <- mc_loglik(
  brts        = tr$brts,
  pars        = c(0.1, 0.5, -0.02, 0.01),
  sample_size = 25,
  maxN        = 250,
  max_missing = 1e4,
  max_lambda  = 500,
  lower_bound = rep(-1e6, 4),
  upper_bound = rep( 1e6, 4),
  xtol_rel    = 1e-3,
  num_threads = 1
)
ll$fhat      # Monte Carlo log-likelihood estimate
ll$weights   # log importance weights (logf - logg)
```

---

## Inference module

`estimate_rates()` is the unified entry point for parameter inference. It accepts a simulated or empirical tree and dispatches to the chosen method.

### Methods

| Method | Function | When to use |
|--------|----------|-------------|
| `"em"` | MCEM | Fine-grained refinement; converges to a local optimum |
| `"de"` | `emphasis_de()` | Broad exploration; good seed for EM |

### Two-stage workflow

```r
library(emphasis)
set.seed(42)
sim <- simulate_tree(pars = c(0.6, 0.005, 0.15, -0.003), max_t = 8, model = "pd")

lb <- c(0.01, -1, -1, -1)
ub <- c(2,     1,  1,  1)

# Stage 1 — DE: broad exploration
fit_de <- estimate_rates(sim$tes, method = "de",
  lower_bound = lb, upper_bound = ub,
  control = list(num_iterations = 15, num_points = 50))

# Stage 2 — EM: refine from DE estimate
fit_em <- estimate_rates(sim$tes, method = "em",
  lower_bound = lb, upper_bound = ub,
  init_pars = fit_de$pars,
  control = list(sample_size = 200, tol = 0.05))
fit_em$pars
fit_em$loglik
```

### Empirical trees

```r
data(bird.orders, package = "ape")
extant <- prune_to_extant(bird.orders)
fit <- estimate_rates(extant, method = "de",
  lower_bound = c(0.01, -0.5, -0.5, -0.5),
  upper_bound = c(2,     0.5,  0.5,  0.5))
fit$pars
```

---

## Diagnostics — survival probability via GAM

`train_GAM()` estimates the probability that a clade survives to the present,
P(survival | θ), by fitting a binomial GAM to forward simulations run with
`max_tries = 0` (so extinct runs are retained rather than retried).

```r
library(emphasis)
set.seed(1)

# Simulate a training grid — max_tries = 0 keeps extinct outcomes
pars_mat <- cbind(beta_0  = runif(500, 0.3, 1.5),
                  beta_N  = -0.02,            # fixed competition coefficient
                  gamma_0 = runif(500, 0.05, 0.6),
                  gamma_N =  0.0)

sims <- simulate_tree(pars = pars_mat, max_t = 10,
                      model = "dd", max_tries = 0, useDDD = FALSE)

# Fit binomial GAM: survived ~ s(beta_0) + s(beta_N) + s(gamma_0) + s(gamma_N)
gam_fit <- train_GAM(sims, pars_mat, model = "dd")
summary(gam_fit)

# Predict at new parameter values
newpars <- data.frame(beta_0 = 0.8, beta_N = -0.02,
                      gamma_0 = 0.2, gamma_N = 0.0)
predict_survival(gam_fit, newpars)
```

A full replication of the paper experiments (LDD and LPD models, heatmaps,
and empirical validation) is in
`inst/scripts/paper_replication.R`.

---

## Authors & Contact

**Author & Maintainer**
Francisco Richter — <richtf@usi.ch>

**Contributors**
Thijs Janzen — <t.janzen@rug.nl>
Hanno Hildenbrandt — <h.hildenbrandt@rug.nl>

## License

GPL-3. See [LICENSE](LICENSE) for details.

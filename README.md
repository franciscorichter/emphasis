# emphasis — Evolutionary Modeling for Phylogenetic Inference

[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](./LICENSE)
![Lifecycle: Experimental](https://img.shields.io/badge/lifecycle-experimental-orange)

## Overview

`emphasis` (Evolutionary Modeling Platform for Adaptive Speciation, Inference, and Simulation) is an R package for studying diversification of species and the ecological drivers of macro-evolution. It combines fast C++ back-ends with R-friendly interfaces for:

- stochastic simulation of phylogenies under general diversification models,
- Monte Carlo expectation-maximisation (MCEM) parameter estimation,
- differential-evolution (DE) searches for good starting values, and
- tooling to explore non-homogeneous Poisson/exponential processes.

The package is organised around a small set of verbs — **simulate**, **augment**, **estimate**, **diagnose** — each backed by an R helper in `R/` and a matching Rd page.

## Feature highlights

- **Simulation**: `simulate_tree()` covers all models from a single entry point, with a binary model vector selecting which covariates are active.
- **Augmentation**: `augmentPD()` and `AugmentMultiplePhyloPD()` produce augmented trees for downstream inference.
- **Estimation**: `emphasis_de()` / `emphasis_de_factorial()` and `estimate_rates()` implement robust inference workflows.
- **Utilities**: Non-homogeneous exponential samplers (`generateNonHomogeneousExp()`, `nhExpRand()`), gradient helpers, and GAM-based diagnostics.

Each high-level routine accepts plain R objects (e.g. `phylo`, numeric vectors) and delegates heavy lifting to `Rcpp` / `RcppParallel` code paths.

## Modular architecture

`emphasis` is structured around three modules that mirror a typical diversification analysis pipeline.

| Module | Purpose | Key ingredients |
|--------|---------|-----------------|
| **Module A — Simulation** | Generate phylogenies and synthetic corpora. Covers forward simulation and conditional augmentation of observed trees. | `simulate_tree()`, `augmentPD()`, `AugmentMultiplePhyloPD()` |
| **Module B — Inference** | Estimate diversification parameters with EM or differential evolution. | `estimate_rates()`, `estimate_rates_control()`, `emphasis_de()` |
| **Module C — Diagnostics** | Assess fit and communicate results through GAM diagnostics and summary dashboards. | `train_GAM()` |

## Installation

```r
# Install from CRAN (when available)
# install.packages("emphasis")

# Or install the development version from GitHub
# install.packages("devtools")
devtools::install_github("franciscorichter/emphasis")
```

## Simulation module

`simulate_tree()` is the single entry point for all forward and conditional tree simulation.

### Rate model

All diversification models in `emphasis` are special cases of a single general per-lineage rate formulation:

```
lambda(s,t) = max(0, beta_0  + beta_N * N(t) + beta_P * P(t) + beta_E * E(s,t))
mu(s,t)     = max(0, gamma_0 + gamma_N * N(t) + gamma_P * P(t) + gamma_E * E(s,t))
```

where:
- **N(t)** — current species richness
- **P(t)** — total pendant branch-length of alive lineages (phylogenetic diversity)
- **E(s,t)** — pendant edge of lineage *s* (time since its last divergence); makes rates **lineage-specific**

The default link function is identity (clamped to 0); both rates can depend on any combination of N, P, E.

### Model specification

The `model` argument controls which covariates are active via a 3-element binary vector `c(use_N, use_P, use_E)`. Named shortcuts are also accepted:

| String | Binary vector | Active covariates |
|--------|--------------|-------------------|
| `"cr"` | `c(0, 0, 0)` | none — constant rate |
| `"dd"` | `c(1, 0, 0)` | N only — diversity dependence |
| `"pd"` | `c(0, 1, 0)` | P only — phylogenetic diversity dependence |
| `"ep"` | `c(0, 0, 1)` | E only — evolutionary pendant |

Mixed models such as `c(1, 1, 0)` (N + P) or `c(1, 1, 1)` (full) are fully supported.

Models are nested: CR ⊂ any single-covariate model ⊂ mixed models ⊂ full model.

### Parameter vector

`pars` always follows the order **lambda parameters first, then mu parameters**, each in N–P–E order. Only parameters for active covariates are included:

```
pars = c(beta_0, [beta_N], [beta_P], [beta_E],
         gamma_0, [gamma_N], [gamma_P], [gamma_E])
```

Length is always `2 + 2 * sum(model)`:

| Model | Length | `pars` layout |
|-------|--------|----------------|
| `"cr"` — `c(0,0,0)` | 2 | `c(beta_0, gamma_0)` |
| `"dd"` — `c(1,0,0)` | 4 | `c(beta_0, beta_N, gamma_0, gamma_N)` |
| `"pd"` — `c(0,1,0)` | 4 | `c(beta_0, beta_P, gamma_0, gamma_P)` |
| `"ep"` — `c(0,0,1)` | 4 | `c(beta_0, beta_E, gamma_0, gamma_E)` |
| `c(1,1,0)` | 6 | `c(beta_0, beta_N, beta_P, gamma_0, gamma_N, gamma_P)` |
| `c(1,1,1)` | 8 | `c(beta_0, beta_N, beta_P, beta_E, gamma_0, gamma_N, gamma_P, gamma_E)` |

### Single-tree draws

```r
library(emphasis)
set.seed(123)

# CR — constant rate
cr <- simulate_tree(pars = c(0.5, 0.1), max_t = 5, model = "cr")
cr$status   # "done" / "extinct" / "too_large"

# DD — diversity dependence on N
dd <- simulate_tree(pars = c(0.8, -0.02, 0.2, -0.005), max_t = 5, model = "dd")

# PD — phylogenetic diversity dependence
pd <- simulate_tree(pars = c(0.6, 0.005, 0.15, -0.003), max_t = 5, model = "pd")

# EP — evolutionary pendant (per-lineage rates)
ep <- simulate_tree(pars = c(0.5, 0.01, 0.1, -0.003), max_t = 5, model = "ep")

# Mixed: N + P dependence  c(1,1,0), 6 parameters
np <- simulate_tree(pars = c(0.5, -0.01, 0.005, 0.15, -0.002, 0.001),
                    max_t = 5, model = c(1, 1, 0))

# Full model  c(1,1,1), 8 parameters
full <- simulate_tree(pars = c(0.5, -0.01, 0.005, 0.02, 0.15, -0.002, 0.001, -0.005),
                      max_t = 5, model = c(1, 1, 1))
```

The return value is always a list with:

| Element | Content |
|---------|---------|
| `tes`   | Extant-only phylogeny (`phylo`) |
| `tas`   | Full phylogeny including extinct lineages (`phylo` or `NULL`) |
| `L`     | L-table matrix in DDD format |
| `brts`  | Branching times of the extant tree |
| `status` | `"done"`, `"extinct"`, or `"too_large"` |
| `model` | Resolved binary vector `c(use_N, use_P, use_E)` |
| `pars`  | Parameter vector as supplied |

### Conditional simulation — augmenting an extant tree

`simulate_tree()` also serves as the conditional-simulation entry point. Pass `tree =` an observed extant tree and it stochastically fills in the missing extinct lineages, returning a `tas` phylo that includes both the original extant branches and the simulated extinct ones.

```r
set.seed(42)
sim <- simulate_tree(pars = c(0.6, 0.005, 0.15, -0.003), max_t = 5, model = "pd")

# Conditional simulation: augment an existing extant tree
aug <- simulate_tree(tree = sim, pars = c(0.6, 0.005, 0.15, -0.003), model = "pd")
aug$tes    # same extant-only phylo
aug$tas    # full phylo with stochastically drawn extinct lineages

# Also works with a plain phylo object as input
aug2 <- simulate_tree(tree = sim$tes, pars = c(0.6, 0.005, 0.15, -0.003), model = "pd")
```

| `tree` argument | Mode | Produces |
|-----------------|------|----------|
| `NULL` (default) | Forward simulation | new tree from scratch |
| `simulate_tree()` result, `phylo`, or brts vector | Conditional simulation | extant tree + stochastic extinct lineages |

## Inference module

`estimate_rates()` is the unified entry point for parameter inference. It accepts a simulated or empirical tree, prunes it to extant tips if needed, and dispatches to the chosen learning method.

### Methods

| Method | Function | When to use |
|--------|----------|-------------|
| `"em"` | MCEM | Fine-grained refinement; converges to a local optimum |
| `"de"` | `emphasis_de()` | Broad exploration; good seed for EM |

A robust two-stage workflow uses DE to find the basin and EM to refine:

```r
library(emphasis)
set.seed(42)
sim <- simulate_tree(pars = c(0.6, 0.005, 0.15, -0.003), max_t = 8, model = "pd")

lb <- c(0.01, -1, -1, -1)
ub <- c(2,     1,  1,  1)

# Stage 1 — DE: broad exploration
fit_de <- estimate_rates(sim, method = "de",
  lower_bound = lb, upper_bound = ub,
  control = list(num_iterations = 15, num_points = 50))
fit_de$pars

# Stage 2 — EM: refine from DE estimate
fit_em <- estimate_rates(sim, method = "em",
  lower_bound = lb, upper_bound = ub,
  init_pars = fit_de$pars,
  control = list(sample_size = 200, tol = 0.05))
fit_em$pars
fit_em$loglik
```

### Tuning

`estimate_rates_control()` shows all defaults and is the recommended way to build a `control` list:

```r
estimate_rates_control("em")   # view EM defaults
estimate_rates_control("de")   # view DE defaults

ctl <- estimate_rates_control("em")
ctl$sample_size <- 500
fit <- estimate_rates(sim, lower_bound = lb, upper_bound = ub,
                      init_pars = fit_de$pars, control = ctl)
```

### Empirical trees

Pass a `phylo` object directly. Use `prune_to_extant()` first if the tree contains extinct lineages:

```r
data(bird.orders, package = "ape")
extant <- prune_to_extant(bird.orders)
fit <- estimate_rates(extant, method = "de",
  lower_bound = c(0.01, -0.5, -0.5, -0.5),
  upper_bound = c(2,     0.5,  0.5,  0.5))
fit$pars
```

## Core modules

| Module | Purpose | Key functions |
|--------|---------|---------------|
| `simulate.R` | Unified simulator entry point | `simulate_tree()` |
| `inference.R` | Unified inference entry point | `estimate_rates()`, `estimate_rates_control()`, `prune_to_extant()` |
| `augment.R` | Batch augmentation of missing lineages | `augmentPD()`, `AugmentMultiplePhyloPD()` |
| `generate.R` | Non-homogeneous process samplers | `generateNonHomogeneousExp()`, `nhExpRand()`, `rate_t()`, `ExponentialRate()` |
| `de.R` | Differential-evolution back-end | `emphasis_de()`, `emphasis_de_factorial()` |
| `emphasis.R` | MCEM back-end | `mcEM_step()` |
| `gam.R` | Diagnostics over simulation sweeps | `train_GAM()` |

See `man/` for in-depth documentation on each function.

## Diagnostics & loop control

```r
results <- AugmentMultiplePhyloPD(
  phylo = bird.orders,
  n_trees = 20,
  mu_interval = c(0.05, 0.2),
  lambda_interval = c(0.1, 0.6),
  betaN_interval = c(-0.5, 0.1),
  betaP_interval = c(-0.2, 0.2)
)

gam_fit <- train_GAM(results)
summary(gam_fit)
```

## Documentation & vignettes

- Function docs live in `man/*.Rd` (generated via roxygen). Start with `?emphasis-package` for a package overview.
- Worked examples are available via `browseVignettes("emphasis")` after installing the package.
- The `Documentation/` directory contains PDF summaries of internal research experiments.

## Contributing

Pull requests and issues are welcome! A good contribution typically includes:

1. A reproducible example or failing test.
2. Updated documentation (`roxygen2` comments + README snippets).
3. Benchmark notes if the change affects performance-critical code.

## Authors & Contact

**Author & Maintainer**
Francisco Richter — <richtf@usi.ch>

**Contributors**
Thijs Janzen — <t.janzen@rug.nl>
Hanno Hildenbrandt — <h.hildenbrandt@rug.nl>

## License

GPL-3. See [LICENSE](LICENSE) for details.

---

For more information, see the package documentation and vignettes, or contact the maintainer.

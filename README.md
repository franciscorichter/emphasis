# emphasis: Evolutionary Modeling on Phylogenetic Applications

[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](./LICENSE)
![Lifecycle: Experimental](https://img.shields.io/badge/lifecycle-experimental-orange)

## Overview

`emphasis` is an R package for studying diversification of species and the ecological drivers of macro-evolution. It combines fast C++ back-ends with R-friendly interfaces for:

- stochastic simulation of phylogenies with density and PD dependence,
- Monte Carlo expectation-maximisation (MCEM) parameter estimation,
- differential-evolution (DE) searches for good starting values, and
- tooling to explore non-homogeneous Poisson/exponential processes.

The package is organised around a small set of verbs — **simulate**, **augment**, **estimate**, **diagnose** — each backed by an R helper in `R/` and a matching Rd page.

## Feature highlights

- **Simulation**: `sim_tree_pd_*`, `generatePhyloPD()` and friends quickly explore scenarios.
- **Augmentation**: `augmentPD()` and `AugmentMultiplePhyloPD()` produce augmented trees for downstream inference.
- **Estimation**: `emphasis_de()` / `emphasis_de_factorial()` and `mcEM_step()` implement robust inference workflows.
- **Utilities**: Non-homogeneous exponential samplers (`generateNonHomogeneousExp()`, `nhExpRand()`), gradient helpers, GAM-based diagnostics, and plotting utilities.

Each high-level routine accepts plain R objects (e.g. `phylo`, numeric vectors) and delegates heavy lifting to `Rcpp` / `RcppParallel` code paths.

## Modular architecture

`emphasis` is structured around three modules that mirror a typical diversification analysis pipeline. The README doubles as the landing page for documentation, so each module is summarised below with entry points and canonical workflows.

| Module | Purpose | Key ingredients |
|--------|---------|-----------------|
| **Module A — Preprocess** | Clean and enrich trees before inference. Handles data ingestion, tip pruning, crown-age alignment, augmentation of extinct lineages, and generation of synthetic training data. | `prune_to_extant()`, `simulate_tree(tree = ...)`, `generatePhyloPD()`, `AugmentMultiplePhyloPD()`, non-homogeneous samplers |
| **Module B — Inference** | Estimate diversification parameters with EM or differential evolution. Interfaces directly with Module A outputs. | `estimate_rates()`, `estimate_rates_control()`, `mcEM_step()`, `emphasis_de()` |
| **Module C — Goodness of fit** | Diagnose fits, compare fast/slow simulators, and run benchmark visualisations plus GAM diagnostics over augmented sweeps. | `train_GAM()`, benchmarking scripts under `analysis/`, summary plots in `results/benchmarks/` |

### Module A focus — preprocessing workflows

Module A is responsible for preparing trees so inference runs smoothly:

1. **Ingest / harmonise** trees: `prune_to_extant(tree)` removes extinct tips to obtain an extant-only backbone. `generatePhyloPD()` creates curated simulation sets with known parameters.
2. **Augment missing history** when conditioning on empirical data: `simulate_tree(tree = extant_tree, pars = ...)` fills in extinct branches (single draw). Use `AugmentMultiplePhyloPD()` for batches if you plan to fit GAM diagnostics later.
3. **Standardise metadata**: store outputs in the list format returned by `simulate_tree()` (`tes`, `tas`, `L`, `brts`, `status`) so downstream modules can consume them without extra conversion.
4. **Export training corpora**: write the `generatePhyloPD()` results to disk (e.g. via `readr::write_csv` or RDS). These corpora drive Module B parameter searches and Module C goodness-of-fit experiments (benchmark scripts already read from `results/benchmarks/`).

The remainder of this README walks through Modules B and C using these preprocessed objects, so the entire documentation story is visible on this GitHub landing page.

## Installation

```r
# Install from CRAN (when available)
# install.packages("emphasis")

# Or install the development version from GitHub
# install.packages("devtools")
devtools::install_github("franciscorichter/emphasis")
```

## Simulation module

The `simulate.R` and `generate.R` helpers cover everything from single-tree draws to large training sets. Use them to explore parameter spaces before moving on to augmentation/estimation.

### Model hierarchy

All models in `emphasis` are special cases of a single general per-lineage speciation rate:

```
lambda_i = max(0, lambda0 + betaN*N + betaP*(P/N) + betaE*d_i)
```

where **N** is species richness, **P = Σ d_j** is total phylogenetic diversity, and **d_i = t − start_date_i** is the pendant edge of lineage *i* (time since it last diverged). The extinction rate **μ** is uniform across lineages. Models are nested by progressively zeroing coefficients:

| Model | `betaN` | `betaP` | `betaE` | Mechanism |
|-------|---------|---------|---------|-----------|
| **CR** — constant rate | 0 | 0 | 0 | no dependence |
| **DD** — diversity dependent | free | 0 | 0 | clade-level richness |
| **PD** — phylogenetic diversity | free | free | 0 | richness + clade-level PD |
| **EP** — evolutionary pendant | free | free | free | richness + clade PD + per-lineage pendant |

CR ⊂ DD ⊂ PD ⊂ EP: each row is a special case of the one below it.

### Single-tree draws

`simulate_tree()` is the unified entry point for all four models.

```r
library(emphasis)

set.seed(123)

# CR — constant rate (betaN = betaP = betaE = 0)
cr_tree <- simulate_tree(pars = c(0.1, 0.4), max_t = 3, model = "cr")
cr_tree$status

# DD — diversity dependent (betaP = betaE = 0): richness N drives speciation
dd_tree <- simulate_tree(pars = c(0.1, 0.4, -0.05, 0), max_t = 5, model = "pd",
                         max_lin = 1e4, max_tries = 10)
dd_tree$status

# PD — phylogenetic diversity (betaE = 0): richness + clade-level PD
pd_tree <- simulate_tree(pars = c(0.1, 0.4, -0.05, 0.02), max_t = 5, model = "pd")
pd_tree$status

# EP — evolutionary pendant: full model, per-lineage rates via d_i
ep_tree <- simulate_tree(pars = c(0.1, 0.5, -0.02, 0.01, 0.05), max_t = 5, model = "ep")
ep_tree$status
length(ep_tree$tes$tip.label)
```

If you need direct access to specialized diagnostics (extinction sweeps, scenario grids), use the lower-level helpers documented in the reference manual; the high-level interface remains `simulate_tree()`.

### Conditional simulation — augmenting an extant tree

`simulate_tree()` also serves as the conditional-simulation entry point. Pass `tree =` an observed extant tree and it stochastically fills in the missing extinct lineages under the PD model, returning a `tas` phylo that includes both the original extant branches and the simulated extinct ones.

```r
set.seed(42)
sim <- simulate_tree(pars = c(0.1, 0.5, -0.02, 0.01), max_t = 5)
sim$status          # "done"
sim$tes             # extant-only phylo
sim$tas             # full phylo with extinct branches (from forward sim)

# Conditional simulation: augment an existing extant tree
aug <- simulate_tree(tree = sim, pars = c(0.1, 0.5, -0.02, 0.01))
aug$tes             # same extant-only phylo as sim$tes
aug$tas             # full phylo with stochastically drawn extinct lineages
length(aug$tas$tip.label) >= length(aug$tes$tip.label)  # TRUE

# Also works with a plain phylo object as input
aug2 <- simulate_tree(tree = sim$tes, pars = c(0.1, 0.5, -0.02, 0.01))
aug2$status
```

The two simulation modes share the same interface and return identical list shapes:

| `tree` argument | Mode | Produces |
|-----------------|------|----------|
| `NULL` (default) | Forward simulation | new tree from scratch (`tes`, `tas`, `L`, `brts`) |
| `simulate_tree()` result, `phylo`, or brts vector | Conditional simulation | extant tree + stochastic extinct lineages |

### Dataset factories & non-homogeneous processes

```r
training_set <- generatePhyloPD(
  n_trees = 30,
  mu_interval = c(0.05, 0.2),
  lambda_interval = c(0.2, 0.8),
  betaN_interval = c(-0.2, 0.1),
  betaP_interval = c(-0.1, 0.1),
  max_lin = 1e5,
  max_tries = 10
)

# Build a custom non-homogeneous thinning routine
cov_func <- function(t) 1 + sin(t)
rate <- function(t) ExponentialRate(matrix(c(cov_func(t)), ncol = 1), c(0.1, 0.5))
nh_draws <- nhExpRand(10, rate_func = rate, now = 0, tMax = 50)
```

`generatePhyloPD()` keeps every simulated tree plus the parameters that generated it, making it ideal for benchmarking inference pipelines or fitting surrogate models.

## Inference module

`estimate_rates()` is the unified entry point for parameter inference. It accepts a simulated or empirical tree, prunes it to extant tips if needed, and dispatches to the chosen learning method.

### The `estimate_rates()` workflow

```
simulate_tree()  →  prune_to_extant()  →  estimate_rates()
     ↓                    ↓                      ↓
 full tree          extant-only tree        parameter estimates
 (sim$tas)          (sim$tes)              (mu, lambda0, betaN, betaP)
```

`estimate_rates()` accepts any of these inputs directly — a `simulate_tree()` result, a `phylo` object, or a branching-time vector — so the prune step is usually implicit.

### Methods

The `method` argument selects the optimisation back-end:

| Method | Function | When to use |
|--------|----------|-------------|
| `"em"` | `mcEM_step()` | Fine-grained refinement; converges to a local optimum |
| `"de"` | `emphasis_de()` | Broad exploration; good seed for EM |

A robust two-stage workflow uses DE to find the basin and EM to refine:

```r
library(emphasis)

set.seed(42)
sim <- simulate_tree(pars = c(0.1, 0.5, -0.02, 0.01), max_t = 8, model = "pd")

lb <- c(0.01, 0.01, -1, -1)
ub <- c(1,    2,    1,  1)

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

# Override a subset
ctl <- estimate_rates_control("em")
ctl$sample_size <- 500
ctl$burnin      <- 50
fit <- estimate_rates(sim, lower_bound = lb, upper_bound = ub,
                     init_pars = fit_de$pars, control = ctl)
```

### Empirical trees

Pass a `phylo` object directly. Use `prune_to_extant()` first if the tree contains extinct lineages:

```r
data(bird.orders, package = "ape")
extant <- prune_to_extant(bird.orders)   # no-op if already extant-only
fit <- estimate_rates(extant, method = "de",
  lower_bound = c(0.01, 0.01, -0.5, -0.5),
  upper_bound = c(1,    2,     0.5,  0.5))
fit$pars
```

Use `get_required_sampling_size()` to adapt MCEM sample sizes, and `AugmentMultiplePhyloPD()` + `train_GAM()` to diagnose bias across simulated sweeps.

## Core modules

| Module | Purpose | Key functions |
|--------|---------|---------------|
| `simulate.R` | Unified simulator entry point | `simulate_tree()` |
| `inference.R` | Unified inference entry point | `estimate_rates()`, `estimate_rates_control()`, `prune_to_extant()` |
| `generate.R` | Dataset factories & non-homogeneous processes | `generatePhyloPD()`, `generateNonHomogeneousExp()`, `nhExpRand()`, `rate_t()`, `ExponentialRate()` |
| `augment.R` | Batch augmentation of missing lineages | `augmentPD()`, `AugmentMultiplePhyloPD()` |
| `de.R` | Differential-evolution back-end | `emphasis_de()`, `emphasis_de_factorial()` |
| `emphasis.R` | MCEM back-end | `mcEM_step()` |
| `gam.R` | Diagnostics over simulation sweeps | `train_GAM()` |
| `utils.R` | Shared helpers & plotting utilities | branching-time helpers, safe wrappers |

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

required_n <- get_required_sampling_size(results$loglik_estimation, tol = 0.05)
```

Combine diagnostics with simulation sweeps to spot parameter regions that need more intensive MCEM sampling.

## Documentation & vignettes

- Function docs live in `man/*.Rd` (generated via roxygen). Start with `?emphasis-package` for a package overview.
- Worked examples are available via `browseVignettes("emphasis")` after installing the package. They cover simulation set-ups, MCEM tuning, and diagnostics.
- The `Documentation/` directory contains PDF summaries of internal research experiments for additional context.

## Contributing

Pull requests and issues are welcome! A good contribution typically includes:

1. A reproducible example or failing test.
2. Updated documentation (`roxygen2` comments + README snippets).
3. Benchmark notes if the change affects performance-critical code.
## Citing emphasis
If you use `emphasis` in your research, please cite:

```
Francisco Richter, Thijs Janzen, Hanno Hildenbrandt. emphasis: Evolutionary Modeling on Phylogenetic Applications with Simulations and Importance Sampling. R package version 0.4.
```

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

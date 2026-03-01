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

### Single-tree draws

```r
library(emphasis)

set.seed(123)
pars <- c(mu = 0.1, lambda0 = 0.4, betaN = -0.05, betaP = 0.02)
single_tree <- sim_tree_pd_cpp(
  pars = pars,
  max_t = 5,
  max_lin = 1e4,
  max_tries = 5,
  useDDD = TRUE
)

single_tree$status   # "done", "extinct", or "too_large"
if (!is.null(single_tree$tes)) {
  tip_ids <- seq_len(min(5, length(single_tree$tes$tip.label)))
  print(single_tree$tes$tip.label[tip_ids])
}
```

Prefer a pure-R fallback? Swap to `sim_tree_pd_R(pars, max_t)` for quick debugging or pedagogical purposes.

### Extinction summaries

```r
trajectories <- sim_tree_is_extinct_pd(
  pars = pars,
  max_t = 12,
  num_repl = 100,
  max_lin = 1e4
)

subset(trajectories, break_condition != "none")
```

The tibble reports extinction flags, extinction times, phylogenetic diversity, and stopping reasons for each replicate.

### Scenario grids

```r
grid <- sim_tree_pd_grid(
  mu_vec = seq(0.05, 0.2, length.out = 3),
  lambda_vec = seq(0.2, 0.6, length.out = 3),
  b_n_vec = c(-0.3, 0),
  b_p_vec = c(-0.1, 0.1),
  max_t = 8,
  num_repl = 25,
  max_N = 5e4
)
head(grid)
```

Use grids to map parameter regions that hit extinction, explode, or stay in feasible regimes.

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

## Estimation module

Once you have simulated data (or a real phylogeny), the estimation module takes over. A typical workflow is DE seeding followed by MCEM refinement:

```r
library(emphasis)
data(bird.orders, package = "ape")
brts <- ape::branching.times(bird.orders)

# 1. Differential-evolution seed search
de_seed <- emphasis_de(
  brts = brts,
  num_iterations = 5,
  num_points = 50,
  max_missing = 1e4,
  sd_vec = rep(0.4, 4),
  lower_bound = c(0, 0, -0.5, -0.5),
  upper_bound = c(1, 2, 0.5, 0.5),
  max_lambda = 500
)

# 2. MCEM refinement
seed_pars <- colMeans(de_seed$min_pars)
fit <- mcEM_step(
  brts = brts,
  pars = seed_pars,
  sample_size = 250,
  soc = 2,
  max_missing = 1e4,
  max_lambda = 500,
  lower_bound = c(-Inf, -Inf, -Inf, -Inf),
  upper_bound = c(Inf, Inf, Inf, Inf),
  xtol = 1e-3,
  tol = 0.1,
  burnin = 20,
  num_threads = 0,
  verbose = TRUE
)

str(fit)
```

Use `get_required_sampling_size()` to adapt MCEM sample sizes, and `AugmentMultiplePhyloPD()` + `train_GAM()` to diagnose bias across simulated sweeps.

## Core modules

| Module | Purpose | Key functions |
|--------|---------|---------------|
| `simulate.R` | Core simulators for PD-dependent diversification | `sim_tree_pd_R()`, `sim_tree_pd_cpp()`, `sim_tree_is_extinct_pd()`, `sim_tree_pd_grid()` |
| `generate.R` | Dataset factories & non-homogeneous processes | `generatePhyloPD()`, `generateNonHomogeneousExp()`, `nhExpRand()`, `rate_t()`, `ExponentialRate()` |
| `augment.R` | Augmentation of missing lineages | `augmentPD()`, `AugmentMultiplePhyloPD()` |
| `de.R` | Differential-evolution seed searches | `emphasis_de()`, `emphasis_de_factorial()` |
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
- Francisco Richter (<richtf@usi.ch>)
- Thijs Janzen (<t.janzen@rug.nl>)
- Hanno Hildenbrandt (<h.hildenbrandt@rug.nl>)

## License
GPL-3. See [LICENSE](LICENSE) for details.

---

For more information, see the package documentation and vignettes, or contact the maintainer.

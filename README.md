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

## Quick start

```r
library(emphasis)

data(bird.orders, package = "ape")

# 1. Generate MCEM-ready branching times
brts <- ape::branching.times(bird.orders)

# 2. Seed the optimisation with differential evolution
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

# 3. Run MCEM using the best DE iterate
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

The snippet highlights the recommended workflow: start with DE to obtain stable seeds, then call `mcEM_step()` (or your own wrapper) for final refinement.

## Getting started guides

Load the package and explore the main entry points:

```r
library(emphasis)

# Example: Fit a diversification model to a phylogenetic tree
# (see vignettes for detailed workflows)
data(bird.orders, package = "ape")
model <- list(pars = c(lambda0 = 0.1, betaN = 0.01, betaP = 0.01))
result <- emphasis(bird.orders, model)
print(result)
```

## Core modules

| Module | Purpose | Key functions |
|--------|---------|---------------|
| `simulate.R` | Sample PD-dependent phylogenies or extinction outcomes | `sim_tree_pd_R()`, `sim_tree_pd_cpp()`, `sim_tree_is_extinct_pd()`, `sim_tree_pd_grid()` |
| `generate.R` | Convenience wrappers for building training sets or non-homogeneous processes | `generatePhyloPD()`, `generateNonHomogeneousExp()`, `nhExpRand()`, `rate_t()`, `ExponentialRate()` |
| `augment.R` | Augment trees with missing lineages | `augmentPD()`, `AugmentMultiplePhyloPD()` |
| `de.R` | Differential-evolution samplers for initial parameter search | `emphasis_de()`, `emphasis_de_factorial()` |
| `gam.R` | GAM-based diagnostics of simulation sweeps | `train_GAM()` |
| `utils.R` | Helper utilities used across the package | plotting helpers, safe wrappers |

See `man/` for in-depth documentation on each function.

## Simulation & estimation workflows

### Simulation sweeps

```r
grid <- sim_tree_pd_grid(
  mu_vec = seq(0.05, 0.2, length.out = 4),
  lambda_vec = seq(0.1, 0.5, length.out = 4),
  b_n_vec = c(-0.4, -0.2),
  b_p_vec = c(-0.1, 0),
  max_t = 10,
  num_repl = 50,
  max_N = 1e4
)

head(grid)
```

### GAM diagnostics

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

### MCEM loop control helpers

`mcEM_step()` returns a list with the MCEM trace, latest parameter estimates, iteration count, and standard error of the log-likelihood estimate. Combine it with `get_required_sampling_size()` to adaptively increase sample sizes.

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

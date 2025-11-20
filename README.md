# emphasis: Evolutionary Modeling on Phylogenetic Applications

[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](./LICENSE)
![Lifecycle: Experimental](https://img.shields.io/badge/lifecycle-experimental-orange)

## Overview

`emphasis` is an R package providing statistical tools for the study of diversification of species and the ecological drivers of macro-evolution. It implements advanced simulation and inference methods for phylogenetic analysis, including non-homogeneous Poisson processes, importance sampling, and expectation-maximization (EM) algorithms for model fitting.

## Features
- Simulation of phylogenetic trees under various diversification models
- Estimation of model parameters using EM approaches
- Support for non-homogeneous Poisson processes
- Parallel computation via RcppParallel
- Integration with the `ape` package for phylogenetic data structures

## Installation

```r
# Install from CRAN (when available)
# install.packages("emphasis")

# Or install the development version from GitHub
# install.packages("devtools")
devtools::install_github("franciscorichter/emphasis")
```

## Getting Started

Load the package and check out the main function:

```r
library(emphasis)

# Example: Fit a diversification model to a phylogenetic tree
# (see vignettes for detailed workflows)
data(bird.orders, package = "ape")
model <- list(pars = c(lambda0 = 0.1, betaN = 0.01, betaP = 0.01))
result <- emphasis(bird.orders, model)
print(result)
```

## Main Functions
- `emphasis()`: Main interface for fitting diversification models using EM.
- `simulate_evolution()`: Simulate phylogenetic trees under various models.
- `loglikelihood()`, `augmentPD()`, `generateNonHomogeneousExp()`, `nhExpRand()`: Utility and simulation functions.

## Simulation Methods

The emphasis package provides several simulation tools for evolutionary and phylogenetic modeling:

- `simulate_evolution()`: Simulate phylogenetic trees under a variety of diversification models (stub; see vignettes for examples).
- `generateNonHomogeneousExp()`: Simulate events from a non-homogeneous Poisson process (stub).
- `nhExpRand()`: Sample event times from a non-homogeneous exponential process (stub).

**Example Usage:**
```r
# Example (pseudo-code, see vignettes for real usage):
model <- list(pars = c(lambda0 = 0.1, betaN = 0.01, betaP = 0.01))
sim_tree <- simulate_evolution(model)
```

## Estimation Methods

emphasis supports several estimation approaches for model fitting:

### MCEM Estimation
- **Function:** `estimation_mcem()`
- **Description:** Fits a diversification model to a phylogenetic tree using the Monte Carlo Expectation-Maximization algorithm.
- **Example:**
```r
result <- estimation_mcem(phylo, model)
```

### MLE Estimation
- **Function:** `estimation_mle()`
- **Description:** Fits a diversification model using Maximum Likelihood Estimation.
- **Example:**
```r
result <- estimation_mle(phylo, model)
```

## Vignettes & Examples

- [Simulation](vignettes/Simulation.html): Introduction to simulating non-homogeneous Poisson processes and diversification models.
- [Estimation](vignettes/documentation_estimation.html): Parameter estimation in phylogenetic diversification models.

Run the vignettes in RStudio or with:
```r
browseVignettes("emphasis")
```

**How to contribute:**
- Fork the repo and submit pull requests
- Open issues for bugs or feature requests
- Add new models, methods, or vignettes!
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

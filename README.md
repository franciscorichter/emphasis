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

See the [reference manual](man/) and [vignettes](vignettes/) for all available functions and detailed documentation.

## Vignettes & Examples

- [Simulation](vignettes/Simulation.html): Introduction to simulating non-homogeneous Poisson processes and diversification models.
- [Estimation](vignettes/documentation_estimation.html): Parameter estimation in phylogenetic diversification models.

Run the vignettes in RStudio or with:
```r
browseVignettes("emphasis")
```

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

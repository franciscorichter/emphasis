\name{emphasis-package}
\alias{emphasis-package}
\docType{package}
\title{
  Evolutionary Modeling on Phylogenetic Applications with Simulations and Importance Sampling
}
\description{
  Tools for simulating phylogenetic trees, augmenting missing lineages, and fitting diversification models using Monte Carlo EM, Cross-Entropy Method, and importance sampling.
}
\details{
  The package revolves around four workflows:
  \itemize{
    \item \strong{Simulation}: `simulate_tree()` for forward and conditional simulation.
    \item \strong{Augmentation}: `mc_loglik()` for Monte Carlo log-likelihood via importance sampling.
    \item \strong{Estimation}: `estimate_rates()`, `emphasis_cem()` (CEM back-end).
    \item \strong{Diagnostics}: `train_GAM()` plus plotting utilities.
  }
  Heavy computations are delegated to C++ via `Rcpp` / `RcppParallel` for performance.
}
\author{
Francisco Richter <richtf@usi.ch>

Contributors: Thijs Janzen, Hanno Hildenbrandt
}
\keyword{ package }
\seealso{
  \code{vignette("Simulation", package = "emphasis")},
  \code{vignette("documentation_estimation", package = "emphasis")}
}
\examples{
  \dontrun{
    data(bird.orders, package = "ape")
    result <- emphasis(bird.orders, model = list(pars = c(0.1, 0.2, -0.1, 0)))
    str(result)
  }
}

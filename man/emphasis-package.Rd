\name{emphasis-package}
\alias{emphasis-package}
\docType{package}
\title{
  Evolutionary Modeling on Phylogenetic Applications with Simulations and Importance Sampling
}
\description{
  Tools for simulating phylogenetic trees, augmenting missing lineages, and fitting diversification models using Monte Carlo EM, differential evolution, and importance sampling.
}
\details{
  The package revolves around four workflows:
  \itemize{
    \item \strong{Simulation}: `sim_tree_pd_*`, `generatePhyloPD()` and helpers.
    \item \strong{Augmentation}: `augmentPD()` and `AugmentMultiplePhyloPD()`.
    \item \strong{Estimation}: `emphasis_de()`, `emphasis_de_factorial()`, `mcEM_step()`.
    \item \strong{Diagnostics}: `train_GAM()` plus plotting utilities in `utils.R`.
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

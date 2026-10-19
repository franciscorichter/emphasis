#' emphasis: Monte Carlo Inference for Phylogenetic Diversification Models
#'
#' Birth-death models of species diversification whose speciation and
#' extinction rates depend on covariates the process generates: the number of
#' lineages \eqn{N}, the age-imbalance \eqn{D}, the evolutionary
#' distinctiveness \eqn{ED} and the number of ancestral splits \eqn{K}.  A
#' reconstructed phylogeny records only the extant species; the likelihood is
#' integrated over the unobserved extinct (and, at incomplete sampling,
#' unsampled) lineages by importance sampling on augmented trees.
#'
#' The main functions:
#' \itemize{
#'   \item \code{\link{simulate_tree}}: forward simulation under a model, and
#'     augmentation of an observed tree.
#'   \item \code{\link{estimate_rates}}: maximum-likelihood estimation by Monte
#'     Carlo expectation-maximization (\code{"mcem"}), the cross-entropy
#'     method (\code{"cem"}) or a GAM fitted to a likelihood surface
#'     (\code{"gam"}).
#'   \item \code{\link{compare_models}}: AIC comparison of fits that accounts
#'     for their Monte Carlo error.
#'   \item \code{\link{lineage_table}}, \code{\link{covariate_path}}: the
#'     per-lineage design of a tree and a selection path over the covariates.
#'   \item \code{\link{diagnose_mcem}}, \code{\link{diagnose_cem}},
#'     \code{\link{diagnose_gam}}: convergence and importance-sampling
#'     diagnostics.
#' }
#'
#' Computation is in C++ through \pkg{Rcpp} and \pkg{RcppParallel}; every
#' function runs on one thread unless \code{num_threads} is set.
#'
#' @references Richter F, Haegeman B, Etienne RS, Wit EC (2020). Introducing a
#'   general class of species diversification models for phylogenetic trees.
#'   \emph{Statistica Neerlandica} 74(3), 261--274.
#'   \doi{10.1111/stan.12205}
#' @keywords internal
"_PACKAGE"

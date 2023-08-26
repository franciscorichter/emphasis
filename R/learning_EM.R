#' function to perform one step of the E-M algorithm
#' @param brts vector of branching times
#' @param init_pars vector of initial parameter files
#' @param sample_size number of samples
#' @param maxN maximum number of failed trees
#' @param soc number of lineages at the root/crown (1/2)
#' @param max_missing maximum number of species missing
#' @param max_lambda maximum speciation rate
#' @param lower_bound vector of lower bound values for optimization, should
#' be equal in length to the vector of init_pars
#' @param upper_bound vector of upper bound values for optimization, should
#' be equal in length to the vector of init_pars
#' @param xtol_rel relative tolerance for optimization
#' @param num_threads number of threads used.
#' @return a list with the following components: 
#' \itemize{
#'  \item{trees}{list of trees}
#'  \item{rejected}{number of rejected trees}
#'  \item{rejected_overruns}{number of trees rejected due to too large size}
#'  \item{rejected_lambda}{number of trees rejected due to lambda errors}
#'  \item{rejected_zero_weights}{number of trees rejected to zero weight}
#'  \item{time_elapsed}{time used}
#'  \item{weights}{vector of weights}
#'  \item{fhat}{vector of fhat values}
#'  \item{logf}{vector of logf values}
#'  \item{logg}{vector of logg values}
#' }
#' @export
e_cpp <- function(brts, init_pars, sample_size, maxN, soc, max_missing, max_lambda, lower_bound, upper_bound, xtol_rel, num_threads=1) {
  .Call('_emphasis_rcpp_mce', PACKAGE = 'emphasis', brts, init_pars, sample_size, maxN, soc, max_missing, max_lambda, lower_bound, upper_bound, xtol_rel, num_threads)
}

augmentPD <- function(phylo, pars, maxN, max_missing, lower_bound, upper_bound) {
    num_threads=1
    sample_size = 1
    init_pars = pars
    soc = 2 
    max_lambda = max_missing*100
    xtol_rel = 0.00001
    brts = ape::branching.times(phylo)
    # brts = sort(max(brts) - brts)
    # brts = brts[-1]
  .Call('_emphasis_rcpp_mce', PACKAGE = 'emphasis', brts, init_pars, sample_size, maxN, soc, max_missing, max_lambda, lower_bound, upper_bound, xtol_rel, num_threads)
}

rcpp_mce_grid <- function(pars_R, brts, sample_size, maxN, soc, max_missing, max_lambda, lower_bound, upper_bound, xtol_rel, num_threads) {
  .Call('_emphasis_rcpp_mce_grid', PACKAGE = 'emphasis', pars_R, brts, sample_size, maxN, soc, max_missing, max_lambda, lower_bound, upper_bound, xtol_rel, num_threads)
}

rcpp_mce_grid_factorial <- function(pars_R, brts, sample_size, maxN, soc, max_missing, max_lambda, lower_bound, upper_bound, xtol_rel, num_threads) {
  .Call('_emphasis_rcpp_mce_grid_factorial', PACKAGE = 'emphasis', pars_R, brts, sample_size, maxN, soc, max_missing, max_lambda, lower_bound, upper_bound, xtol_rel, num_threads)
}

#' function to perform one step of the E-M algorithm
#' @param brts vector of branching times
#' @param init_pars vector of initial parameter files
#' @param sample_size number of samples
#' @param maxN maximum number of failed trees
#' @param soc number of lineages at the root/crown (1/2)
#' @param max_missing maximum number of species missing
#' @param max_lambda maximum speciation rate
#' @param lower_bound vector of lower bound values for optimization, should
#' be equal in length to the vector of init_pars
#' @param upper_bound vector of upper bound values for optimization, should
#' be equal in length to the vector of init_pars
#' @param xtol_rel relative tolerance for optimization
#' @param num_threads number of threads used.
#' @param copy_trees if set to true, the trees generated are returned as well
#' @param rconditional R function that evaluates the GAM function.
#' @return a list with the following components: 
#' \itemize{
#'  \item{trees}{list of trees}
#'  \item{rejected}{number of rejected trees}
#'  \item{rejected_overruns}{number of trees rejected due to too large size}
#'  \item{rejected_lambda}{number of trees rejected due to lambda errors}
#'  \item{rejected_zero_weights}{number of trees rejected to zero weight}
#'  \item{estimates}{vector of estimates}
#'  \item{nlopt}{nlopt status}
#'  \item{fhat}{vector of fhat values}
#'  \item{time}{time elapsed}
#'  \item{weights}{vector of weights}
#' }
#' @export
em_cpp <- function(brts, init_pars, sample_size, maxN, soc, max_missing, max_lambda, lower_bound, upper_bound, xtol_rel, num_threads, copy_trees, rconditional = NULL) {
  .Call('_emphasis_rcpp_mcem', PACKAGE = 'emphasis', brts, init_pars, sample_size, maxN, soc, max_missing, max_lambda, lower_bound, upper_bound, xtol_rel, num_threads, copy_trees, rconditional)
}

#' function to perform one step of the E-M algorithm
#' @param e_step result of e_step function, as a list
#' @param init_pars vector of initial parameter values
#' @param plugin string indicating plugin used, currently available: 'rpd1' and
#' 'rpd5c'
#' @param lower_bound vector of lower bound values for optimization, should
#' be equal in length to the vector of init_pars
#' @param upper_bound vector of upper bound values for optimization, should
#' be equal in length to the vector of init_pars 
#' @param xtol_rel relative tolerance for optimization
#' @param num_threads number of threads used.
#' @param rconditional R function that evaluates the GAM function.
#' @return list with the following entries:
#' \itemize{
#'  \item{estimates}{vector of estimates}
#'  \item{nlopt}{nlopt status}
#'  \item{time}{used computation time} 
#' }
#' @export
m_cpp <- function(e_step, init_pars, plugin, lower_bound, upper_bound, xtol_rel, num_threads, rconditional = NULL) {
  .Call('_emphasis_rcpp_mcm', PACKAGE = 'emphasis', e_step, init_pars, plugin, lower_bound, upper_bound, xtol_rel, num_threads, rconditional)
}
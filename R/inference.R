#' Prune extinct tips from a phylogenetic tree
#'
#' Drops tips that did not survive to the present (i.e. tips whose root-to-tip
#' distance is shorter than the maximum, within a small numerical tolerance).
#'
#' @param phy A \code{phylo} object, possibly containing extinct lineages.
#' @param tol Tolerance for deciding whether a tip is at the present.
#'   Default \code{1e-8}.
#' @return A \code{phylo} object with only extant tips retained.
#' @examples
#' \dontrun{
#' sim <- simulate_tree(c(0.1, 0.5, -0.02, 0.01), max_t = 5, model = "pd")
#' extant <- prune_to_extant(sim$tas)
#' }
#' @export
prune_to_extant <- function(phy, tol = 1e-8) {
  if (!inherits(phy, "phylo")) stop("'phy' must be a phylo object.")
  heights  <- ape::node.depth.edgelength(phy)
  tip_h    <- heights[seq_len(ape::Ntip(phy))]
  extinct  <- which(abs(tip_h - max(tip_h)) > tol)
  if (length(extinct) > 0) phy <- ape::drop.tip(phy, phy$tip.label[extinct])
  phy
}

#' Default control parameters for estimate_rates
#'
#' Returns a list of default tuning parameters for \code{\link{estimate_rates}}.
#' Override individual elements and pass the modified list via the
#' \code{control} argument.
#'
#' @param method Which method's defaults to return: \code{"em"} or \code{"de"}.
#' @param n_pars Number of parameters (used to set the default \code{sd_vec}
#'   for DE when \code{NULL}).
#' @return A named list of control parameters.
#'
#' Shared parameters (both methods):
#' \describe{
#'   \item{\code{soc}}{Stem (1) or crown (2) age condition. Default \code{2}.}
#'   \item{\code{max_missing}}{Max missing lineages tolerated during
#'     augmentation. Default \code{1e4}.}
#'   \item{\code{max_lambda}}{Max speciation rate for augmentation.
#'     Default \code{500}.}
#'   \item{\code{num_threads}}{Threads for parallel computation. Default
#'     \code{1}.}
#' }
#'
#' EM-specific parameters:
#' \describe{
#'   \item{\code{sample_size}}{MC sample size per iteration. Default \code{200}.}
#'   \item{\code{xtol}}{Relative tolerance for the M-step optimiser.
#'     Default \code{1e-3}.}
#'   \item{\code{tol}}{Standard error of log-likelihood below which EM
#'     is considered converged. Default \code{0.1}.}
#'   \item{\code{burnin}}{Number of burn-in iterations before convergence
#'     is tested. Default \code{20}.}
#' }
#'
#' DE-specific parameters:
#' \describe{
#'   \item{\code{num_iterations}}{Number of DE iterations. Default \code{20}.}
#'   \item{\code{num_points}}{Particle population size. Default \code{50}.}
#'   \item{\code{sd_vec}}{Initial perturbation SDs. \code{NULL} (default)
#'     sets each to \code{(upper - lower) / 4} at call time.}
#'   \item{\code{maxN}}{Max augmentation attempts per particle. Default
#'     \code{10}.}
#'   \item{\code{disc_prop}}{Fraction of best particles retained per
#'     iteration. Default \code{0.5}.}
#' }
#' @export
estimate_rates_control <- function(method = c("em", "de"), n_pars = 4) {
  method <- match.arg(method)
  common <- list(
    soc         = 2L,
    max_missing = 1e4,
    max_lambda  = 500,
    num_threads = 1L
  )
  if (method == "em") {
    c(common, list(
      sample_size = 200L,
      xtol        = 1e-3,
      tol         = 0.1,
      burnin      = 20L
    ))
  } else {
    c(common, list(
      num_iterations = 20L,
      num_points     = 50L,
      sd_vec         = NULL,
      maxN           = 10L,
      disc_prop      = 0.5
    ))
  }
}

# --------------------------------------------------------------------------- #
#  Internal helpers                                                            #
# --------------------------------------------------------------------------- #

#' @keywords internal
.extract_brts <- function(tree) {
  if (is.numeric(tree)) {
    return(sort(tree, decreasing = TRUE))
  }
  if (inherits(tree, "phylo")) {
    return(sort(ape::branching.times(tree), decreasing = TRUE))
  }
  # simulate_tree() result — prefer $brts (already sorted from DDD::L2brts)
  if (is.list(tree) && !inherits(tree, "phylo")) {
    if (!is.null(tree$brts) && is.numeric(tree$brts))
      return(sort(tree$brts, decreasing = TRUE))
    if (!is.null(tree$tes) && inherits(tree$tes, "phylo"))
      return(sort(ape::branching.times(tree$tes), decreasing = TRUE))
    if (!is.null(tree$tas) && inherits(tree$tas, "phylo"))
      return(sort(ape::branching.times(prune_to_extant(tree$tas)), decreasing = TRUE))
    stop("List passed to estimate_rates() must contain '$brts', '$tes', or '$tas'.")
  }
  stop("'tree' must be a phylo object, a simulate_tree() result, or a numeric branching-time vector.")
}

#' @keywords internal
.par_names <- function(n) {
  nms <- c("mu", "lambda0", "betaN", "betaP", "betaE")
  if (n <= length(nms)) nms[seq_len(n)] else paste0("par", seq_len(n))
}

#' @keywords internal
.run_em <- function(brts, init_pars, lower_bound, upper_bound, ctrl) {
  raw <- mcEM_step(
    brts        = brts,
    pars        = init_pars,
    sample_size = ctrl$sample_size,
    soc         = ctrl$soc,
    max_missing = ctrl$max_missing,
    max_lambda  = ctrl$max_lambda,
    lower_bound = lower_bound,
    upper_bound = upper_bound,
    xtol        = ctrl$xtol,
    tol         = ctrl$tol,
    burnin      = ctrl$burnin,
    num_threads = ctrl$num_threads,
    verbose     = ctrl$verbose
  )
  fhat_vec <- raw$mcem$fhat
  loglik   <- if (any(is.finite(fhat_vec)))
    utils::tail(fhat_vec[is.finite(fhat_vec)], 1L) else NA_real_
  list(pars = raw$pars, loglik = loglik, details = raw)
}

#' @keywords internal
.run_de <- function(brts, lower_bound, upper_bound, ctrl) {
  sd_vec <- ctrl$sd_vec
  if (is.null(sd_vec)) sd_vec <- (upper_bound - lower_bound) / 4
  raw <- emphasis_de(
    brts           = brts,
    num_iterations = ctrl$num_iterations,
    num_points     = ctrl$num_points,
    max_missing    = ctrl$max_missing,
    sd_vec         = sd_vec,
    lower_bound    = lower_bound,
    upper_bound    = upper_bound,
    maxN           = ctrl$maxN,
    max_lambda     = ctrl$max_lambda,
    disc_prop      = ctrl$disc_prop,
    verbose        = ctrl$verbose,
    num_threads    = ctrl$num_threads
  )
  # min_loglik stores -fhat of best particle (vals = -fhat in emphasis_de)
  loglik <- if (length(raw$minloglik) > 0)
    -utils::tail(raw$minloglik, 1L) else NA_real_
  list(pars = raw$obtained_estim, loglik = loglik, details = raw)
}

# --------------------------------------------------------------------------- #
#  Public API                                                                  #
# --------------------------------------------------------------------------- #

#' Estimate diversification rates from a phylogenetic tree
#'
#' Unified entry point for parameter inference under PD-dependent
#' diversification. Accepts a simulated or empirical tree, extracts branching
#' times (pruning to extant tips if needed), and dispatches to the chosen
#' optimisation back-end.
#'
#' @param tree One of:
#'   \itemize{
#'     \item a \code{phylo} object (extant tips only recommended),
#'     \item the list returned by \code{\link{simulate_tree}} (\code{$tes} is
#'       used automatically, falling back to a pruned \code{$tas}), or
#'     \item a numeric vector of branching times (crown-age ordering).
#'   }
#' @param method Optimisation method:
#'   \describe{
#'     \item{\code{"em"}}{Monte Carlo Expectation-Maximisation (default).
#'       Iterates E- and M-steps until the log-likelihood standard error
#'       falls below \code{control$tol}. Requires \code{init_pars}.}
#'     \item{\code{"de"}}{Differential Evolution. Initialises a particle
#'       population randomly within the bounds and evolves it over
#'       \code{control$num_iterations} iterations. Does not need
#'       \code{init_pars}.}
#'   }
#' @param lower_bound Numeric vector of lower parameter bounds
#'   (\code{c(mu, lambda0, betaN, betaP)} for the PD model).
#' @param upper_bound Numeric vector of upper parameter bounds (same order).
#' @param init_pars Starting parameter vector. Required for \code{"em"};
#'   ignored for \code{"de"}. If \code{NULL} with \code{"em"}, the midpoint
#'   of the bounds is used.
#' @param control Named list of tuning parameters. Use
#'   \code{\link{estimate_rates_control}} to inspect and modify defaults.
#' @param verbose Print progress messages (\code{TRUE}/\code{FALSE}).
#' @return A list with:
#'   \describe{
#'     \item{\code{pars}}{Named numeric vector of parameter estimates.}
#'     \item{\code{loglik}}{Final log-likelihood estimate.}
#'     \item{\code{method}}{Method used (\code{"em"} or \code{"de"}).}
#'     \item{\code{details}}{Full output from the back-end function
#'       (\code{mcEM_step} or \code{emphasis_de}).}
#'   }
#' @examples
#' \dontrun{
#' set.seed(42)
#' sim <- simulate_tree(c(0.1, 0.5, -0.02, 0.01), max_t = 5, model = "pd")
#'
#' # EM estimation (starts from parameter midpoint)
#' fit_em <- estimate_rates(sim, method = "em",
#'   lower_bound = c(0, 0, -1, -1),
#'   upper_bound = c(1, 2,  1,  1))
#' fit_em$pars
#'
#' # DE estimation (no init_pars needed)
#' fit_de <- estimate_rates(sim, method = "de",
#'   lower_bound = c(0, 0, -1, -1),
#'   upper_bound = c(1, 2,  1,  1),
#'   control = list(num_iterations = 10, num_points = 30))
#' fit_de$pars
#' }
#' @export
estimate_rates <- function(tree,
                           method      = c("em", "de"),
                           lower_bound,
                           upper_bound,
                           init_pars   = NULL,
                           control     = list(),
                           verbose     = FALSE) {
  method <- match.arg(method)

  if (!is.numeric(lower_bound) || !is.numeric(upper_bound))
    stop("'lower_bound' and 'upper_bound' must be numeric vectors.")
  if (length(lower_bound) != length(upper_bound))
    stop("'lower_bound' and 'upper_bound' must have the same length.")

  brts <- .extract_brts(tree)

  if (is.null(init_pars) && method == "em")
    init_pars <- (lower_bound + upper_bound) / 2

  defaults <- estimate_rates_control(method, n_pars = length(lower_bound))
  ctrl     <- utils::modifyList(defaults, control)
  ctrl$verbose <- isTRUE(verbose) || isTRUE(ctrl$verbose)

  raw <- switch(method,
    em = .run_em(brts, init_pars, lower_bound, upper_bound, ctrl),
    de = .run_de(brts, lower_bound, upper_bound, ctrl)
  )

  pars <- stats::setNames(as.numeric(raw$pars),
                          .par_names(length(lower_bound)))
  list(pars = pars, loglik = raw$loglik, method = method, details = raw$details)
}

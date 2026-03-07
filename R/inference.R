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
#' @keywords internal
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
#' @param method Which method's defaults to return: \code{"mcem"} or \code{"cem"}.
#' @param n_pars Number of parameters (used to set the default \code{sd_vec}
#'   for CEM when \code{NULL}).
#' @return A named list of control parameters.
#'
#' Shared parameters (both methods):
#' \describe{
#'   \item{\code{lower_bound}}{Numeric vector of lower parameter bounds.
#'     Length must be \code{2 + 2*sum(model)}. \strong{Required} (no default).}
#'   \item{\code{upper_bound}}{Numeric vector of upper parameter bounds
#'     (same order). \strong{Required} (no default).}
#'   \item{\code{verbose}}{Print progress messages. Default \code{FALSE}.}
#'   \item{\code{max_missing}}{Max extinct lineages per augmented tree.
#'     Trees that require more than this many extinct lineages are discarded.
#'     Default \code{1e4}.}
#'   \item{\code{num_threads}}{Threads for parallel computation. Default
#'     \code{1}.}
#' }
#'
#' MCEM-specific parameters:
#' \describe{
#'   \item{\code{sampling}}{Sampling scheme. Currently only
#'     \code{"dynamic_fresh"} is implemented.}
#'   \item{\code{num_trees}}{Augmented trees per EM iteration. Default
#'     \code{200}. Alias: \code{sample_size}.}
#'   \item{\code{xtol}}{Relative tolerance for the M-step optimiser.
#'     Default \code{1e-3}.}
#'   \item{\code{max_iter}}{Maximum EM iterations. Default \code{200}.}
#'   \item{\code{maxN}}{Total augmentation attempts per E-step (accepted
#'     plus rejected). Must exceed \code{num_trees}; increase for models
#'     with high rejection rates (e.g. PD/DD with strongly negative
#'     covariate slope). Default \code{2000}.}
#'   \item{\code{tol}}{Parameter-stability convergence threshold. MCEM
#'     is considered converged when the largest relative parameter change
#'     (scaled by the search range) is below \code{tol} for \code{patience}
#'     consecutive iterations. Default \code{1e-3}.}
#'   \item{\code{patience}}{Number of consecutive iterations with parameter
#'     change below \code{tol} required to declare convergence. Default
#'     \code{3}.}
#' }
#'
#' CEM-specific parameters:
#' \describe{
#'   \item{\code{max_iter}}{Maximum CEM iterations (hard cap). Default \code{20}.}
#'   \item{\code{num_particles}}{Particle population size (parameter vectors
#'     per iteration). Default \code{50}. Alias: \code{num_points}.}
#'   \item{\code{sd_vec}}{Initial perturbation SDs. \code{NULL} (default)
#'     sets each to \code{(upper - lower) / 4} at call time.}
#'   \item{\code{maxN}}{Max augmentation attempts per particle evaluation.
#'     Default \code{10}.}
#'   \item{\code{num_trees}}{Augmented trees simulated per particle per
#'     iteration. Default \code{1}. Alias: \code{sample_size}.}
#'   \item{\code{shared_trees}}{If \code{FALSE} (default), each particle is
#'     evaluated on its own simulated trees.  If \code{TRUE}, trees are pooled
#'     across all particles and every particle's \code{fhat} is re-evaluated
#'     over the full pool (lower variance; bias correction not yet implemented).}
#'   \item{\code{disc_prop}}{Elite fraction retained per iteration (Cross-Entropy
#'     Method). Default \code{0.5}.}
#'   \item{\code{bias_correct}}{Apply moment-based bias correction to the IS
#'     log-likelihood. Default \code{FALSE}.}
#'   \item{\code{num_trees}}{When \code{> 1}, bootstrap variance of the
#'     log-likelihood is computed automatically (B = 200 replicates) from
#'     the IS weights already drawn.  Minimum recommendation: \code{5}.
#'     \code{loglik_var} in the fit is \code{NA} when \code{num_trees = 1}.}
#' }
#' @keywords internal
estimate_rates_control <- function(method = c("mcem", "cem", "gam"), n_pars = 4) {
  method <- match.arg(method)
  common <- list(
    lower_bound = NULL,
    upper_bound = NULL,
    verbose     = FALSE,
    max_missing = 1e4,
    num_threads = 1L
  )
  if (method == "mcem") {
    c(common, list(
      sampling    = "dynamic_fresh",
      sample_size = 200L,       # alias: num_trees
      max_iter    = 200L,
      maxN        = 2000L,      # total augmentation attempts; must exceed sample_size
      xtol        = 1e-3,
      tol         = 1e-3,
      patience    = 3L
    ))
  } else if (method == "cem") {
    c(common, list(
      max_iter     = 20L,
      num_particles = 50L,      # alias: num_points
      sd_vec       = NULL,
      maxN         = 10L,
      num_trees    = 1L,        # alias: sample_size
      shared_trees = FALSE,
      disc_prop    = 0.5,
      tol          = 1e-4,
      patience     = 5L,
      bias_correct = FALSE
    ))
  } else {
    # GAM method
    c(common, list(
      grid_points  = 20L,       # per-parameter grid resolution
      sample_size  = 200L,      # IS trees per grid point
      maxN         = NULL,      # auto-set from sample_size
      max_lambda   = 500,
      spline_type  = "univariate"
    ))
  }
}

# --------------------------------------------------------------------------- #
#  Internal helpers                                                            #
# --------------------------------------------------------------------------- #

#' Resolve control-parameter aliases
#'
#' Supports both old names (\code{num_points}, \code{sample_size}) and new
#' names (\code{num_particles}, \code{num_trees}).  User-supplied values take
#' precedence; if both old and new are present, the new name wins.
#' Internally the code always reads the OLD names for backwards compatibility
#' with \code{emphasis_cem} and the C++ layer.
#' @keywords internal
.resolve_control_aliases <- function(ctrl, method, user_ctrl = list()) {
  # Sync alias pairs.  If the user explicitly supplied an old name, it

  # overrides the default new name (and vice-versa).  When both old and
  # new names come from defaults only, the new name is authoritative.
  .sync <- function(ctrl, new_nm, old_nm, user_ctrl) {
    user_new <- !is.null(user_ctrl[[new_nm]])
    user_old <- !is.null(user_ctrl[[old_nm]])
    if (user_old && !user_new) {
      # User supplied old name only -> propagate to new
      ctrl[[new_nm]] <- ctrl[[old_nm]]
    } else if (user_new && !user_old) {
      # User supplied new name only -> propagate to old
      ctrl[[old_nm]] <- ctrl[[new_nm]]
    } else {
      # Both from defaults, or both from user -> new name wins; sync to old
      ctrl[[old_nm]] <- ctrl[[new_nm]]
    }
    ctrl
  }
  if (method == "cem") {
    ctrl <- .sync(ctrl, "num_particles", "num_points",  user_ctrl)
    ctrl <- .sync(ctrl, "num_trees",     "sample_size", user_ctrl)
  }
  if (method == "mcem") {
    ctrl <- .sync(ctrl, "num_trees", "sample_size", user_ctrl)
  }
  ctrl
}

#' @keywords internal
.extract_brts <- function(tree) {
  if (is.numeric(tree)) {
    return(sort(tree, decreasing = TRUE))
  }
  if (inherits(tree, "phylo")) {
    return(sort(ape::branching.times(tree), decreasing = TRUE))
  }
  # simulate_tree() result: try $tes, then $tas (pruned to extant)
  if (is.list(tree) && !inherits(tree, "phylo")) {
    if (!is.null(tree$tes) && inherits(tree$tes, "phylo"))
      return(sort(ape::branching.times(tree$tes), decreasing = TRUE))
    if (!is.null(tree$tas) && inherits(tree$tas, "phylo"))
      return(sort(ape::branching.times(prune_to_extant(tree$tas)), decreasing = TRUE))
    stop("List passed to estimate_rates() must contain '$tes' or '$tas'.")
  }
  stop("'tree' must be a phylo object, a simulate_tree() result, or a numeric branching-time vector.")
}

#' @keywords internal
.par_names <- function(model_bin) {
  lam <- c("beta_0",
           if (model_bin[1]) "beta_N",
           if (model_bin[2]) "beta_P",
           if (model_bin[3]) "beta_E")
  mu <- c("gamma_0",
          if (model_bin[1]) "gamma_N",
          if (model_bin[2]) "gamma_P",
          if (model_bin[3]) "gamma_E")
  c(lam, mu)
}

#' @keywords internal
.contract_pars <- function(pars8, model_bin) {
  active <- which(model_bin == 1L)
  lam <- c(pars8[1L], pars8[active + 1L])
  mu  <- c(pars8[5L], pars8[active + 5L])
  c(lam, mu)
}

#' @keywords internal
.run_mcem <- function(brts, init_pars, lower_bound, upper_bound, ctrl, model = c(0L, 0L, 0L), link = 0L) {
  raw <- switch(ctrl$sampling,
    dynamic_fresh   = .mcem_dynamic_fresh(
      brts        = brts,
      pars        = init_pars,
      sample_size = ctrl$sample_size,
      maxN        = ctrl$maxN,
      max_missing = ctrl$max_missing,
      lower_bound = lower_bound,
      upper_bound = upper_bound,
      max_iter    = ctrl$max_iter,
      xtol        = ctrl$xtol,
      tol         = ctrl$tol,
      patience    = ctrl$patience,
      num_threads = ctrl$num_threads,
      verbose     = ctrl$verbose,
      model       = model,
      link        = link
    ),
    dynamic_recycle = stop("'dynamic_recycle' sampling not yet implemented."),
    dynamic_mixed   = stop("'dynamic_mixed' sampling not yet implemented."),
    stop("Unknown sampling scheme: ", ctrl$sampling)
  )
  fhat_vec <- raw$mcem$fhat
  loglik   <- if (any(is.finite(fhat_vec)))
    utils::tail(fhat_vec[is.finite(fhat_vec)], 1L) else NA_real_
  loglik_var <- if (!is.null(raw$loglik_var)) raw$loglik_var else NA_real_
  list(pars = raw$pars, loglik = loglik, loglik_var = loglik_var, details = raw)
}

#' @keywords internal
.run_cem <- function(brts, lower_bound, upper_bound, ctrl, model = c(0L, 0L, 0L), link = 0L) {
  sd_vec <- ctrl$sd_vec
  if (is.null(sd_vec)) sd_vec <- (upper_bound - lower_bound) / 4
  raw <- emphasis_cem(
    brts         = brts,
    max_iter     = ctrl$max_iter,
    num_points   = ctrl$num_points,
    max_missing  = ctrl$max_missing,
    sd_vec       = sd_vec,
    lower_bound  = lower_bound,
    upper_bound  = upper_bound,
    maxN         = ctrl$maxN,
    sample_size  = ctrl$sample_size,
    shared_trees = ctrl$shared_trees,
    disc_prop    = ctrl$disc_prop,
    tol          = ctrl$tol,
    patience     = ctrl$patience,
    bias_correct = ctrl$bias_correct,
    verbose      = ctrl$verbose,
    num_threads  = ctrl$num_threads,
    model        = model,
    link         = link
  )
  if (identical(raw$converged, "all_failed"))
    warning("CEM failed: all particles were rejected at every attempt. ",
            "Try increasing num_particles (e.g. >= 20) or widening the ",
            "parameter bounds.")
  # Use the final iteration's best fhat, not the global max over all iterations.
  # The global max picks lucky single-tree spikes; the final iteration reflects
  # the converged population.
  loglik <- if (length(raw$best_loglik) > 0)
    raw$best_loglik[length(raw$best_loglik)] else NA_real_
  list(pars = raw$obtained_estim, loglik = loglik,
       loglik_var = raw$loglik_var, details = raw)
}

#' @keywords internal
.run_gam <- function(brts, lower_bound, upper_bound, ctrl,
                     model = c(0L, 0L, 0L), link = 0L) {
  if (!requireNamespace("mgcv", quietly = TRUE))
    stop("Package 'mgcv' is required for method=\"gam\". ",
         "Install it with install.packages(\"mgcv\").")

  model_bin <- as.integer(model)
  link_int  <- as.integer(link)
  par_names <- .par_names(model_bin)
  n_par     <- length(par_names)
  compact_lb <- .contract_pars(lower_bound, model_bin)
  compact_ub <- .contract_pars(upper_bound, model_bin)

  # Build parameter grid (uniform grid per dimension)
  gp <- as.integer(ctrl$grid_points)
  grid_list <- lapply(seq_len(n_par), function(j) {
    seq(compact_lb[j], compact_ub[j], length.out = gp)
  })
  names(grid_list) <- par_names
  pars_grid <- as.matrix(expand.grid(grid_list))

  if (ctrl$verbose)
    cat(sprintf("  GAM grid: %d points (%d per parameter, %d parameters)\n",
                nrow(pars_grid), gp, n_par))

  # Evaluate IS log-likelihood at each grid point
  surface <- estimate_likelihood_surface(
    tree        = brts,
    pars_mat    = pars_grid,
    model       = model_bin,
    sample_size = as.integer(ctrl$sample_size),
    link        = if (link_int == 0L) "linear" else "exponential",
    max_missing = ctrl$max_missing,
    max_lambda  = ctrl$max_lambda,
    maxN        = ctrl$maxN,
    num_threads = ctrl$num_threads,
    verbose     = ctrl$verbose
  )

  n_valid <- sum(is.finite(surface$fhat))
  if (n_valid < 10L)
    stop("Only ", n_valid, " grid points have finite fhat; ",
         "need at least 10. Try increasing grid_points or sample_size.")

  if (ctrl$verbose)
    cat(sprintf("  GAM surface: %d/%d valid grid points\n",
                n_valid, nrow(surface)))

  # Fit GAM to the surface
  gam_fit <- train_likelihood_GAM(
    surface, par_names = par_names,
    spline_type = ctrl$spline_type
  )

  # Optimize to find MLE
  mle <- find_MLE(gam_fit,
                  lower_bound = compact_lb,
                  upper_bound = compact_ub,
                  par_names   = par_names)

  pars8 <- .expand_pars(mle$pars, model_bin)

  details <- list(
    surface     = surface,
    gam_fit     = gam_fit,
    mle_optim   = mle$optim,
    convergence = mle$convergence,
    par_names   = par_names,
    grid_points = gp
  )

  list(pars = pars8, loglik = mle$loglik, loglik_var = NA_real_,
       details = details)
}

# --------------------------------------------------------------------------- #
#  Public API                                                                  #
# --------------------------------------------------------------------------- #

#' Estimate diversification rates from a phylogenetic tree
#'
#' Unified entry point for parameter inference under the general
#' diversification model. Accepts a simulated or empirical tree, extracts
#' branching times (pruning to extant tips if needed), and dispatches to
#' the chosen optimisation back-end.
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
#'     \item{\code{"mcem"}}{Monte Carlo Expectation-Maximisation (default).
#'       Iterates E- and M-steps until parameters stabilise (max relative
#'       change < \code{tol} for \code{patience} consecutive iterations)
#'       or \code{max_iter} is reached. Requires \code{init_pars}.}
#'     \item{\code{"cem"}}{Monte Carlo Cross-Entropy Method. Initialises a
#'       particle population randomly within the bounds and evolves it until
#'       convergence (annealing exhausted, plateau, or \code{control$max_iter}
#'       reached). Does not need \code{init_pars}.}
#'     \item{\code{"gam"}}{GAM-based MLE. Evaluates the IS log-likelihood
#'       on a grid spanning the bounds, fits a smooth GAM surface, and
#'       optimises it via L-BFGS-B. One-shot (non-iterative); no
#'       \code{init_pars} needed. Requires the \pkg{mgcv} package.}
#'   }
#' @param model Model specification. Accepts:
#'   \itemize{
#'     \item a formula such as \code{~ N}, \code{~ N + PD}, \code{~ PD + EP},
#'     \item a string shortcut (\code{"cr"}, \code{"dd"}, \code{"pd"},
#'       \code{"ep"}), or
#'     \item a length-3 binary integer vector \code{c(use_N, use_P, use_E)}.
#'   }
#'   Default \code{"cr"} (constant rate).
#' @param init_pars Starting parameter vector. Required for \code{"mcem"};
#'   ignored for \code{"cem"}. If \code{NULL} with \code{"mcem"}, the
#'   midpoint of the bounds is used.
#' @param control Named list of tuning parameters. Must include
#'   \code{lower_bound} and \code{upper_bound}. Use
#'   \code{\link{estimate_rates_control}} to inspect other defaults.
#'   Key entries:
#'   \describe{
#'     \item{\code{lower_bound}}{Numeric vector of lower parameter bounds.
#'       Length must be \code{2 + 2*sum(model)}.
#'       Order: \code{c(beta_0, [beta_N], [beta_P], [beta_E],
#'       gamma_0, [gamma_N], [gamma_P], [gamma_E])}.}
#'     \item{\code{upper_bound}}{Numeric vector of upper parameter bounds
#'       (same order).}
#'     \item{\code{verbose}}{Print progress messages. Default \code{FALSE}.}
#'   }
#' @param link Link function for rate computation: \code{"linear"} (default)
#'   uses \code{max(0, eta)}; \code{"exponential"} uses \code{exp(eta)}.
#' @return A list with class \code{"emphasis_fit"} containing:
#'   \describe{
#'     \item{\code{pars}}{Named numeric vector of parameter estimates.}
#'     \item{\code{loglik}}{Final log-likelihood estimate.}
#'     \item{\code{loglik_var}}{Bootstrap variance of \code{loglik} due to
#'       Monte Carlo noise. Populated automatically when
#'       \code{control$num_trees >= 2} (B = 200 replicates); \code{NA}
#'       when \code{num_trees = 1}.}
#'     \item{\code{n_pars}}{Number of free parameters (used for AIC).}
#'     \item{\code{AIC}}{Akaike Information Criterion: \code{-2 * loglik + 2 * n_pars}.}
#'     \item{\code{method}}{Method used (\code{"mcem"} or \code{"cem"}).}
#'     \item{\code{model}}{Resolved binary model vector.}
#'     \item{\code{details}}{Full output from the back-end function
#'       (\code{.mcem_dynamic_fresh} or \code{emphasis_cem}).}
#'   }
#' @examples
#' \dontrun{
#' set.seed(42)
#' # Constant-rate model
#' sim <- simulate_tree(c(0.5, 0.1), max_t = 5, model = "cr")
#' fit <- estimate_rates(sim, method = "mcem", model = "cr",
#'   control = list(lower_bound = c(0, 0), upper_bound = c(2, 1)))
#' fit$pars
#'
#' # Diversity-dependent model
#' sim_dd <- simulate_tree(c(0.5, -0.005, 0.1, 0), max_t = 8, model = "dd")
#' fit_dd <- estimate_rates(sim_dd, method = "mcem", model = "dd",
#'   control = list(lower_bound = c(0.1, -0.1, 0, -0.01),
#'                  upper_bound = c(2, 0.01, 0.5, 0.01)))
#' fit_dd$pars
#' }
#' @export
estimate_rates <- function(tree,
                           method    = c("mcem", "cem", "gam"),
                           model     = "cr",
                           init_pars = NULL,
                           control   = list(),
                           link      = "linear") {
  method    <- match.arg(method)

  model_bin <- .resolve_model(model)
  link_int  <- .resolve_link(link)

  defaults <- estimate_rates_control(method)
  ctrl     <- utils::modifyList(defaults, control)
  # Resolve parameter aliases (old names still work)
  ctrl <- .resolve_control_aliases(ctrl, method, user_ctrl = control)

  lower_bound <- ctrl$lower_bound
  upper_bound <- ctrl$upper_bound

  if (is.null(lower_bound) || is.null(upper_bound))
    stop("'lower_bound' and 'upper_bound' must be supplied in 'control'.\n",
         "  e.g. control = list(lower_bound = c(...), upper_bound = c(...))")
  if (!is.numeric(lower_bound) || !is.numeric(upper_bound))
    stop("'lower_bound' and 'upper_bound' must be numeric vectors.")
  if (length(lower_bound) != length(upper_bound))
    stop("'lower_bound' and 'upper_bound' must have the same length.")

  expected_n <- 2L + 2L * sum(model_bin)
  if (length(lower_bound) != expected_n)
    stop(.pars_error_msg(model_bin, expected_n))

  brts <- .extract_brts(tree)

  if (is.null(init_pars) && method == "mcem") {
    init_pars <- (lower_bound + upper_bound) / 2
    # For linear link with DD/PD/EP models, the midpoint of the slope bounds
    # can push lambda to zero (lambda = max(0, beta_0 + beta_X * X)).
    # Start slopes near zero so lambda stays positive during early E-steps.
    if (link_int == 0L && any(model_bin != 0L)) {
      n_pars <- length(init_pars)
      slope_idx <- which(model_bin == 1L)           # which covariates are active
      lam_slopes <- 1L + slope_idx                   # compact positions for beta_X
      mu_slopes  <- (n_pars / 2L) + 1L + slope_idx  # compact positions for gamma_X
      all_slopes <- c(lam_slopes, mu_slopes)
      all_slopes <- all_slopes[all_slopes <= n_pars]
      # Clamp slopes: start at 0 if 0 is within bounds, else nearest bound
      for (j in all_slopes) {
        init_pars[j] <- max(lower_bound[j], min(upper_bound[j], 0))
      }
    }
  }

  # Warn if PD/EP model with linear link — IS collapse is very likely
  if (link_int == 0L && (model_bin[2] == 1L || model_bin[3] == 1L) &&
      any(lower_bound < 0)) {
    model_str <- .model_label(model_bin)
    message(
      "Note: ", model_str, " model with linear link can suffer IS collapse ",
      "when the covariate drives lambda to zero.\n",
      "  Consider using link=\"exponential\" for more stable estimation."
    )
  }

  # Expand compact pars/bounds to full 8-element layout for C++
  lb8 <- .expand_pars(lower_bound, model_bin)
  ub8 <- .expand_pars(upper_bound, model_bin)
  ip8 <- if (!is.null(init_pars)) .expand_pars(init_pars, model_bin) else NULL

  if (ctrl$verbose) {
    model_str <- .model_label(model_bin)
    cat(sprintf(
      "[estimate_rates] method=%s  model=%s  n_brts=%d  link=%s\n",
      toupper(method), model_str, length(brts),
      if (link_int == 0L) "linear" else "exponential"
    ))
    if (method == "cem") {
      ss   <- ctrl$num_trees
      np   <- ctrl$num_particles
      mode <- if (isTRUE(ctrl$shared_trees)) "shared (Mode 2)" else "independent (Mode 1)"
      if (ss == 1L) {
        fhat_desc <- paste0(
          "  fhat           : logf - log_q  (single IS log-weight; high variance)\n",
          "                   loglik_var = NA  -- use num_trees >= 5 for variance\n")
      } else {
        fhat_desc <- sprintf(paste0(
          "  fhat           : log mean exp(logf - log_q)  over %d trees\n",
          "                   loglik_var bootstrapped automatically (B=200)\n"), ss)
        if (ss < 5L)
          fhat_desc <- paste0(fhat_desc,
            "                   NOTE: num_trees < 5; recommend >= 5 for stable variance\n")
      }
      cat(sprintf(paste0(
        "  particles/iter : %d\n",
        "  trees/particle : %d  =>  ~%d augmented trees per iteration\n"),
        np, ss, np * ss))
      cat(fhat_desc)
      cat(sprintf("  tree mode      : %s\n  max_iter=%d\n",
                  mode, ctrl$max_iter))
    }
    if (method == "mcem") {
      cat(sprintf("  num_trees=%d  max_iter=%d  tol=%.1e  patience=%d  xtol=%.1e\n",
                  ctrl$num_trees, ctrl$max_iter, ctrl$tol, ctrl$patience, ctrl$xtol))
      if (ctrl$num_trees >= 2L)
        cat("  loglik_var   : bootstrapped automatically (B=200)\n")
      else
        cat("  loglik_var   : NA (need num_trees >= 2)\n")
    }
    if (method == "gam") {
      cat(sprintf("  grid_points=%d  sample_size=%d  spline=%s\n",
                  ctrl$grid_points, ctrl$sample_size, ctrl$spline_type))
    }
  }

  raw <- switch(method,
    mcem = .run_mcem(brts, ip8, lb8, ub8, ctrl, model = model_bin, link = link_int),
    cem  = .run_cem(brts, lb8, ub8, ctrl, model = model_bin, link = link_int),
    gam  = .run_gam(brts, lb8, ub8, ctrl, model = model_bin, link = link_int)
  )

  # Contract 8-element result back to compact
  compact_pars <- .contract_pars(as.numeric(raw$pars), model_bin)
  pars <- stats::setNames(compact_pars, .par_names(model_bin))
  n_pars <- length(pars)
  aic <- if (is.finite(raw$loglik)) -2 * raw$loglik + 2 * n_pars else NA_real_
  loglik_var <- if (!is.null(raw$loglik_var)) raw$loglik_var else NA_real_

  if (ctrl$verbose) {
    cat(sprintf("[estimate_rates] loglik=%.4f  AIC=%.4f  pars: %s\n",
                raw$loglik, aic,
                paste(names(pars), "=", round(pars, 4), collapse="  ")))
  }

  result <- list(pars = pars, loglik = raw$loglik, loglik_var = loglik_var,
                 n_pars = n_pars, AIC = aic,
                 method = method, model = model_bin, details = raw$details)
  class(result) <- "emphasis_fit"
  result
}


#' Print method for emphasis_fit objects
#'
#' @param x An \code{emphasis_fit} object returned by \code{\link{estimate_rates}}.
#' @param ... Additional arguments (ignored).
#' @export
print.emphasis_fit <- function(x, ...) {
  model_str <- .model_label(x$model)
  cat("emphasis fit (", x$method, ", model = ", model_str, ")\n\n", sep = "")
  cat("Parameters:\n")
  print(round(x$pars, 6))
  cat("\nLog-likelihood:", round(x$loglik, 4))
  if (!is.na(x$loglik_var))
    cat(sprintf("  (MC se: %.4f)", sqrt(x$loglik_var)))
  cat("\nAIC:           ", round(x$AIC, 4), "\n")
  cat("n_pars:        ", x$n_pars, "\n")
  invisible(x)
}


#' Compare fitted emphasis models via AIC
#'
#' Given two or more \code{\link{estimate_rates}} results, returns a summary
#' table sorted by AIC. The model with the lowest AIC is preferred.
#'
#' When all fits carry a \code{loglik_var} (bootstrap variance, available
#' when \code{control$num_trees > 1}), pairwise Gaussian tests of equal AIC
#' are appended.  For each pair \eqn{(i, j)},
#' \deqn{T_{ij} = \widehat{\rm AIC}_i - \widehat{\rm AIC}_j, \qquad
#'   \mathrm{Var}(T_{ij}) = 4[\mathrm{Var}(\hat\ell_i) +
#'   \mathrm{Var}(\hat\ell_j)]}
#' and \eqn{p = 2\Phi(-|T_{ij}|/\sqrt{\mathrm{Var}(T_{ij})})} under
#' the null of equal AIC.
#'
#' @param ... Two or more \code{emphasis_fit} objects (results of
#'   \code{\link{estimate_rates}}). Optionally named; unnamed fits are
#'   auto-labelled from their model specification.
#' @return A data frame with columns \code{model}, \code{n_pars},
#'   \code{loglik}, \code{AIC}, \code{delta_AIC}, and \code{AICw}.
#'   If bootstrap variances are present, an additional \code{loglik_se}
#'   column is included and pairwise p-values are printed as a message.
#' @examples
#' \dontrun{
#' fit_cr <- estimate_rates(tree, model = "cr",
#'   lower_bound = c(0, 0), upper_bound = c(2, 1))
#' fit_dd <- estimate_rates(tree, model = "dd",
#'   lower_bound = c(0, -0.1, 0, -0.01),
#'   upper_bound = c(2, 0.01, 0.5, 0.01))
#' compare_models(CR = fit_cr, DD = fit_dd)
#' }
#' @keywords internal
compare_models <- function(...) {
  fits <- list(...)
  if (length(fits) < 2L)
    stop("compare_models() requires at least two fitted models.")

  nms <- names(fits)
  if (is.null(nms)) nms <- rep("", length(fits))
  for (i in seq_along(fits)) {
    if (nms[i] == "") nms[i] <- .model_label(fits[[i]]$model)
  }

  get_var <- function(f) if (!is.null(f$loglik_var)) f$loglik_var else NA_real_
  vars <- vapply(fits, get_var, numeric(1L))

  tab <- data.frame(
    model     = nms,
    n_pars    = vapply(fits, `[[`, 0L, "n_pars"),
    loglik    = vapply(fits, `[[`, 0.0, "loglik"),
    AIC       = vapply(fits, `[[`, 0.0, "AIC"),
    stringsAsFactors = FALSE
  )

  # Include MC standard errors if available
  if (any(!is.na(vars))) {
    tab$loglik_se <- sqrt(vars)
  }

  tab <- tab[order(tab$AIC, na.last = TRUE), ]
  tab$delta_AIC <- tab$AIC - min(tab$AIC, na.rm = TRUE)
  w <- exp(-0.5 * tab$delta_AIC)
  tab$AICw <- w / sum(w, na.rm = TRUE)
  rownames(tab) <- NULL

  # Pairwise AIC Gaussian test when all fits have bootstrap variance
  if (all(!is.na(vars)) && length(fits) >= 2L) {
    n <- length(fits)
    pmat <- matrix(NA_real_, n, n, dimnames = list(nms, nms))
    for (i in seq_len(n - 1L)) {
      for (j in seq(i + 1L, n)) {
        T_ij   <- tab$AIC[tab$model == nms[i]] - tab$AIC[tab$model == nms[j]]
        var_T  <- 4 * (vars[i] + vars[j])
        if (is.finite(var_T) && var_T > 0) {
          pmat[i, j] <- pmat[j, i] <-
            2 * stats::pnorm(-abs(T_ij) / sqrt(var_T))
        }
      }
    }
    attr(tab, "pairwise_p") <- pmat
  }

  tab
}


#' @keywords internal
.model_label <- function(model_bin) {
  covs <- c("N", "PD", "EP")[which(model_bin == 1L)]
  if (length(covs) == 0L) return("CR")
  paste(covs, collapse = " + ")
}


# --------------------------------------------------------------------------- #
#  Three-model covariate comparison                                            #
# --------------------------------------------------------------------------- #

#' Fit and compare the three single-covariate diversification models
#'
#' Fits three 4-parameter models -- diversity dependence (DD, covariate N),
#' phylogenetic diversity dependence (PD, covariate P), and evolutionary
#' pendant (EP, covariate E) -- to the same tree using a two-stage
#' MCDE -> MCEM workflow, then ranks them by AIC and computes AIC weights.
#'
#' All three models share the same parameter layout:
#' \code{c(beta_0, beta_X, gamma_0, gamma_X)}, so a single pair of bounds
#' applies to all three.
#'
#' @param tree Ultrametric tree. Accepts a \code{phylo} object, a
#'   \code{\link{simulate_tree}} result, or a numeric branching-time vector.
#' @param lower_bound Length-4 numeric vector of lower parameter bounds:
#'   \code{c(beta_0, beta_X, gamma_0, gamma_X)}. Default
#'   \code{c(0.1, -0.5, 0, -0.5)}.
#' @param upper_bound Length-4 numeric vector of upper parameter bounds.
#'   Default \code{c(3, 0.5, 1, 0.5)}.
#' @param link Link function: \code{"linear"} (default) or
#'   \code{"exponential"} (EP not supported with exponential).
#' @param control Named list with optional sub-lists \code{cem} and
#'   \code{mcem} to override tuning defaults for each stage. Bounds and
#'   verbose are taken from the outer arguments; no need to repeat them here.
#'   Example: \code{list(cem = list(max_iter = 20, num_particles = 50),
#'                       mcem = list(num_trees = 300, tol = 0.05))}.
#' @param verbose Print progress messages. Default \code{FALSE}.
#' @return A list of class \code{"model_selection"} with:
#'   \describe{
#'     \item{\code{summary}}{Data frame with columns \code{model},
#'       \code{n_pars}, \code{loglik}, \code{AIC}, \code{delta_AIC},
#'       \code{AICw}, sorted by AIC.}
#'     \item{\code{best_model}}{Name of the best-supported model
#'       (\code{"dd"}, \code{"pd"}, or \code{"ep"}).}
#'     \item{\code{best_pars}}{Named parameter estimates of the best model.}
#'     \item{\code{fits}}{Named list of the three \code{emphasis_fit}
#'       objects (\code{$dd}, \code{$pd}, \code{$ep}).}
#'   }
#' @seealso \code{\link{estimate_rates}}, \code{\link{compare_models}}
#' @examples
#' \dontrun{
#' set.seed(42)
#' sim <- simulate_tree(pars = c(0.8, -0.02, 0.2, 0), model = "dd", max_t = 8)
#' sel <- select_diversification_model(sim,
#'   lower_bound = c(0.1, -0.5, 0, -0.5),
#'   upper_bound = c(3,    0.5, 1,  0.5))
#' sel
#' sel$best_model   # "dd"
#' sel$best_pars
#' }
#' @keywords internal
select_diversification_model <- function(tree,
                                         lower_bound = c(0.1, -0.5, 0, -0.5),
                                         upper_bound = c(3,    0.5, 1,  0.5),
                                         link        = "linear",
                                         control     = list(),
                                         verbose     = FALSE) {
  if (length(lower_bound) != 4L || length(upper_bound) != 4L)
    stop("'lower_bound' and 'upper_bound' must each have length 4: ",
         "c(beta_0, beta_X, gamma_0, gamma_X).")

  defaults <- list(
    cem  = list(max_iter = 15, num_particles = 40),
    mcem = list(num_trees = 200, tol = 1e-3)
  )
  ctrl <- utils::modifyList(defaults, control)

  # Inject bounds and verbose into each stage's control list
  bounds <- list(lower_bound = lower_bound, upper_bound = upper_bound,
                 verbose = verbose)
  ctrl$cem  <- utils::modifyList(bounds, ctrl$cem)
  ctrl$mcem <- utils::modifyList(bounds, ctrl$mcem)

  model_names <- c("dd", "pd", "ep")
  fits <- vector("list", 3L)
  names(fits) <- model_names

  for (m in model_names) {
    if (verbose) message("Fitting model: ", toupper(m), " ...")

    fit_de <- tryCatch(
      estimate_rates(tree, method = "cem", model = m,
                     link = link, control = ctrl$cem),
      error = function(e) {
        warning("MCDE failed for model '", m, "': ", conditionMessage(e))
        NULL
      }
    )

    init <- if (!is.null(fit_de)) fit_de$pars else (lower_bound + upper_bound) / 2

    fits[[m]] <- tryCatch(
      estimate_rates(tree, method = "mcem", model = m,
                     init_pars = init, link = link, control = ctrl$mcem),
      error = function(e) {
        warning("MCEM failed for model '", m, "': ", conditionMessage(e))
        NULL
      }
    )
  }

  # Drop any models that failed entirely
  fits <- Filter(Negate(is.null), fits)
  if (length(fits) < 2L)
    stop("Fewer than two models converged; cannot perform model selection.")

  tab <- do.call(compare_models, fits)

  # Best model = row with delta_AIC == 0
  best_key <- tab$model[tab$delta_AIC == 0 & !is.na(tab$delta_AIC)][1L]
  if (is.na(best_key) || !best_key %in% names(fits))
    best_key <- names(fits)[1L]

  result <- list(
    summary    = tab,
    best_model = unname(best_key),
    best_pars  = fits[[best_key]]$pars,
    fits       = fits
  )
  class(result) <- "model_selection"
  result
}


#' Print method for model_selection objects
#'
#' @param x A \code{model_selection} object from
#'   \code{\link{select_diversification_model}}.
#' @param ... Additional arguments (ignored).
#' @export
print.model_selection <- function(x, ...) {
  cat("Diversification model selection\n")
  cat("================================\n\n")
  tab <- x$summary
  tab$loglik    <- round(tab$loglik,    3)
  tab$AIC       <- round(tab$AIC,       3)
  tab$delta_AIC <- round(tab$delta_AIC, 3)
  tab$AICw      <- round(tab$AICw,      4)
  print(tab, row.names = FALSE)
  cat("\nBest model:", toupper(x$best_model), "\n")
  cat("Parameters:\n")
  print(round(x$best_pars, 6))
  invisible(x)
}

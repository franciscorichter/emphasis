#' Train a survival-probability GAM on forward simulation results
#'
#' Fits a binomial GAM to estimate P(survival | parameters) from forward
#' simulations run with \code{max_tries = 0}.  Each simulation row has a
#' binary outcome: 1 if \code{status == "done"}, 0 if \code{status ==
#' "extinct"} or \code{"too_large"}.
#'
#' @param simulations List of \code{simulate_tree()} outputs (one per row of
#'   \code{pars_mat}), produced with \code{max_tries = 0}.
#' @param pars_mat Numeric matrix whose rows are the parameter vectors used to
#'   generate \code{simulations}.  Column names are used as predictor names; if
#'   absent, names are auto-generated from \code{model}.
#' @param model String or length-3 binary integer vector matching the model
#'   used in simulation.  Used only when \code{pars_mat} has no column names.
#'   Default \code{"cr"}.
#' @param spline_type \code{"univariate"} (default) fits independent smooth
#'   terms \code{s(x)} for each varying predictor.  \code{"bivariate"} fits
#'   tensor-product smooths \code{te(x, y)} for every pair of varying
#'   predictors, which can capture interaction effects.
#' @return A \code{gam} object (binomial family) representing
#'   \eqn{P(\text{survival} \mid \theta)}.
#' @seealso \code{\link{simulate_tree}}, \code{\link{predict_survival}}
#' @export
#' @examples
#' \dontrun{
#' set.seed(1)
#' pars_mat <- cbind(beta_0  = runif(200, 0.3, 1.2),
#'                   gamma_0 = runif(200, 0.05, 0.5))
#' sims <- simulate_tree(pars = pars_mat, max_t = 8,
#'                       model = "cr", max_tries = 0, useDDD = FALSE)
#' gam_uni <- train_GAM(sims, pars_mat, model = "cr")
#' gam_biv <- train_GAM(sims, pars_mat, model = "cr", spline_type = "bivariate")
#' }
train_GAM <- function(simulations, pars_mat, model = "cr",
                      spline_type = c("univariate", "bivariate")) {
  spline_type <- match.arg(spline_type)
  survived <- sapply(simulations, function(s) as.integer(s$status == "done"))

  pars_df <- as.data.frame(pars_mat)
  if (is.null(colnames(pars_mat)) || any(colnames(pars_mat) == "")) {
    colnames(pars_df) <- .par_names(.resolve_model(model))
  }

  sim_data  <- cbind(pars_df, survived = survived)
  par_names <- colnames(pars_df)

  # Drop constant columns -- mgcv cannot fit a smooth on zero-variance predictors
  varying <- par_names[sapply(par_names,
                              function(p) length(unique(pars_df[[p]])) > 1L)]
  if (length(varying) == 0L) stop("All parameter columns are constant; nothing to fit.")

  if (spline_type == "univariate" || length(varying) < 2L) {
    smooth_terms <- paste0("s(", varying, ")", collapse = " + ")
  } else {
    # Bivariate: tensor-product smooth for every pair of varying predictors
    pairs        <- utils::combn(varying, 2L, simplify = FALSE)
    smooth_terms <- paste(
      vapply(pairs, function(p) paste0("te(", p[1L], ", ", p[2L], ")"),
             character(1L)),
      collapse = " + "
    )
  }

  form <- stats::as.formula(paste("survived ~", smooth_terms))
  cat("Training", spline_type, "survival GAM on", nrow(sim_data),
      "simulations (survival rate:", round(mean(survived), 3), ")...\n")
  mgcv::gam(form, data = sim_data, family = stats::binomial())
}


#' Predict survival probability from a trained GAM
#'
#' Convenience wrapper around \code{predict.gam(..., type = "response")}.
#'
#' @param gam_fit A \code{gam} object returned by \code{\link{train_GAM}}.
#' @param newpars Numeric matrix or data frame of parameter values to predict
#'   at.  Column names must match those used in \code{train_GAM}.
#' @return Numeric vector of predicted survival probabilities.
#' @export
predict_survival <- function(gam_fit, newpars) {
  stats::predict(gam_fit, newdata = as.data.frame(newpars), type = "response")
}


# =========================================================================== #
#  Automatic bound detection                                                    #
# =========================================================================== #

#' Automatically determine parameter bounds via forward simulation
#'
#' Starts from a safe center point (low rates, zero covariates) and expands
#' outward along each parameter axis using bisection to find the feasibility
#' boundary --- where trees transition from surviving with sensible size to
#' going extinct or exploding.  Much more efficient than blanket LHS sampling
#' over a wide box.  Optionally trains a survival GAM over the detected region.
#'
#' @param tree A \code{phylo} object or anything accepted by
#'   \code{\link{estimate_rates}}.
#' @param model Model specification: \code{"cr"}, \code{"dd"}, etc.
#' @param link \code{"linear"} or \code{"exponential"}.
#' @param n_test Number of replicate simulations per candidate point
#'   to assess feasibility (default 5).
#' @param bisect_steps Number of bisection steps per axis direction
#'   (default 8; gives ~1/256 resolution).
#' @param margin Fractional expansion of detected bounds (default 0.5).
#' @param tip_range Multiplier range for acceptable tip counts relative
#'   to the observed tree.  Default \code{c(0.1, 10)}.
#' @param train_surv_gam If \code{TRUE} (default), train a survival GAM
#'   on the detected region.
#' @param num_threads Parallel threads for batch simulation (default 1).
#' @param verbose Print progress (default \code{TRUE}).
#' @return A list with components:
#'   \describe{
#'     \item{lower_bound}{Numeric vector of lower bounds (compact layout).}
#'     \item{upper_bound}{Numeric vector of upper bounds (compact layout).}
#'     \item{survival_gam}{A \code{gam} object if trained, else \code{NULL}.}
#'     \item{center}{The safe center point used as starting location.}
#'   }
#' @export
#' @examples
#' \dontrun{
#' library(ape)
#' data(bird.orders)
#' ab <- auto_bounds(bird.orders, model = "dd", link = "exponential")
#' fit <- estimate_rates(bird.orders, method = "gam", model = "dd",
#'          link = "exponential",
#'          control = list(lower_bound = ab$lower_bound,
#'                         upper_bound = ab$upper_bound),
#'          cond = ab$survival_gam)
#' }
auto_bounds <- function(tree, model = "cr", link = "linear",
                        n_test = 5L, bisect_steps = 8L,
                        margin = 0.5,
                        tip_range = c(0.1, 10),
                        train_surv_gam = TRUE,
                        num_threads = 1L,
                        verbose = TRUE) {
  model_bin <- .resolve_model(model)
  link_int  <- .resolve_link(link)
  n_pars    <- 2L + 2L * sum(model_bin)
  pnames    <- .par_names(model_bin)

  # Crown age and target tips from observed tree
  if (inherits(tree, "phylo")) {
    max_t  <- max(ape::branching.times(tree))
    n_tips <- ape::Ntip(tree)
  } else if (is.numeric(tree)) {
    max_t  <- max(tree)
    n_tips <- length(tree) + 1L
  } else if (is.list(tree) && !is.null(tree$tes)) {
    max_t  <- max(ape::branching.times(tree$tes))
    n_tips <- ape::Ntip(tree$tes)
  } else {
    stop("'tree' must be a phylo, branching times, or simulate_tree result.")
  }

  brts    <- .extract_brts(tree)
  tip_lo  <- max(2, floor(n_tips * tip_range[1]))
  tip_hi  <- ceiling(n_tips * tip_range[2])
  max_lin <- as.integer(max(20 * n_tips, 500))

  # ── Safe center: observed net rate, zero covariates ────────
  r_hat  <- log(max(n_tips, 2) / 2) / max(max_t, 0.1)
  mu_hat <- 0.2 * r_hat
  lam_hat <- r_hat + mu_hat

  if (link_int == 0L) {
    center_lam <- c(lam_hat)
    center_mu  <- c(mu_hat)
  } else {
    center_lam <- c(log(max(lam_hat, 1e-4)))
    center_mu  <- c(log(max(mu_hat, 1e-6)))
  }
  active <- which(model_bin == 1L)
  center <- c(center_lam, rep(0, length(active)),
              center_mu,  rep(0, length(active)))
  names(center) <- pnames

  # Verify center is feasible
  if (verbose) cat(sprintf(
    "auto_bounds: %s link=%s (%d pars)\n",
    paste(model_bin, collapse = ","),
    if (link_int == 0L) "linear" else "exponential",
    n_pars
  ))
  if (verbose) cat(sprintf(
    "  center: %s\n",
    paste(round(center, 4), collapse = ", ")
  ))

  center_ok <- .test_feasibility(
    center, model, link, max_t, max_lin,
    n_test = 10L, tip_lo = tip_lo, tip_hi = tip_hi,
    num_threads = num_threads
  )

  if (!center_ok) {
    # Try a grid of centers near r_hat
    if (verbose) cat("  center infeasible, searching...\n")
    center <- .find_feasible_center(
      model_bin, link_int, max_t, n_tips, model, link,
      max_lin, tip_lo, tip_hi, num_threads
    )
    if (is.null(center)) {
      warning("auto_bounds: could not find a feasible center point.")
      wide <- .wide_bounds(model_bin, link_int, max_t, n_tips)
      return(list(
        lower_bound = wide$lb, upper_bound = wide$ub,
        survival_gam = NULL, center = NULL
      ))
    }
    names(center) <- pnames
    if (verbose) cat(sprintf(
      "  found center: %s\n",
      paste(round(center, 4), collapse = ", ")
    ))
  } else {
    if (verbose) cat("  center feasible\n")
  }

  # -- Collect all feasible points to determine bounds ---------
  wide <- .wide_bounds(model_bin, link_int, max_t, n_tips)
  feasible_pts <- list(center)  # center is known feasible

  # Phase 1: Axis-aligned bisection (original method)
  if (verbose) cat("  Phase 1: axis-aligned bisection\n")
  for (j in seq_len(n_pars)) {
    if (verbose) cat(sprintf("    axis %s: ", pnames[j]))
    hi <- .bisect_boundary(center, j, center[j], wide$ub[j],
                           model, link, max_t, max_lin,
                           n_test, bisect_steps, tip_lo, tip_hi, num_threads)
    lo <- .bisect_boundary(center, j, center[j], wide$lb[j],
                           model, link, max_t, max_lin,
                           n_test, bisect_steps, tip_lo, tip_hi, num_threads)
    if (verbose) cat(sprintf("[%.4f, %.4f]\n", lo, hi))
    # Record the feasible boundary points
    pt_hi <- center; pt_hi[j] <- hi; feasible_pts[[length(feasible_pts) + 1L]] <- pt_hi
    pt_lo <- center; pt_lo[j] <- lo; feasible_pts[[length(feasible_pts) + 1L]] <- pt_lo
  }

  # Phase 2: Compensatory diagonals
  # For each active covariate, explore the direction where the intercept
  # and the covariate coefficient compensate each other, keeping the
  # predicted rate at the observed covariate value roughly constant.
  # This is general: eta = intercept + coeff * X_obs = const means
  # d(intercept) = -X_obs * d(coeff), i.e. direction (1, -X_obs) in
  # the (intercept, coeff) subspace.
  obs_covs <- .observed_covariates(brts, model_bin)

  if (length(obs_covs) > 0L) {
    if (verbose) cat("  Phase 2: compensatory diagonals\n")

    for (cov_info in obs_covs) {
      i_int  <- cov_info$intercept_idx
      i_coef <- cov_info$coeff_idx
      X_obs  <- cov_info$X_obs

      if (verbose) cat(sprintf("    %s (X_obs=%.2f):\n", cov_info$name, X_obs))

      # Explore multiple compensatory slopes:
      # The exact compensatory slope is X_obs, but we also explore
      # shallower slopes (X_obs/2, X_obs/5, 1) to push the intercept further.
      # Each slope s means: direction = (1, -s) in the (intercept, coeff) subspace.
      slopes <- unique(c(X_obs, X_obs / 2, X_obs / 5, X_obs / 20, 1, 0.1))

      for (s in slopes) {
        direction <- rep(0, n_pars)
        direction[i_int]  <- 1
        direction[i_coef] <- -s
        direction <- direction / sqrt(sum(direction^2))

        hi_pt <- .bisect_direction(center, direction, wide, 1,
                                   model, link, max_t, max_lin,
                                   n_test, bisect_steps, tip_lo, tip_hi, num_threads)
        lo_pt <- .bisect_direction(center, direction, wide, -1,
                                   model, link, max_t, max_lin,
                                   n_test, bisect_steps, tip_lo, tip_hi, num_threads)

        if (!is.null(hi_pt)) feasible_pts[[length(feasible_pts) + 1L]] <- hi_pt
        if (!is.null(lo_pt)) feasible_pts[[length(feasible_pts) + 1L]] <- lo_pt

        if (verbose) {
          hi_str <- if (!is.null(hi_pt)) sprintf("(%.4f, %.4f)", hi_pt[i_int], hi_pt[i_coef]) else "---"
          lo_str <- if (!is.null(lo_pt)) sprintf("(%.4f, %.4f)", lo_pt[i_int], lo_pt[i_coef]) else "---"
          cat(sprintf("      slope=%.1f: %s .. %s\n", s, lo_str, hi_str))
        }
      }

      # Same for the extinction side
      i_int_mu  <- cov_info$mu_intercept_idx
      i_coef_mu <- cov_info$mu_coeff_idx

      for (s in slopes) {
        dir_mu <- rep(0, n_pars)
        dir_mu[i_int_mu]  <- 1
        dir_mu[i_coef_mu] <- -s
        dir_mu <- dir_mu / sqrt(sum(dir_mu^2))

        hi_mu <- .bisect_direction(center, dir_mu, wide, 1,
                                   model, link, max_t, max_lin,
                                   n_test, bisect_steps, tip_lo, tip_hi, num_threads)
        lo_mu <- .bisect_direction(center, dir_mu, wide, -1,
                                   model, link, max_t, max_lin,
                                   n_test, bisect_steps, tip_lo, tip_hi, num_threads)
        if (!is.null(hi_mu)) feasible_pts[[length(feasible_pts) + 1L]] <- hi_mu
        if (!is.null(lo_mu)) feasible_pts[[length(feasible_pts) + 1L]] <- lo_mu
      }
    }
  }

  # Phase 3: Random directions
  n_random <- max(4L, n_pars)
  if (verbose) cat(sprintf("  Phase 3: %d random directions\n", n_random))
  set.seed(12345L)  # reproducible
  for (k in seq_len(n_random)) {
    direction <- stats::rnorm(n_pars)
    direction <- direction / sqrt(sum(direction^2))
    hi_pt <- .bisect_direction(center, direction, wide, 1,
                               model, link, max_t, max_lin,
                               n_test, bisect_steps, tip_lo, tip_hi, num_threads)
    lo_pt <- .bisect_direction(center, direction, wide, -1,
                               model, link, max_t, max_lin,
                               n_test, bisect_steps, tip_lo, tip_hi, num_threads)
    if (!is.null(hi_pt)) feasible_pts[[length(feasible_pts) + 1L]] <- hi_pt
    if (!is.null(lo_pt)) feasible_pts[[length(feasible_pts) + 1L]] <- lo_pt
  }

  # -- Compute bounds as component-wise min/max of all feasible points --
  feas_mat <- do.call(rbind, feasible_pts)
  lb <- apply(feas_mat, 2, min)
  ub <- apply(feas_mat, 2, max)

  # Add margin
  span <- ub - lb
  span <- pmax(span, 1e-6)
  lb <- pmax(lb - margin * span, wide$lb)
  ub <- pmin(ub + margin * span, wide$ub)

  if (verbose) {
    cat(sprintf("  %d feasible points collected\n", nrow(feas_mat)))
    cat("  final bounds:\n")
    for (j in seq_len(n_pars)) {
      cat(sprintf("    %s: [%.4f, %.4f]\n", pnames[j], lb[j], ub[j]))
    }
  }

  # -- IS feasibility diagnostic (informational only) ----------
  if (verbose) {
    cat("  IS feasibility check...\n")
    .diagnose_is_bounds(center, lb, ub, brts, model_bin, link_int, verbose)
  }

  # ── Survival GAM on detected region ────────────────────────
  surv_gam <- NULL
  if (train_surv_gam) {
    if (verbose) cat("  training survival GAM...\n")
    n_gam   <- 500L
    gam_mat <- .lhs_sample(n_gam, lb, ub)
    colnames(gam_mat) <- pnames

    gam_sims <- simulate_tree(
      pars = gam_mat, max_t = max_t,
      model = model, link = link,
      max_tries = 0,
      max_lin = max_lin,
      num_threads = num_threads
    )

    surv_gam <- train_GAM(
      gam_sims$simulations, gam_mat, model = model
    )
  }

  list(
    lower_bound  = stats::setNames(lb, pnames),
    upper_bound  = stats::setNames(ub, pnames),
    survival_gam = surv_gam,
    center       = center
  )
}


# Test whether a parameter vector produces feasible trees.
# Returns TRUE if >= 50% of n_test sims survive with sensible tips.
.test_feasibility <- function(pars, model, link, max_t, max_lin,
                              n_test = 5L, tip_lo = 2, tip_hi = 500,
                              num_threads = 1L) {
  pars_mat <- matrix(rep(pars, n_test), nrow = n_test, byrow = TRUE)
  sims <- simulate_tree(
    pars = pars_mat, max_t = max_t,
    model = model, link = link,
    max_tries = 0, max_lin = max_lin,
    num_threads = num_threads
  )
  ntips <- sapply(sims$simulations, function(s) {
    if (s$status == "done" && !is.null(s$tes))
      length(s$tes$tip.label) else 0L
  })
  sum(ntips >= tip_lo & ntips <= tip_hi) >= ceiling(n_test / 2)
}


# Search for a feasible center by trying a grid of rate values.
.find_feasible_center <- function(model_bin, link_int, max_t, n_tips,
                                  model, link, max_lin, tip_lo, tip_hi,
                                  num_threads) {
  r_hat  <- log(max(n_tips, 2) / 2) / max(max_t, 0.1)
  active <- which(model_bin == 1L)

  # Try a grid of baseline rate multipliers
  mults <- c(1, 0.5, 2, 0.25, 3, 0.1, 5)
  mu_fracs <- c(0.2, 0.0, 0.5, 0.1)

  for (m in mults) {
    for (mf in mu_fracs) {
      lam <- r_hat * m / (1 - mf)
      mu  <- lam * mf
      if (link_int == 0L) {
        cand <- c(max(lam, 1e-4), rep(0, length(active)),
                  max(mu, 0), rep(0, length(active)))
      } else {
        cand <- c(log(max(lam, 1e-4)), rep(0, length(active)),
                  log(max(mu, 1e-6)), rep(0, length(active)))
      }
      ok <- .test_feasibility(
        cand, model, link, max_t, max_lin,
        n_test = 5L, tip_lo = tip_lo, tip_hi = tip_hi,
        num_threads = num_threads
      )
      if (ok) return(cand)
    }
  }
  NULL
}


# Bisection to find the boundary of feasibility along one axis.
# Starts from center[j] (feasible) and probes toward `target`.
# Returns the last feasible value found.
.bisect_boundary <- function(center, j, safe, target,
                             model, link, max_t, max_lin,
                             n_test, steps, tip_lo, tip_hi,
                             num_threads) {
  lo <- safe
  hi <- target

  for (s in seq_len(steps)) {
    mid <- (lo + hi) / 2
    cand <- center
    cand[j] <- mid

    ok <- .test_feasibility(
      cand, model, link, max_t, max_lin,
      n_test = n_test, tip_lo = tip_lo, tip_hi = tip_hi,
      num_threads = num_threads
    )

    if (ok) {
      lo <- mid   # feasible — push further
    } else {
      hi <- mid   # infeasible — pull back
    }
  }
  lo  # last known feasible value
}





# Diagnostic: report which boundary corners have IS issues.
# Does NOT modify bounds -- just prints warnings.
.diagnose_is_bounds <- function(center, lb, ub, brts, model_bin, link_int,
                                 verbose = TRUE) {
  n_pars <- length(center)
  any_issues <- FALSE

  for (j in seq_len(n_pars)) {
    for (side in c("lo", "hi")) {
      boundary <- if (side == "lo") lb[j] else ub[j]
      test_pt <- center
      test_pt[j] <- boundary

      pars8 <- .expand_pars(test_pt, model_bin)
      is_ok <- tryCatch({
        raw <- augment_trees(
          brts = brts, pars = pars8,
          sample_size = 3L, maxN = 100L,
          max_missing = 500L, max_lambda = 1e6,
          num_threads = 1L,
          model = as.integer(model_bin),
          link  = as.integer(link_int)
        )
        if (length(raw$trees) == 0L) FALSE
        else {
          logf <- eval_logf(pars8, raw$trees,
                            model = as.integer(model_bin),
                            link  = as.integer(link_int))
          any(is.finite(logf$logf) & is.finite(logf$logg))
        }
      }, error = function(e) FALSE)

      if (!is_ok) {
        if (verbose) cat(sprintf(
          "    IS warning: %s %s bound (%.4f) - augmentation may struggle here\n",
          names(center)[j], side, boundary
        ))
        any_issues <- TRUE
      }
    }
  }

  if (verbose && !any_issues) cat("    all boundary probes OK\n")
  invisible(any_issues)
}


# Compute typical covariate values at the observed tree.
# Returns a list of lists, one per active covariate, each with:
#   intercept_idx, coeff_idx, mu_intercept_idx, mu_coeff_idx, X_obs, name
.observed_covariates <- function(brts, model_bin) {
  n_tips <- length(brts) + 1L
  max_t  <- max(brts)
  active <- which(model_bin == 1L)
  if (length(active) == 0L) return(list())

  # Typical covariate values at the present for the observed tree:
  #   N: number of tips
  #   P: total pendant branch length (sum of pendant ages at present)
  #   E: mean pendant age
  # For pendant ages at the present: each tip's pendant age is
  # (crown_age - most_recent_branching_time_on_its_path).
  # A rough estimate: mean pendant age ~ crown_age / n_tips (for balanced trees)
  # or use the actual branching times.

  # Pendant ages at the present: for each of the n_tips tips,
  # the pendant age is the time from the most recent internal node to the present.
  # The n_tips smallest branching times give the parent nodes of tips.
  sorted_brts <- sort(brts)
  # The youngest n_tips-1 branching times each contribute 2 pendant branches
  # minus the ones that are internal. Approximate: mean pendant age.
  if (length(sorted_brts) >= 1L) {
    pendant_ages <- sorted_brts[seq_len(min(n_tips, length(sorted_brts)))]
    mean_pendant <- mean(pendant_ages)
  } else {
    mean_pendant <- max_t / 2
  }

  result <- list()
  # Parameter layout: [beta_0, (beta_N), (beta_P), (beta_E), gamma_0, ...]
  # Speciation intercept is always index 1
  # Extinction intercept is at 1 + length(active) + 1 = length(active) + 2
  n_lam_pars <- 1L + length(active)
  mu_int_idx <- n_lam_pars + 1L

  for (k in seq_along(active)) {
    cov_type <- active[k]  # 1=N, 2=P, 3=E
    coeff_idx    <- 1L + k           # position in lambda block
    mu_coeff_idx <- mu_int_idx + k   # position in mu block

    X_obs <- switch(cov_type,
      n_tips,                                  # N
      n_tips * mean_pendant,                   # P ~ N * mean_E
      mean_pendant                             # E
    )
    cov_name <- switch(cov_type, "N", "P", "E")

    result[[length(result) + 1L]] <- list(
      intercept_idx    = 1L,
      coeff_idx        = coeff_idx,
      mu_intercept_idx = mu_int_idx,
      mu_coeff_idx     = mu_coeff_idx,
      X_obs            = X_obs,
      name             = cov_name
    )
  }
  result
}


# Bisect along an arbitrary direction in parameter space.
# Starts at center, moves in `sign * direction` until infeasible.
# Returns the last feasible point, or NULL if center itself fails along this direction.
# Clamps all probes to the wide bounds box.
.bisect_direction <- function(center, direction, wide, sign = 1,
                               model, link, max_t, max_lin,
                               n_test, steps, tip_lo, tip_hi,
                               num_threads) {
  n_pars <- length(center)
  # Find maximum step size before hitting any wide bound
  max_step <- Inf
  for (j in seq_len(n_pars)) {
    d <- sign * direction[j]
    if (d > 1e-12) {
      max_step <- min(max_step, (wide$ub[j] - center[j]) / d)
    } else if (d < -1e-12) {
      max_step <- min(max_step, (wide$lb[j] - center[j]) / d)
    }
  }
  if (!is.finite(max_step) || max_step <= 0) return(NULL)

  safe_step <- 0
  far_step  <- max_step

  for (s in seq_len(steps)) {
    mid_step <- (safe_step + far_step) / 2
    cand <- center + sign * mid_step * direction
    # Clamp to wide bounds (safety)
    cand <- pmax(cand, wide$lb)
    cand <- pmin(cand, wide$ub)

    ok <- .test_feasibility(cand, model, link, max_t, max_lin,
                            n_test = n_test, tip_lo = tip_lo, tip_hi = tip_hi,
                            num_threads = num_threads)
    if (ok) {
      safe_step <- mid_step
    } else {
      far_step <- mid_step
    }
  }

  if (safe_step > 0) {
    pt <- center + sign * safe_step * direction
    pmax(pmin(pt, wide$ub), wide$lb)
  } else {
    NULL
  }
}


# Wide initial bounds, scaled by crown age and tree size.
#
# Key constraint: expected tips ~ 2*exp(r*T) where r = lambda - mu.
# For observed N tips and crown age T: r_hat ~ log(N/2)/T.
# Upper bound on net rate: r_max such that 2*exp(r_max*T) ~ 50*N.
# This prevents the probe from wasting time on explosive combos.
#
# When covariates are active (DD/PD/EP), the intercept can be much higher
# than the CR feasibility limit because covariate slopes regulate growth.
# E.g. DD with lambda_0=0.7 and beta_N=-0.03 has K~23, perfectly feasible.
# We multiply lam_hi by 5x for covariate models to avoid excluding the
# true parameter region.
.wide_bounds <- function(model_bin, link_int, max_t, n_tips) {
  active <- which(model_bin == 1L)
  T <- max(max_t, 0.1)
  N <- max(n_tips, 2)

  r_hat  <- log(N / 2) / T                # observed net rate
  r_max  <- log(50 * N / 2) / T           # net rate for 50x tips
  lam_hi <- max(r_max + 0.5 * r_hat, 0.1) # allow some extinction

  # Covariate models can have much higher intercepts (slopes regulate growth)
  if (length(active) > 0L) lam_hi <- lam_hi * 5

  lam_lo <- max(0.1 * r_hat, 1e-4)
  cov_hi <- max(0.3, 3 * abs(r_hat))

  if (link_int == 0L) {
    # Linear link: rates are direct (non-negative)
    lb_lam <- c(lam_lo, rep(-cov_hi, length(active)))
    ub_lam <- c(lam_hi, rep(cov_hi, length(active)))
    lb_mu  <- c(0.0, rep(-cov_hi, length(active)))
    ub_mu  <- c(lam_hi * 0.9, rep(cov_hi, length(active)))
  } else {
    # Exponential link: log-scale
    lb_lam <- c(log(lam_lo), rep(-cov_hi, length(active)))
    ub_lam <- c(log(lam_hi), rep(cov_hi, length(active)))
    lb_mu  <- c(log(lam_lo) - 3, rep(-cov_hi, length(active)))
    ub_mu  <- c(log(lam_hi), rep(cov_hi, length(active)))
  }

  list(lb = c(lb_lam, lb_mu), ub = c(ub_lam, ub_mu))
}


# Latin Hypercube Sample: n points in [lb, ub]
.lhs_sample <- function(n, lb, ub) {
  d <- length(lb)
  mat <- matrix(0, nrow = n, ncol = d)
  for (j in seq_len(d)) {
    perm <- sample.int(n)
    mat[, j] <- lb[j] + (perm - stats::runif(n)) / n * (ub[j] - lb[j])
  }
  mat
}


# =========================================================================== #
#  Use case 2: GAM-based MLE                                                   #
# =========================================================================== #

#' Estimate log-likelihood at a grid of parameter values via importance sampling
#'
#' For each row of \code{pars_mat}, augments \code{sample_size} trees via IS
#' and computes the IS log-likelihood estimate (fhat).  Returns a data frame
#' ready for \code{\link{train_likelihood_GAM}}.
#'
#' @param tree A \code{phylo} object, branching-time vector, or anything
#'   accepted by \code{\link{estimate_rates}}.
#' @param pars_mat Numeric matrix whose rows are compact parameter vectors
#'   (same layout as \code{pars} in \code{\link{estimate_rates}}).
#' @param model Model specification: \code{"cr"}, \code{"dd"}, \code{"pd"},
#'   \code{"ep"}, or a length-3 binary integer vector.
#' @param sample_size Number of augmented trees per grid point (default 200).
#' @param link \code{"linear"} or \code{"exponential"}.
#' @param max_missing Maximum missing lineages per augmentation (default 1e4).
#' @param max_lambda Maximum speciation rate (default 500).
#' @param maxN Total augmentation attempts per grid point.  If \code{NULL},
#'   set automatically.
#' @param num_threads Number of threads for the C++ E-step.
#' @param bias_correct Logical; if \code{TRUE}, use the moment-based
#'   bias-corrected IS estimator (truncated at order \code{K}).
#'   Default \code{FALSE} uses the standard log-mean-exp estimator.
#' @param verbose Print progress (default \code{TRUE}).
#' @return A data frame with columns for each parameter plus \code{fhat}
#'   (IS log-likelihood) and \code{n_trees} (number of valid trees obtained).
#'   Rows where IS failed have \code{fhat = NA}.
#' @seealso \code{\link{train_likelihood_GAM}}, \code{\link{find_MLE}}
#' @keywords internal
estimate_likelihood_surface <- function(tree, pars_mat, model = "cr",
                                        sample_size = 200L,
                                        link = "linear",
                                        max_missing = 1e4,
                                        max_lambda = 500,
                                        maxN = NULL,
                                        num_threads = 1L,
                                        bias_correct = FALSE,
                                        verbose = TRUE,
                                        max_time = NULL) {
  model_bin <- .resolve_model(model)
  link_int  <- .resolve_link(link)
  brts      <- .extract_brts(tree)

  if (!is.matrix(pars_mat)) stop("'pars_mat' must be a matrix.")
  n_pts <- nrow(pars_mat)

  pars_df <- as.data.frame(pars_mat)
  if (is.null(colnames(pars_mat)) || any(colnames(pars_mat) == "")) {
    colnames(pars_df) <- .par_names(model_bin)
  }

  fhat    <- rep(NA_real_, n_pts)
  n_trees <- rep(0L, n_pts)

  # Worker function: evaluate IS log-likelihood at one grid point
  eval_one_point <- function(i) {
    pars8 <- .expand_pars(pars_mat[i, ], model_bin)
    maxN_i <- if (is.null(maxN)) max(2000L, 200L * as.integer(sample_size))
              else as.integer(maxN)
    raw <- tryCatch(
      augment_trees(
        brts = brts, pars = pars8,
        sample_size = as.integer(sample_size), maxN = maxN_i,
        max_missing = as.integer(max_missing),
        max_lambda = as.numeric(max_lambda),
        num_threads = if (use_r_parallel) 1L else as.integer(num_threads),
        model = as.integer(model_bin), link = as.integer(link_int)
      ),
      error = function(e) NULL
    )
    if (!is.null(raw)) {
      logf_i <- eval_logf(pars8, raw$trees,
                          model = as.integer(model_bin),
                          link = as.integer(link_int))
      list(
        fhat = .is_fhat(logf_i$logf, logf_i$logg,
                        bias_correct = bias_correct,
                        n_zero_weight = .n0(raw$rejected_zero_weights)),
        n_trees = length(logf_i$logf)
      )
    } else {
      list(fhat = NA_real_, n_trees = 0L)
    }
  }

  t0_grid <- proc.time()[3]

  # Adaptive: parallelize at R level when there's enough work per point
  # (sample_size >= 10 and enough points to fill cores).
  # Otherwise pass num_threads to C++ TBB for within-point parallelism.
  use_r_parallel <- (num_threads > 1L &&
                     .Platform$OS.type == "unix" &&
                     sample_size >= 10L &&
                     n_pts >= num_threads)

  if (use_r_parallel) {
    if (verbose) cat(sprintf(
      "  IS grid: %d points x %d trees (%d cores)...\n",
      n_pts, sample_size, num_threads
    ))
    results_list <- parallel::mclapply(
      seq_len(n_pts), eval_one_point,
      mc.cores = num_threads
    )
    for (i in seq_len(n_pts)) {
      res <- results_list[[i]]
      if (!inherits(res, "try-error") && !is.null(res)) {
        fhat[i]    <- res$fhat
        n_trees[i] <- res$n_trees
      }
    }
  } else {
    # Sequential with progress bar
    if (verbose) {
      pb <- progress::progress_bar$new(
        format = "  IS grid [:bar] :current/:total (:percent) eta: :eta",
        total = n_pts, clear = FALSE)
    }
    for (i in seq_len(n_pts)) {
      res <- eval_one_point(i)
      fhat[i]    <- res$fhat
      n_trees[i] <- res$n_trees
      if (verbose) pb$tick()
      .check_time_budget(t0_grid, i, n_pts, max_time,
                         label = "grid point")
    }
  }

  if (verbose) cat(sprintf(
    "  grid done: %.1fs, %d/%d valid\n",
    proc.time()[3] - t0_grid, sum(is.finite(fhat)), n_pts
  ))

  cbind(pars_df, fhat = fhat, n_trees = n_trees)
}


#' Fit a GAM to the IS log-likelihood surface
#'
#' Takes the output of \code{\link{estimate_likelihood_surface}} and fits a
#' smooth GAM approximation to the log-likelihood as a function of the model
#' parameters.
#'
#' @param surface Data frame returned by \code{\link{estimate_likelihood_surface}}.
#'   Rows with \code{fhat = NA} are dropped automatically.
#' @param par_names Character vector of column names in \code{surface} to use
#'   as predictors.  If \code{NULL} (default), all columns except \code{fhat}
#'   and \code{n_trees} are used.
#' @param spline_type \code{"univariate"} (default) or \code{"bivariate"}.
#'   See \code{\link{train_GAM}} for details.
#' @return A \code{gam} object (Gaussian family) representing
#'   \eqn{\hat\ell(\theta) \approx \log L(\theta \mid \text{data})}.
#' @seealso \code{\link{estimate_likelihood_surface}}, \code{\link{find_MLE}}
#' @keywords internal
train_likelihood_GAM <- function(surface,
                                 par_names = NULL,
                                 spline_type = c("univariate", "bivariate")) {
  spline_type <- match.arg(spline_type)
  surface <- surface[is.finite(surface$fhat), , drop = FALSE]
  if (nrow(surface) < 10L)
    stop("Too few valid likelihood estimates (", nrow(surface),
         "); need at least 10. Try a larger grid or more IS samples.")

  if (is.null(par_names)) {
    par_names <- setdiff(colnames(surface), c("fhat", "n_trees"))
  }

  varying <- par_names[vapply(par_names,
    function(p) length(unique(surface[[p]])) > 1L, logical(1L))]
  if (length(varying) == 0L) stop("All parameter columns are constant.")

  # Cap basis dimension at the number of unique values per predictor
  n_unique <- vapply(varying,
    function(p) length(unique(surface[[p]])), integer(1L))

  if (spline_type == "univariate" || length(varying) < 2L) {
    smooth_terms <- vapply(varying, function(p) {
      k <- min(n_unique[[p]], 10L)
      if (k < 3L) k <- 3L
      sprintf("s(%s, k=%d)", p, k)
    }, character(1L))
    smooth_terms <- paste(smooth_terms, collapse = " + ")
  } else {
    pairs <- utils::combn(varying, 2L, simplify = FALSE)
    smooth_terms <- paste(
      vapply(pairs, function(p) {
        k1 <- min(n_unique[[p[1L]]], 5L)
        k2 <- min(n_unique[[p[2L]]], 5L)
        sprintf("te(%s, %s, k=c(%d,%d))", p[1L], p[2L], max(3L, k1), max(3L, k2))
      }, character(1L)),
      collapse = " + ")
  }

  form <- stats::as.formula(paste("fhat ~", smooth_terms))
  cat("Training", spline_type, "likelihood GAM on", nrow(surface),
      "grid points...\n")
  mgcv::gam(form, data = surface, family = stats::gaussian())
}


#' Find MLE by optimizing a fitted likelihood GAM
#'
#' Uses \code{\link[stats]{optim}} (L-BFGS-B) to maximize the GAM-predicted
#' log-likelihood surface.
#'
#' @param gam_fit A \code{gam} object returned by
#'   \code{\link{train_likelihood_GAM}}.
#' @param lower_bound Numeric vector of lower bounds for each parameter.
#' @param upper_bound Numeric vector of upper bounds for each parameter.
#' @param start Optional numeric vector of starting values.  If \code{NULL},
#'   the grid point with the highest predicted value is used.
#' @param par_names Character vector of parameter names matching the GAM
#'   predictors.  If \code{NULL}, inferred from the GAM formula.
#' @return A list with components:
#'   \describe{
#'     \item{pars}{Named numeric vector of MLE estimates.}
#'     \item{loglik}{Predicted log-likelihood at the MLE.}
#'     \item{convergence}{Convergence code from \code{optim} (0 = success).}
#'     \item{optim}{Full \code{optim} output.}
#'   }
#' @seealso \code{\link{train_likelihood_GAM}},
#'   \code{\link{estimate_likelihood_surface}}
#' @keywords internal
find_MLE <- function(gam_fit, lower_bound, upper_bound,
                     start = NULL, par_names = NULL) {
  if (is.null(par_names)) {
    # Use var.summary (robust for GAMs with s()/te() terms)
    # all.vars() can pick up spurious names like 'k' from s(x, k=10)
    par_names <- names(gam_fit$var.summary)
  }
  if (length(lower_bound) != length(par_names) ||
      length(upper_bound) != length(par_names))
    stop("'lower_bound' and 'upper_bound' must match the number of GAM predictors (",
         length(par_names), ").")

  neg_loglik <- function(pars) {
    nd <- as.data.frame(as.list(stats::setNames(pars, par_names)))
    -stats::predict(gam_fit, newdata = nd, type = "response")
  }

  if (is.null(start)) {
    train_data <- gam_fit$model
    best_idx   <- which.max(stats::fitted(gam_fit))
    start      <- as.numeric(train_data[best_idx, par_names])
  }

  opt <- stats::optim(
    par = start, fn = neg_loglik,
    method = "L-BFGS-B",
    lower = lower_bound, upper = upper_bound
  )

  list(
    pars        = stats::setNames(opt$par, par_names),
    loglik      = -opt$value,
    convergence = opt$convergence,
    optim       = opt
  )
}


# .par_names() is defined in inference.R -- shared by both modules

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
#' Explores parameter space by drawing random parameter vectors from a wide
#' initial region, simulating trees at each, and identifying the subregion
#' that produces "sensible" trees (survived, reasonable size).  Optionally
#' trains a survival GAM over the feasible region.
#'
#' @param tree A \code{phylo} object or anything accepted by
#'   \code{\link{estimate_rates}}.  Used to determine crown age and
#'   target tree size.
#' @param model Model specification: \code{"cr"}, \code{"dd"}, etc.
#' @param link \code{"linear"} or \code{"exponential"}.
#' @param n_probe Number of parameter vectors to simulate (default 500).
#' @param margin Fractional expansion of the detected feasible range
#'   (default 0.2 = 20\% on each side).
#' @param tip_range Multiplier range for acceptable tip counts relative
#'   to the observed tree.  Default \code{c(0.1, 10)}: accept trees
#'   with 10\% to 10x the observed number of tips.
#' @param train_surv_gam If \code{TRUE} (default), train a survival GAM
#'   on the feasible region using a second round of simulations.
#' @param verbose Print progress (default \code{TRUE}).
#' @return A list with components:
#'   \describe{
#'     \item{lower_bound}{Numeric vector of lower bounds (compact layout).}
#'     \item{upper_bound}{Numeric vector of upper bounds (compact layout).}
#'     \item{survival_gam}{A \code{gam} object (if \code{train_surv_gam = TRUE}),
#'       or \code{NULL}.}
#'     \item{n_feasible}{Number of feasible parameter vectors found.}
#'     \item{n_probed}{Total parameter vectors simulated.}
#'     \item{feasible_pars}{Matrix of feasible parameter vectors.}
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
                        n_probe = 500L, margin = 0.2,
                        tip_range = c(0.1, 10),
                        train_surv_gam = TRUE,
                        verbose = TRUE) {
  model_bin <- .resolve_model(model)
  link_int  <- .resolve_link(link)
  n_pars    <- 2L + 2L * sum(model_bin)
  pnames    <- .par_names(model_bin)

  # Crown age and target tips from observed tree
  if (inherits(tree, "phylo")) {
    max_t   <- max(ape::branching.times(tree))
    n_tips  <- ape::Ntip(tree)
  } else if (is.numeric(tree)) {
    max_t   <- max(tree)
    n_tips  <- length(tree) + 1L
  } else if (is.list(tree) && !is.null(tree$tes)) {
    max_t   <- max(ape::branching.times(tree$tes))
    n_tips  <- ape::Ntip(tree$tes)
  } else {
    stop("'tree' must be a phylo, branching times, or simulate_tree result.")
  }

  tip_lo <- max(2, floor(n_tips * tip_range[1]))
  tip_hi <- ceiling(n_tips * tip_range[2])

  # ── Wide initial bounds (model + link + crown-age dependent) ─
  wide <- .wide_bounds(model_bin, link_int, max_t, n_tips)

  # ── Phase 1: probe with LHS ────────────────────────────────
  if (verbose) cat(sprintf(
    "auto_bounds: probing %d points (%s, link=%s, %d pars)...\n",
    n_probe, paste(model_bin, collapse = ","),
    if (link_int == 0L) "linear" else "exponential", n_pars
  ))

  pars_mat <- .lhs_sample(n_probe, wide$lb, wide$ub)
  colnames(pars_mat) <- pnames

  max_lin <- as.integer(max(20 * n_tips, 500))
  sims <- simulate_tree(
    pars = pars_mat, max_t = max_t,
    model = model, link = link,
    max_tries = 0,
    max_lin = max_lin
  )

  # ── Classify outcomes ──────────────────────────────────────
  status <- sapply(sims$simulations, function(s) s$status)
  ntips  <- sapply(sims$simulations, function(s) {
    if (s$status == "done" && !is.null(s$tes)) {
      length(s$tes$tip.label)
    } else {
      0L
    }
  })

  feasible <- (status == "done") & (ntips >= tip_lo) & (ntips <= tip_hi)
  n_feasible <- sum(feasible)

  if (verbose) cat(sprintf(
    "  survived: %d/%d, sensible tips [%.0f, %.0f]: %d\n",
    sum(status == "done"), n_probe, tip_lo, tip_hi, n_feasible
  ))

  if (n_feasible < 5L) {
    warning(
      "auto_bounds: only ", n_feasible,
      " feasible points found. Returning wide bounds. ",
      "Try increasing n_probe or relaxing tip_range."
    )
    return(list(
      lower_bound  = wide$lb,
      upper_bound  = wide$ub,
      survival_gam = NULL,
      n_feasible   = n_feasible,
      n_probed     = n_probe,
      feasible_pars = pars_mat[feasible, , drop = FALSE]
    ))
  }

  # ── Narrow bounds to feasible region + margin ──────────────
  feas_mat <- pars_mat[feasible, , drop = FALSE]
  col_min  <- apply(feas_mat, 2, min)
  col_max  <- apply(feas_mat, 2, max)
  span     <- col_max - col_min
  span     <- pmax(span, 1e-6)  # avoid zero-width

  lb <- pmax(col_min - margin * span, wide$lb)
  ub <- pmin(col_max + margin * span, wide$ub)

  if (verbose) {
    cat("  bounds:\n")
    for (j in seq_len(n_pars)) {
      cat(sprintf("    %s: [%.4f, %.4f]\n", pnames[j], lb[j], ub[j]))
    }
  }

  # ── Phase 2: survival GAM on narrowed region ───────────────
  surv_gam <- NULL
  if (train_surv_gam) {
    if (verbose) cat("  training survival GAM on narrowed region...\n")
    n_gam   <- max(500L, n_probe)
    gam_mat <- .lhs_sample(n_gam, lb, ub)
    colnames(gam_mat) <- pnames

    gam_sims <- simulate_tree(
      pars = gam_mat, max_t = max_t,
      model = model, link = link,
      max_tries = 0,
      max_lin = max_lin
    )

    surv_gam <- train_GAM(gam_sims$simulations, gam_mat, model = model)
  }

  list(
    lower_bound   = stats::setNames(lb, pnames),
    upper_bound   = stats::setNames(ub, pnames),
    survival_gam  = surv_gam,
    n_feasible    = n_feasible,
    n_probed      = n_probe,
    feasible_pars = feas_mat
  )
}


# Wide initial bounds, scaled by crown age and tree size.
#
# Key constraint: expected tips ~ 2*exp(r*T) where r = lambda - mu.
# For observed N tips and crown age T: r_hat ~ log(N/2)/T.
# Upper bound on net rate: r_max such that 2*exp(r_max*T) ~ 50*N.
# This prevents the probe from wasting time on explosive combos.
.wide_bounds <- function(model_bin, link_int, max_t, n_tips) {
  active <- which(model_bin == 1L)
  T <- max(max_t, 0.1)
  N <- max(n_tips, 2)

  r_hat  <- log(N / 2) / T                # observed net rate
  r_max  <- log(50 * N / 2) / T           # net rate for 50x tips
  lam_hi <- max(r_max + 0.5 * r_hat, 0.1) # allow some extinction
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

  if (verbose) {
    pb <- progress::progress_bar$new(
      format = "  IS grid [:bar] :current/:total (:percent) eta: :eta",
      total = n_pts, clear = FALSE)
  }

  t0_grid <- proc.time()[3]
  for (i in seq_len(n_pts)) {
    pars8 <- .expand_pars(pars_mat[i, ], model_bin)
    maxN_i <- if (is.null(maxN)) max(2000L, 200L * as.integer(sample_size))
              else as.integer(maxN)
    raw <- tryCatch(
      augment_trees(
        brts = brts, pars = pars8,
        sample_size = as.integer(sample_size), maxN = maxN_i,
        max_missing = as.integer(max_missing),
        max_lambda = as.numeric(max_lambda),
        num_threads = as.integer(num_threads),
        model = as.integer(model_bin), link = as.integer(link_int)
      ),
      error = function(e) NULL
    )
    if (!is.null(raw)) {
      logf_i <- eval_logf(pars8, raw$trees,
                          model = as.integer(model_bin),
                          link = as.integer(link_int))
      fhat[i]    <- .is_fhat(logf_i$logf, logf_i$logg,
                             bias_correct = bias_correct,
                             n_zero_weight = .n0(raw$rejected_zero_weights))
      n_trees[i] <- length(logf_i$logf)
    }
    if (verbose) pb$tick()
    .check_time_budget(t0_grid, i, n_pts, max_time, label = "grid point")
  }

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

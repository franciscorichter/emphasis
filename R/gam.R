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

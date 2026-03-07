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


# .par_names() is defined in inference.R -- shared by both modules

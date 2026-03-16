#' Full estimation pipeline: auto_bounds -> GAM -> CEM -> MCEM
#'
#' Runs a multi-stage estimation workflow with automatic bounds detection,
#' exploratory fits (GAM and CEM), and optional MCEM refinement.  Each stage
#' is timed and logged.  If MCEM fails, the best exploratory result is
#' returned.
#'
#' @section Stages:
#' \enumerate{
#'   \item \strong{Bounds}: \code{\link{auto_bounds}} detects parameter bounds
#'     and trains a survival GAM for conditioning.
#'   \item \strong{GAM}: Fits a smooth likelihood surface over a grid and
#'     optimises it.  Fast, one-shot.
#'   \item \strong{CEM}: Cross-entropy global search.  Stochastic, no init
#'     needed.  Explores the space differently than GAM.
#'   \item \strong{MCEM}: Monte Carlo EM starting from the best init found
#'     in stages 2--3.  Local refinement for precise estimates.
#' }
#'
#' @param tree A \code{phylo} object, branching-time vector, or
#'   \code{\link{simulate_tree}} result.
#' @param model Model specification: \code{"cr"}, \code{"dd"}, \code{"pd"},
#'   \code{"ep"}, a formula, or a binary vector.
#' @param link \code{"linear"} (default) or \code{"exponential"}.
#' @param stages Character vector of stages to run.  Default
#'   \code{c("bounds", "gam", "cem", "mcem")}.  Remove stages to skip them,
#'   e.g. \code{c("bounds", "gam")} for a quick exploratory fit.
#' @param control Named list with optional sub-lists \code{bounds}, \code{gam},
#'   \code{cem}, \code{mcem} to override defaults for each stage.  Common
#'   parameters (\code{num_threads}, \code{max_missing}, \code{max_time}) can
#'   be set at the top level and are inherited by all stages.
#' @param verbose Logical.  If \code{TRUE} (default), print progress and a
#'   summary log to the console.
#' @return A list of class \code{"emphasis_pipeline"} with:
#'   \describe{
#'     \item{\code{pars}}{Named parameter vector from the best stage.}
#'     \item{\code{loglik}}{Log-likelihood from the best stage.}
#'     \item{\code{AIC}}{AIC from the best stage.}
#'     \item{\code{best_stage}}{Which stage produced the final result.}
#'     \item{\code{log}}{Data frame with one row per stage: stage name,
#'       pars, loglik, AIC, elapsed seconds, status.}
#'     \item{\code{mcem_trace}}{If MCEM ran, a data frame of per-iteration
#'       fhat and parameter values.  \code{NULL} otherwise.}
#'     \item{\code{bounds}}{The \code{auto_bounds} result (bounds + survival GAM).}
#'     \item{\code{fits}}{Named list of individual \code{emphasis_fit} objects.}
#'   }
#' @examples
#' \dontrun{
#' library(ape)
#' data(bird.orders)
#' result <- emphasis_pipeline(bird.orders, model = "dd", link = "linear")
#' result           # prints summary log
#' result$pars      # best parameters
#' result$log       # full stage log
#' }
#' @export
emphasis_pipeline <- function(tree,
                              model   = "cr",
                              link    = "linear",
                              stages  = c("bounds", "gam", "cem", "mcem"),
                              control = list(),
                              verbose = TRUE) {

  model_bin <- .resolve_model(model)
  model_str <- .model_label(model_bin)
  link_int  <- .resolve_link(link)
  link_str  <- if (link_int == 0L) "linear" else "exponential"
  brts      <- .extract_brts(tree)
  n_brts    <- length(brts)

  # Shared defaults
  num_threads <- control$num_threads %||% max(1L, parallel::detectCores() - 1L)
  max_missing <- control$max_missing %||% 1e4
  max_time    <- control$max_time %||% 3600

  if (verbose) cat(sprintf(
    "[pipeline] model=%s  link=%s  n_brts=%d  threads=%d\n",
    model_str, link_str, n_brts, num_threads
  ))

  # Result accumulators
  run_log    <- data.frame()
  fits       <- list()
  ab         <- NULL
  mcem_trace <- NULL

  # Helper: append to log
  log_stage <- function(stage, fit, elapsed, status) {
    pars_str <- if (!is.null(fit) && !is.null(fit$pars))
      paste(names(fit$pars), "=", round(fit$pars, 4), collapse = ", ") else ""
    ll  <- if (!is.null(fit)) fit$loglik else NA_real_
    aic <- if (!is.null(fit)) fit$AIC else NA_real_

    row <- data.frame(
      stage   = stage,
      loglik  = ll,
      AIC     = aic,
      elapsed = elapsed,
      status  = status,
      pars    = pars_str,
      stringsAsFactors = FALSE
    )
    run_log <<- rbind(run_log, row)

    if (verbose) {
      cat(sprintf("  [%s] %s  loglik=%.2f  AIC=%.2f  (%.1fs)\n",
                  stage, status,
                  if (is.finite(ll)) ll else NA,
                  if (is.finite(aic)) aic else NA,
                  elapsed))
      if (nchar(pars_str) > 0) cat(sprintf("    pars: %s\n", pars_str))
    }
  }

  # ── Stage 1: Bounds ──────────────────────────────────────────
  if ("bounds" %in% stages) {
    if (verbose) cat("\n[Stage 1] auto_bounds...\n")
    t0 <- proc.time()[3]
    ab <- tryCatch(
      auto_bounds(tree, model = model, link = link,
                  num_threads = num_threads,
                  margin  = control$bounds$margin %||% 0.5,
                  n_test  = control$bounds$n_test %||% 5L,
                  verbose = verbose),
      error = function(e) {
        if (verbose) cat(sprintf("  bounds FAILED: %s\n", e$message))
        NULL
      }
    )
    elapsed_bounds <- proc.time()[3] - t0

    if (is.null(ab)) {
      log_stage("bounds", NULL, elapsed_bounds, "failed")
      # Cannot continue without bounds
      result <- list(
        pars = NULL, loglik = NA_real_, AIC = NA_real_,
        best_stage = "none", log = run_log, mcem_trace = NULL,
        bounds = NULL, fits = fits
      )
      class(result) <- "emphasis_pipeline"
      return(result)
    }

    log_stage("bounds", NULL, elapsed_bounds, "ok")
    if (verbose) {
      cat(sprintf("    lb: %s\n", paste(round(ab$lower_bound, 4), collapse = ", ")))
      cat(sprintf("    ub: %s\n", paste(round(ab$upper_bound, 4), collapse = ", ")))
    }
  } else {
    # Bounds must be provided in control
    if (is.null(control$lower_bound) || is.null(control$upper_bound))
      stop("When 'bounds' stage is skipped, control must contain 'lower_bound' and 'upper_bound'.")
    ab <- list(
      lower_bound  = control$lower_bound,
      upper_bound  = control$upper_bound,
      survival_gam = control$survival_gam,
      center       = (control$lower_bound + control$upper_bound) / 2
    )
  }

  cond <- ab$survival_gam

  # ── Stage 2: GAM ────────────────────────────────────────────
  if ("gam" %in% stages) {
    if (verbose) cat("\n[Stage 2] GAM surface fit...\n")

    gam_defaults <- list(n_grid = 150, sample_size = 200)
    gam_ctrl <- utils::modifyList(gam_defaults, control$gam %||% list())
    gam_ctrl$lower_bound <- ab$lower_bound
    gam_ctrl$upper_bound <- ab$upper_bound
    gam_ctrl$num_threads <- num_threads
    gam_ctrl$max_missing <- max_missing
    gam_ctrl$max_time    <- max_time
    gam_ctrl$verbose     <- verbose

    t0 <- proc.time()[3]
    gam_fit <- tryCatch(
      estimate_rates(tree, method = "gam", model = model, link = link,
                     cond = cond, control = gam_ctrl),
      error = function(e) {
        if (verbose) cat(sprintf("  GAM FAILED: %s\n", e$message))
        NULL
      }
    )
    elapsed_gam <- proc.time()[3] - t0

    gam_ok <- !is.null(gam_fit) && is.finite(gam_fit$loglik)
    fits$gam <- gam_fit
    log_stage("gam", if (gam_ok) gam_fit else NULL,
              elapsed_gam, if (gam_ok) "ok" else "failed")
  }

  # ── Stage 3: CEM ────────────────────────────────────────────
  if ("cem" %in% stages) {
    if (verbose) cat("\n[Stage 3] CEM global search...\n")

    cem_defaults <- list(max_iter = 20, num_particles = 50, num_trees = 5)
    cem_ctrl <- utils::modifyList(cem_defaults, control$cem %||% list())
    cem_ctrl$lower_bound <- ab$lower_bound
    cem_ctrl$upper_bound <- ab$upper_bound
    cem_ctrl$num_threads <- num_threads
    cem_ctrl$max_missing <- max_missing
    cem_ctrl$max_time    <- max_time
    cem_ctrl$verbose     <- verbose

    t0 <- proc.time()[3]
    cem_fit <- tryCatch(
      estimate_rates(tree, method = "cem", model = model, link = link,
                     cond = cond, control = cem_ctrl),
      error = function(e) {
        if (verbose) cat(sprintf("  CEM FAILED: %s\n", e$message))
        NULL
      }
    )
    elapsed_cem <- proc.time()[3] - t0

    cem_ok <- !is.null(cem_fit) && is.finite(cem_fit$loglik)
    fits$cem <- cem_fit
    log_stage("cem", if (cem_ok) cem_fit else NULL,
              elapsed_cem, if (cem_ok) "ok" else "failed")
  }

  # ── Pick best init for MCEM ─────────────────────────────────
  best_init     <- NULL
  best_init_ll  <- -Inf
  best_init_src <- "midpoint"

  for (stage_name in c("gam", "cem")) {
    f <- fits[[stage_name]]
    if (!is.null(f) && is.finite(f$loglik) && f$loglik > best_init_ll) {
      best_init    <- f$pars
      best_init_ll <- f$loglik
      best_init_src <- stage_name
    }
  }

  if (verbose && !is.null(best_init)) {
    cat(sprintf("\n  Best init from %s (loglik=%.2f)\n", best_init_src, best_init_ll))
  }

  # ── Stage 4: MCEM ───────────────────────────────────────────
  if ("mcem" %in% stages) {
    if (verbose) cat("\n[Stage 4] MCEM refinement...\n")

    mcem_defaults <- list(sample_size = 200, maxN = 5000, max_iter = 200,
                          tol = 1e-3, patience = 3)
    mcem_ctrl <- utils::modifyList(mcem_defaults, control$mcem %||% list())
    mcem_ctrl$lower_bound <- ab$lower_bound
    mcem_ctrl$upper_bound <- ab$upper_bound
    mcem_ctrl$num_threads <- num_threads
    mcem_ctrl$max_missing <- max_missing
    mcem_ctrl$max_time    <- max_time
    mcem_ctrl$verbose     <- verbose

    t0 <- proc.time()[3]
    mcem_fit <- tryCatch(
      estimate_rates(tree, method = "mcem", model = model, link = link,
                     init_pars = best_init, cond = cond, control = mcem_ctrl),
      error = function(e) {
        if (verbose) cat(sprintf("  MCEM FAILED: %s\n", e$message))
        NULL
      }
    )
    elapsed_mcem <- proc.time()[3] - t0

    mcem_ok <- !is.null(mcem_fit) && is.finite(mcem_fit$loglik)
    fits$mcem <- mcem_fit
    log_stage("mcem", if (mcem_ok) mcem_fit else NULL,
              elapsed_mcem, if (mcem_ok) "ok" else "failed")

    # Extract MCEM iteration trace
    if (mcem_ok && !is.null(mcem_fit$details) && !is.null(mcem_fit$details$mcem)) {
      mcem_trace <- mcem_fit$details$mcem
    }

    # Log per-iteration progress if MCEM ran
    if (verbose && !is.null(mcem_trace) && nrow(mcem_trace) > 0) {
      cat(sprintf("    MCEM: %d iterations, final fhat=%.2f\n",
                  nrow(mcem_trace),
                  utils::tail(mcem_trace$fhat[is.finite(mcem_trace$fhat)], 1)))
    }
  }

  # ── Select best result ──────────────────────────────────────
  best_fit   <- NULL
  best_stage <- "none"

  for (stage_name in c("mcem", "cem", "gam")) {
    f <- fits[[stage_name]]
    if (!is.null(f) && is.finite(f$loglik)) {
      if (is.null(best_fit) || f$loglik > best_fit$loglik) {
        best_fit   <- f
        best_stage <- stage_name
      }
    }
  }

  # ── Summary ─────────────────────────────────────────────────
  total_time <- sum(run_log$elapsed)

  if (verbose) {
    cat(sprintf("\n[pipeline] DONE in %.1fs  best=%s  loglik=%.2f  AIC=%.2f\n",
                total_time, best_stage,
                if (!is.null(best_fit)) best_fit$loglik else NA,
                if (!is.null(best_fit)) best_fit$AIC else NA))
  }

  result <- list(
    pars       = if (!is.null(best_fit)) best_fit$pars else NULL,
    loglik     = if (!is.null(best_fit)) best_fit$loglik else NA_real_,
    AIC        = if (!is.null(best_fit)) best_fit$AIC else NA_real_,
    n_pars     = if (!is.null(best_fit)) best_fit$n_pars else NA_integer_,
    best_stage = best_stage,
    method     = best_stage,
    model      = model_bin,
    cond       = !is.null(cond),
    log        = run_log,
    mcem_trace = mcem_trace,
    bounds     = ab,
    fits       = fits
  )
  class(result) <- c("emphasis_pipeline", "emphasis_fit")
  result
}


#' Print method for emphasis_pipeline objects
#'
#' @param x An \code{emphasis_pipeline} result.
#' @param ... Additional arguments (ignored).
#' @export
print.emphasis_pipeline <- function(x, ...) {
  model_str <- .model_label(x$model)
  cond_str <- if (isTRUE(x$cond)) ", conditioned" else ""
  cat("emphasis pipeline (model = ", model_str, cond_str, ")\n", sep = "")
  cat(strrep("=", 50), "\n\n")

  # Stage log
  cat("Stage log:\n")
  for (i in seq_len(nrow(x$log))) {
    r <- x$log[i, ]
    marker <- if (r$stage == x$best_stage) " *" else "  "
    cat(sprintf("  %s%-6s  %6s  loglik=%10s  AIC=%10s  %5.1fs\n",
                marker, r$stage, r$status,
                if (is.finite(r$loglik)) sprintf("%.2f", r$loglik) else "---",
                if (is.finite(r$AIC)) sprintf("%.2f", r$AIC) else "---",
                r$elapsed))
  }

  cat(sprintf("\nBest stage: %s\n", x$best_stage))
  if (!is.null(x$pars)) {
    cat("Parameters:\n")
    print(round(x$pars, 6))
  }
  cat(sprintf("\nLog-likelihood: %.4f\n", x$loglik))
  cat(sprintf("AIC:            %.4f\n", x$AIC))
  cat(sprintf("Total time:     %.1fs\n", sum(x$log$elapsed)))
  invisible(x)
}

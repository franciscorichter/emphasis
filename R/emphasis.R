.mcem_dynamic_fresh <- function(brts,
                      pars,
                      sample_size,
                      maxN,
                      max_missing,
                      lower_bound,
                      upper_bound,
                      max_iter,
                      xtol,
                      tol,
                      patience,
                      num_threads,
                      verbose = FALSE,
                      conditional = NULL,
                      model = c(0L, 0L, 0L),
                      link = 0L) {
  if (inherits(brts, "phylo")) {
    brts <- sort(ape::branching.times(brts), decreasing = TRUE)
  }
  if (!is.numeric(brts)) stop("`brts` must be numeric or a `phylo` object.")
  if (!is.numeric(pars)) stop("`pars` must be numeric.")
  if (!is.null(conditional) && !is.function(conditional))
    stop("`conditional` must be a function or NULL.")

  # Scale parameter changes by the search range for convergence check
  range_vec <- upper_bound - lower_bound
  range_vec[range_vec == 0] <- 1  # fixed params (lb == ub): delta always 0

  streak      <- 0L
  fail_streak <- 0L
  prev_pars   <- pars
  mcem        <- NULL
  last_results <- NULL
  stop_reason <- "max_iter"

  for (i in seq_len(max_iter)) {
    results <- tryCatch(
      em_cpp(brts = brts,
             init_pars = pars,
             sample_size = sample_size,
             maxN = maxN,
             max_missing = max_missing,
             max_lambda = 1e6,
             lower_bound = lower_bound,
             upper_bound = upper_bound,
             xtol_rel = xtol,
             num_threads = num_threads,
             copy_trees = FALSE,
             model = as.integer(model),
             link = as.integer(link),
             rconditional = conditional),
      error = function(e) NULL
    )

    # Handle E-step failure (all trees rejected)
    if (is.null(results)) {
      fail_streak <- fail_streak + 1L
      if (verbose) message(sprintf("Iteration %d: E-step failed (%d consecutive)", i, fail_streak))
      if (fail_streak >= 5L) {
        stop_reason <- "e_step_failure"
        warning("MCEM: 5 consecutive E-step failures; stopping. ",
                "Try widening bounds or increasing max_missing.")
        break
      }
      next
    }
    fail_streak <- 0L

    last_results <- results
    pars <- results$estimates
    delta_max <- max(abs(pars - prev_pars) / range_vec)
    prev_pars <- pars

    # Record iteration
    par_df <- as.data.frame(as.list(stats::setNames(pars, paste0("par", seq_along(pars)))))
    step <- cbind(par_df, data.frame(
      fhat      = results$fhat,
      delta_max = delta_max,
      rejected  = results$rejected,
      num_trees = sample_size,
      time      = results$time
    ))
    mcem <- rbind(mcem, step)

    if (verbose) {
      rej_str <- if (results$rejected > 0L) sprintf("  rej=%d", results$rejected) else ""
      message(sprintf("Iteration %d: fhat=%.4f  delta=%.2e  streak=%d/%d%s",
                      i, results$fhat, delta_max, streak, patience, rej_str))
    }

    # Check convergence: parameter stability
    if (delta_max < tol) {
      streak <- streak + 1L
      if (streak >= patience) {
        stop_reason <- "converged"
        break
      }
    } else {
      streak <- 0L
    }
  }

  # Bootstrap variance from the last iteration's IS weights
  loglik_var <- NA_real_
  if (!is.null(last_results) &&
      length(last_results$logf) >= 2L &&
      all(is.finite(last_results$logf))) {
    loglik_var <- .bootstrap_fhat_var(last_results$logf, last_results$logg,
                                      K = 2L, B = 200L)
  }

  # Store final IS components for diagnostics
  final_IS <- NULL
  if (!is.null(last_results) && length(last_results$logf) > 0L) {
    lw <- last_results$logf - last_results$logg
    final_IS <- list(
      logf  = last_results$logf,
      logg  = last_results$logg,
      lw    = lw,
      fhat  = .is_fhat(last_results$logf, last_results$logg),
      ESS   = .ess_from_lw(lw)
    )
  }

  list(
    mcem        = mcem,
    pars        = pars,
    iterations  = nrow(mcem),
    stop_reason = stop_reason,
    loglik_var  = loglik_var,
    final_IS    = final_IS
  )
}

get_required_sampling_size <- function(M, tol = 0.05) {
  if (!nrow(M)) {
    stop("Input `M` must contain at least one row.")
  }
  hlp <- MASS::rlm(M$fhat ~ I(1 / M$sample_size), weights = M$sample_size)
  ab <- stats::coef(hlp)
  f_r <- ab[1] - tol
  n_r <- ceiling(ab[2] / (f_r - ab[1]))
  n_r
}

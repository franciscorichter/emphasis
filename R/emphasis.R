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
                      link = 0L,
                      max_time = NULL) {
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
  t0_mcem     <- proc.time()[3]

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
        # Try to diagnose why by extracting the last error
        .mcem_warn_estep(brts, pars, lower_bound, upper_bound, model, link)
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

    # Time budget check
    if (!is.null(max_time)) {
      elapsed <- proc.time()[3] - t0_mcem
      if (elapsed > max_time) {
        stop_reason <- "time_budget"
        if (verbose) message(sprintf("Time budget reached (%.0fs > %ds)", elapsed, as.integer(max_time)))
        break
      }
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
    n_rej <- .n0(last_results$rejected) +
              .n0(last_results$rejected_overruns) +
              .n0(last_results$rejected_lambda)
    final_IS <- list(
      logf  = last_results$logf,
      logg  = last_results$logg,
      lw    = lw,
      fhat  = .is_fhat(last_results$logf, last_results$logg,
                        n_zero_weight = .n0(last_results$rejected_zero_weights)),
      ESS   = .ess_from_lw(lw),
      n_rejected = n_rej,
      rejected_zero_weights = .n0(last_results$rejected_zero_weights)
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

# Diagnose why E-step failed and issue an informative warning.
.mcem_warn_estep <- function(brts, pars, lower_bound, upper_bound, model, link) {
  # Run a tiny diagnostic E-step to get rejection breakdown
  err_msg <- tryCatch({
    em_cpp(brts         = brts,
           init_pars    = pars,
           sample_size  = 1L,
           maxN         = 200L,
           max_missing  = 1e4,
           max_lambda   = 1e6,
           lower_bound  = lower_bound,
           upper_bound  = upper_bound,
           xtol_rel     = 1e-3,
           num_threads  = 1L,
           copy_trees   = FALSE,
           model        = as.integer(model),
           link         = as.integer(link))
    NULL
  }, error = function(e) conditionMessage(e))

  zero_w <- if (!is.null(err_msg)) {
    m <- regmatches(err_msg, regexpr("[0-9]+ zero weights", err_msg))
    if (length(m)) as.integer(sub(" zero weights", "", m)) else NA_integer_
  } else NA_integer_

  if (!is.na(zero_w) && zero_w > 50L) {
    warning(
      "MCEM: E-step failed - nearly all augmented trees have zero IS weight ",
      "(rejected_zero_weights=", zero_w, "/200).\n",
      "  Likely cause: speciation rate is zero (lambda=0) under current parameters,\n",
      "  e.g. PD or DD model with large tree and strongly negative covariate slope.\n",
      "  Suggestions:\n",
      "    1. Use link=\"exponential\" (ensures lambda > 0 everywhere).\n",
      "    2. Restrict bounds so the covariate slope cannot drive lambda to zero.\n",
      "    3. Use a simpler model (e.g. CR or DD) for this tree size.",
      call. = FALSE
    )
  } else {
    warning(
      "MCEM: 5 consecutive E-step failures; stopping. ",
      "Try widening bounds or increasing max_missing.",
      call. = FALSE
    )
  }
}

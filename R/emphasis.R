.mcem_dynamic_fresh <- function(brts,
                      pars,
                      sample_size,
                      maxN,
                      max_missing,
                      lower_bound,
                      upper_bound,
                      max_iter,
                      xtol,
                      tol = 1e-2,
                      patience = 3L,
                      num_threads,
                      verbose = FALSE,
                      conditional = NULL,
                      model = c(0L, 0L, 0L),
                      link = 0L,
                      max_time = NULL,
                      rho = 1.0) {
  if (inherits(brts, "phylo")) {
    brts <- sort(ape::branching.times(brts), decreasing = TRUE)
  }
  if (!is.numeric(brts)) stop("`brts` must be numeric or a `phylo` object.")
  if (!is.numeric(pars)) stop("`pars` must be numeric.")
  if (!is.null(conditional) && !is.function(conditional))
    stop("`conditional` must be a function or NULL.")

  # Convergence metric: max_j |theta_k,j - theta_{k-1,j}| / max(|theta_{k-1,j}|, rel_floor).
  # The bound box does not enter; the floor keeps parameters at or near zero
  # from giving an unbounded ratio.
  rel_floor <- 1e-2
  rel_delta <- function(new, old) max(abs(new - old) / pmax(abs(old), rel_floor))

  streak       <- 0L      # consecutive iterations with rel_delta < tol
  fail_streak  <- 0L      # consecutive E-step failures
  n_failed     <- 0L      # E-step failures over the whole run
  n_success    <- 0L      # completed EM iterations
  prev_pars    <- pars    # last successful M-step estimate (init until the first success)
  cur_pars     <- pars    # point the next E-step samples at and the M-step starts from
  par_hist     <- list()  # successful iterates, for the windowed drift
  mcem         <- NULL
  last_results <- NULL
  stop_reason  <- "max_iter"
  t0_mcem      <- proc.time()[3]
  center_pars  <- (lower_bound + upper_bound) / 2
  maxN_cap     <- 50000L

  run_em <- function(at, maxN_now) {
    tryCatch(
      em_cpp(brts = brts,
             init_pars = at,
             sample_size = sample_size,
             maxN = maxN_now,
             max_missing = max_missing,
             max_lambda = 1e6,
             lower_bound = lower_bound,
             upper_bound = upper_bound,
             xtol_rel = xtol,
             num_threads = num_threads,
             copy_trees = FALSE,
             model = as.integer(model),
             link = as.integer(link),
             rho = as.numeric(rho),
             rconditional = conditional),
      error = function(e) NULL
    )
  }

  # One trace row: the iterate, the E-step value, the step metrics and the
  # rejection counters of that E-step. `n_rejected` is the sum final_IS reports.
  trace_row <- function(p, results, delta_max, abs_step, drift, maxN_used, final_estep) {
    par_df <- as.data.frame(as.list(stats::setNames(p, paste0("par", seq_along(p)))))
    n_rej <- .n0(results$rejected) +
             .n0(results$rejected_overruns) +
             .n0(results$rejected_lambda)
    cbind(par_df, data.frame(
      fhat                  = results$fhat,
      delta_max             = delta_max,
      abs_step              = abs_step,
      drift                 = drift,
      rejected              = .n0(results$rejected),
      rejected_overruns     = .n0(results$rejected_overruns),
      rejected_lambda       = .n0(results$rejected_lambda),
      rejected_zero_weights = .n0(results$rejected_zero_weights),
      n_rejected            = n_rej,
      num_trees             = sample_size,
      maxN                  = maxN_used,
      time                  = results$time,
      final_estep           = final_estep
    ))
  }

  for (i in seq_len(max_iter)) {
    maxN_used <- maxN
    results   <- run_em(cur_pars, maxN)

    if (is.null(results)) {
      fail_streak <- fail_streak + 1L
      n_failed    <- n_failed + 1L
      streak      <- 0L

      # Escalation: enlarge the attempt cap. The enlarged value is kept for
      # the rest of the run (never reset after a success, never reduced).
      maxN <- as.integer(max(maxN, min(2 * maxN, maxN_cap)))

      # First failure: retry at the last successful iterate with the larger
      # cap. From the second consecutive failure on, also move the sampling
      # point toward the box centre; prev_pars keeps the last estimate.
      if (fail_streak == 1L) {
        cur_pars <- prev_pars
        how <- "restarting from the last estimate"
      } else {
        cur_pars <- 0.8 * cur_pars + 0.2 * center_pars
        cur_pars <- pmax(cur_pars, lower_bound)
        cur_pars <- pmin(cur_pars, upper_bound)
        how <- "perturbed toward center"
      }

      if (verbose) message(sprintf(
        "Iteration %d: E-step failed (%d consecutive) - recovery: maxN=%d, %s",
        i, fail_streak, maxN, how
      ))

      if (fail_streak >= 8L) {
        stop_reason <- "e_step_failure"
        .mcem_warn_estep(brts, prev_pars, lower_bound, upper_bound, model, link)
        break
      }
    } else {
      fail_streak  <- 0L
      n_success    <- n_success + 1L
      last_results <- results

      new_pars  <- as.numeric(results$estimates)
      abs_step  <- max(abs(new_pars - cur_pars))
      delta_max <- rel_delta(new_pars, cur_pars)
      par_hist[[n_success]] <- new_pars
      # Relative displacement over the last `patience` completed iterations
      drift <- if (n_success > patience)
        rel_delta(new_pars, par_hist[[n_success - patience]]) else NA_real_
      prev_pars <- new_pars
      cur_pars  <- new_pars

      step <- trace_row(new_pars, results, delta_max, abs_step, drift,
                        maxN_used, final_estep = FALSE)
      mcem <- rbind(mcem, step)

      if (verbose) {
        rej_str <- if (step$n_rejected > 0L) sprintf("  rej=%d", step$n_rejected) else ""
        message(sprintf("Iteration %d: fhat=%.4f  delta=%.2e  step=%.2e  streak=%d/%d%s",
                        i, results$fhat, delta_max, abs_step, streak, patience, rej_str))
      }

      # Convergence: `patience` consecutive iterations with rel_delta < tol
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

    # Time budget check (counts failed iterations as well)
    if (!is.null(max_time)) {
      elapsed <- proc.time()[3] - t0_mcem
      if (elapsed > max_time) {
        stop_reason <- "time_budget"
        if (verbose) message(sprintf("Time budget reached (%.0fs > %ds)", elapsed, as.integer(max_time)))
        break
      }
    }
  }

  # Final E-step at the returned iterate, so that fhat, loglik_var and final_IS
  # describe the same point as `pars`. The M-step estimate of this call is
  # discarded. Recorded as the last trace row with final_estep = TRUE.
  final_estep <- FALSE
  if (!is.null(last_results) && stop_reason != "e_step_failure") {
    fin <- run_em(prev_pars, maxN)
    if (!is.null(fin)) {
      last_results <- fin
      final_estep  <- TRUE
      mcem <- rbind(mcem, trace_row(prev_pars, fin, NA_real_, NA_real_, NA_real_,
                                    maxN, final_estep = TRUE))
      if (verbose) message(sprintf("Final E-step at the returned parameters: fhat=%.4f",
                                   fin$fhat))
    } else if (verbose) {
      message("Final E-step at the returned parameters failed; ",
              "fhat is the last iteration's value.")
    }
  }

  # Bootstrap variance from the last E-step's IS weights
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
    pars        = prev_pars,
    iterations  = n_success,
    n_failed    = n_failed,
    stop_reason = stop_reason,
    final_estep = final_estep,
    maxN        = maxN,
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
      "  e.g. D or DD model with large tree and strongly negative covariate slope.\n",
      "  Suggestions:\n",
      "    1. Use link=\"exponential\" (ensures lambda > 0 everywhere).\n",
      "    2. Restrict bounds so the covariate slope cannot drive lambda to zero.\n",
      "    3. Use a simpler model (e.g. CR or DD) for this tree size.",
      call. = FALSE
    )
  } else {
    warning(
      "MCEM: 8 consecutive E-step failures (with adaptive maxN recovery); stopping. ",
      "Try widening bounds, increasing max_missing, or using link=\"exponential\".",
      call. = FALSE
    )
  }
}

#' MCEM with thinning augmentation.
#'
#' E-step and M-step of one iteration are a single \code{em_cpp()} call: the
#' augmented trees are drawn by the thinning sampler, up to \code{maxN}
#' attempts per E-step.
#'
#' Stopping rule, the same statistic \code{.mcem_bdi} uses.  With theta_k the
#' M-step output of iteration k,
#' \deqn{\delta_k = \max_j |\theta_k[j] - \theta_{k-1}[j]| / \max(|\theta_{k-1}[j]|, \epsilon)}
#' and the run stops as "converged" once \eqn{\delta_k < tol} for
#' \code{patience} consecutive iterations.  The search box does not enter, so
#' the same fit stops at the same point whatever bounds the user passes.  The
#' floor \eqn{\epsilon} comes from \code{.rel_floor(brts, link, rel_floor)}
#' and carries the units of the parameters (\code{rel_floor / crown_age} for
#' the rate-valued links, \code{rel_floor} for the log-rates of the
#' exponential link), so the statistic is unchanged when the tree and the
#' rates are expressed in another unit of time.  An iteration whose M-step
#' returned its starting point unchanged carries no information about
#' stability and resets the streak (trace column \code{m_moved}).  The trace
#' also records the absolute step \code{abs_step} and \code{drift}, the same
#' relative displacement measured over the last \code{patience} iterations.
#'
#' Trace.  Row k of \code{mcem} pairs theta_k with the E-step value
#' fhat(theta_\{k-1\}) it was computed from, and with the rejection counters
#' of that E-step.  The likelihood of the returned point is a separate E-step
#' at the final iterate, appended as the last row (\code{m_step = FALSE}, step
#' columns \code{NA}); only that row is evaluated at \code{pars}, and it is
#' not counted in \code{iterations}.  It is attempted whenever at least one
#' iteration succeeded, including after a \code{"e_step_failure"} or
#' \code{"time_budget"} stop, and may itself fail, in which case there is no
#' such row and \code{loglik}, \code{loglik_var} and \code{final_IS} report
#' nothing rather than the previous iterate's values.
#'
#' Rejection columns.  \code{rejected_errors} (unhandled exception),
#' \code{rejected_overruns} (\code{max_missing} exceeded),
#' \code{rejected_lambda} (\code{max_lambda} exceeded) and
#' \code{rejected_nonfinite} (log weight \code{+Inf} or \code{NaN}) are the
#' four disjoint ways a draw is discarded; \code{n_rejected} is their sum, the
#' quantity \code{final_IS$n_rejected} reports, and \code{rejected} is that
#' same total under the name the BDI trace uses.  \code{rejected_zero_weights}
#' is counted apart: those draws completed, have weight zero and stay in the
#' \code{fhat} denominator.
#'
#' @return A list: \code{mcem} (trace, one row per completed iteration plus the
#'   final E-step), \code{pars}, \code{iterations} (number of completed E+M
#'   iterations, not loop passes and not \code{nrow(mcem)}), \code{n_failed}
#'   (E-step failures over the whole run), \code{stop_reason},
#'   \code{final_estep} (did the final E-step at \code{pars} succeed),
#'   \code{maxN} (the ratcheted attempt cap: doubled on every E-step failure
#'   up to 50000, never reset after a success), \code{loglik} (fhat at
#'   \code{pars}, \code{NA} when the final E-step failed), \code{loglik_var},
#'   \code{final_IS}.
#' @keywords internal
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
                      rho = 1.0,
                      rel_floor = 1e-2,
                      stop_rule = c("rel_change", "mc_error"),
                      mc_batches = 5L,
                      mc_z = 1.0,
                      mc_grow = 1.5,
                      max_draws = NULL) {
  stop_rule <- match.arg(stop_rule)
  if (inherits(brts, "phylo")) {
    # .extract_brts also carries the tip starts M and D are measured from.
    brts <- .extract_brts(brts)
  }
  if (!is.numeric(brts)) stop("`brts` must be numeric or a `phylo` object.")
  if (!is.numeric(pars)) stop("`pars` must be numeric.")
  if (!is.null(conditional) && !is.function(conditional))
    stop("`conditional` must be a function or NULL.")

  # Empty when only branching times were passed: M and D then fall back to the
  # crown-age convention, as documented on estimate_rates().
  parent_tip_start <- .pts(brts)
  parent_id        <- .pid(brts)

  # Convergence metric: max_j |theta_k,j - theta_{k-1,j}| / max(|theta_{k-1,j}|, floor).
  # The bound box does not enter; the floor keeps parameters at or near zero
  # from giving an unbounded ratio and carries their units, so the metric is
  # the same under a change of the time unit.
  floor_val <- .rel_floor(brts, link, rel_floor)
  rel_delta <- function(new, old) .rel_change(new, old, floor_val)

  streak       <- 0L      # consecutive iterations with rel_delta < tol
  fail_streak  <- 0L      # consecutive E-step failures
  n_failed     <- 0L      # E-step failures over the whole run
  n_success    <- 0L      # completed EM iterations
  prev_pars    <- pars    # last successful M-step estimate (init until the first success)
  cur_pars     <- pars    # point the next E-step samples at and the M-step starts from
  prev_iterate <- pars    # the iterate the "mc_error" rule measures the step from
  if (is.null(max_draws) || !is.finite(max_draws)) max_draws <- as.integer(8L * sample_size)
  max_draws    <- as.integer(max(max_draws, sample_size))
  par_hist     <- list()  # successful iterates, for the windowed drift
  mcem         <- NULL
  had_success  <- FALSE   # at least one iteration completed, so `pars` is an estimate
  stop_reason  <- "max_iter"
  t0_mcem      <- proc.time()[3]
  center_pars  <- (lower_bound + upper_bound) / 2
  maxN_cap     <- 50000L

  # Monte Carlo error of one M-step estimate, by batch means: the E-step's
  # draws are split into `mc_batches` blocks, the M-step is re-run on each
  # from the full-sample estimate, and the standard error of the full-sample
  # estimate is the between-batch standard deviation over sqrt(B).  This is
  # what the "mc_error" rule compares the EM step against, so that a step is
  # called small only when it is small relative to the noise that produced it
  # (Booth & Hobert 1999; Caffo, Jank & Jones 2005).
  mc_se <- function(results, at) {
    trees <- results$trees
    w     <- as.numeric(results$weights)
    B     <- as.integer(mc_batches)
    if (is.null(trees) || !length(trees) || length(w) != length(trees) || B < 2L) return(NULL)
    n <- length(trees)
    if (n < 2L * B) return(NULL)
    idx <- split(seq_len(n), rep(seq_len(B), length.out = n))
    est <- lapply(idx, function(ii) {
      if (!any(w[ii] > 0)) return(NULL)
      e <- list(trees = trees[ii], weights = w[ii], rejected = 0L, rejected_overruns = 0L,
                rejected_lambda = 0L, rejected_zero_weights = 0L, time = 0, fhat = 0)
      m <- tryCatch(m_cpp(e, at, "rpd5c", lower_bound, upper_bound, xtol, num_threads,
                          model = as.integer(model), link = as.integer(link), rho = as.numeric(rho),
                          rconditional = conditional), error = function(err) NULL)
      if (is.null(m)) NULL else as.numeric(m$estimates)[seq_along(at)]
    })
    est <- do.call(rbind, Filter(Negate(is.null), est))
    if (is.null(est) || nrow(est) < 2L) return(NULL)
    apply(est, 2, stats::sd) / sqrt(nrow(est))
  }

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
             copy_trees = (stop_rule == "mc_error"),
             model = as.integer(model),
             link = as.integer(link),
             rho = as.numeric(rho),
             rconditional = conditional,
             parent_tip_start = parent_tip_start,
             parent_id = parent_id),
      error = function(e) NULL
    )
  }

  # The four disjoint ways the E-step discards a draw. Zero-weight draws are
  # not among them: they completed and stay in the fhat denominator.
  n_rejected_of <- function(results) {
    .n0(results$rejected) +
      .n0(results$rejected_overruns) +
      .n0(results$rejected_lambda) +
      .n0(results$rejected_nonfinite)
  }

  # One trace row: the iterate, the E-step value, the step metrics and the
  # rejection counters of that E-step. Step columns are NA for the final
  # E-step (no M-step), which is flagged `m_step = FALSE`.
  trace_row <- function(p, results, delta_max, abs_step, drift, maxN_used,
                        m_step = TRUE, m_moved = NA) {
    par_df <- as.data.frame(as.list(stats::setNames(p, paste0("par", seq_along(p)))))
    n_rej <- n_rejected_of(results)
    lw <- results$logf - results$logg
    cbind(par_df, data.frame(
      fhat                  = results$fhat,
      delta_max             = delta_max,
      abs_step              = abs_step,
      drift                 = drift,
      m_step                = m_step,
      m_moved               = m_moved,
      rejected              = n_rej,
      rejected_errors       = .n0(results$rejected),
      rejected_overruns     = .n0(results$rejected_overruns),
      rejected_lambda       = .n0(results$rejected_lambda),
      rejected_nonfinite    = .n0(results$rejected_nonfinite),
      rejected_zero_weights = .n0(results$rejected_zero_weights),
      n_rejected            = n_rej,
      num_trees             = if (is.null(results$num_trees)) sample_size
                              else .n0(results$num_trees),
      maxN                  = maxN_used,
      ESS                   = .ess_from_lw(lw),
      # "mc_error" only: the step in units of its own Monte Carlo standard
      # error, and the number of draws the E-step used at that iteration.
      mc_z                  = NA_real_,
      sample_size           = sample_size,
      time                  = results$time
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
      had_success  <- TRUE

      new_pars  <- as.numeric(results$estimates)
      m_moved   <- any(new_pars != cur_pars)
      abs_step  <- max(abs(new_pars - cur_pars))
      delta_max <- rel_delta(new_pars, cur_pars)
      par_hist[[n_success]] <- new_pars
      # Relative displacement over the last `patience` completed iterations
      drift <- if (n_success > patience)
        rel_delta(new_pars, par_hist[[n_success - patience]]) else NA_real_
      prev_pars <- new_pars
      cur_pars  <- new_pars

      step <- trace_row(new_pars, results, delta_max, abs_step, drift,
                        maxN_used, m_step = TRUE, m_moved = m_moved)
      mcem <- rbind(mcem, step)

      if (verbose) {
        rej_str <- if (step$n_rejected > 0L) sprintf("  rej=%d", step$n_rejected) else ""
        message(sprintf("Iteration %d: fhat=%.4f  delta=%.2e  step=%.2e  streak=%d/%d%s%s",
                        i, results$fhat, delta_max, abs_step, streak, patience, rej_str,
                        if (m_moved) "" else "  (M-step returned its start)"))
      }

      # Convergence.
      #
      # "rel_change" (the package's original rule): `patience` consecutive
      # iterations with rel_delta < tol.  It compares the EM step with a fixed
      # tolerance and never looks at the Monte Carlo noise that produced the
      # step, so at a fixed number of draws it stops when three successive
      # random steps happen to be small.
      #
      # "mc_error": the step is compared with its own Monte Carlo standard
      # error (mc_se above).  A step smaller than `mc_z` standard errors in
      # every coordinate is indistinguishable from noise; when that happens
      # the number of draws is raised by `mc_grow` (so the next step is
      # measured more finely) and the run stops only once the draws have
      # reached `max_draws` and the step has stayed inside the noise for
      # `patience` iterations.  The estimate returned is then the mean of the
      # last `patience` iterates, which removes part of the remaining Monte
      # Carlo error.  Booth & Hobert (1999) and Caffo, Jank & Jones (2005).
      #
      # An M-step that returned its starting point unchanged says nothing
      # about stability (the objective may have been undefined there), so it
      # does not count toward patience under either rule.
      if (stop_rule == "rel_change") {
        if (!m_moved) {
          streak <- 0L
        } else if (delta_max < tol) {
          streak <- streak + 1L
          if (streak >= patience) {
            stop_reason <- "converged"
            break
          }
        } else {
          streak <- 0L
        }
      } else {
        se <- mc_se(results, new_pars)
        z  <- if (is.null(se)) NA_real_ else max(abs(new_pars - prev_iterate) / pmax(se, .Machine$double.eps))
        mcem$mc_z[nrow(mcem)] <- z
        mcem$sample_size[nrow(mcem)] <- sample_size
        if (!m_moved || is.na(z)) {
          streak <- 0L
        } else if (z < mc_z) {
          streak <- streak + 1L
          if (sample_size < max_draws) {
            sample_size <- min(as.integer(ceiling(mc_grow * sample_size)), max_draws)
            maxN <- max(maxN, 10L * sample_size)
            streak <- 0L
            if (verbose) message(sprintf("  step within Monte Carlo error (z = %.2f): draws -> %d", z, sample_size))
          } else if (streak >= patience) {
            stop_reason <- "converged"
            break
          }
        } else {
          streak <- 0L
        }
      }
      prev_iterate <- new_pars
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

  # Under "mc_error" the returned estimate is the mean of the last `patience`
  # successful iterates rather than the last one: at the fixed point the
  # iterates fluctuate around the maximiser with the Monte Carlo error the
  # rule measured, and averaging them removes part of it (Caffo, Jank & Jones
  # 2005).  The final E-step below is then run at that average, so the
  # reported likelihood describes the point that is returned.
  if (stop_rule == "mc_error" && had_success && n_success >= 2L) {
    take <- par_hist[seq(max(1L, n_success - patience + 1L), n_success)]
    prev_pars <- colMeans(do.call(rbind, take))
  }

  # Final E-step at the returned iterate, so that fhat, loglik_var and final_IS
  # describe the same point as `pars`. The M-step estimate of this call is
  # discarded. Recorded as the last trace row with m_step = FALSE. It is
  # attempted whenever an iteration succeeded, whatever the stop reason: the
  # alternative is to report the previous iterate's E-step as the likelihood
  # of `pars`. When it fails, loglik, loglik_var and final_IS report nothing.
  final_estep <- FALSE
  fin <- NULL
  if (had_success) {
    fin <- run_em(prev_pars, maxN)
    if (!is.null(fin)) {
      final_estep <- TRUE
      mcem <- rbind(mcem, trace_row(prev_pars, fin, NA_real_, NA_real_, NA_real_,
                                    maxN, m_step = FALSE, m_moved = NA))
      if (verbose) message(sprintf("Final E-step at the returned parameters: fhat=%.4f",
                                   fin$fhat))
    } else if (verbose) {
      message("Final E-step at the returned parameters failed; loglik is NA.")
    }
  }

  loglik     <- NA_real_
  loglik_var <- NA_real_
  final_IS   <- NULL
  if (!is.null(fin)) {
    loglik <- fin$fhat

    # Bootstrap variance from the final E-step's IS weights
    if (length(fin$logf) >= 2L && all(is.finite(fin$logf))) {
      loglik_var <- .bootstrap_fhat_var(fin$logf, fin$logg, K = 2L, B = 200L)
    }

    # IS components at `pars`, for diagnostics
    if (length(fin$logf) > 0L) {
      lw <- fin$logf - fin$logg
      final_IS <- list(
        logf  = fin$logf,
        logg  = fin$logg,
        lw    = lw,
        fhat  = .is_fhat(fin$logf, fin$logg,
                         n_zero_weight = .n0(fin$rejected_zero_weights)),
        ESS   = .ess_from_lw(lw),
        n_rejected = n_rejected_of(fin),
        rejected_zero_weights = .n0(fin$rejected_zero_weights)
      )
    }
  }

  list(
    mcem        = mcem,
    pars        = prev_pars,
    iterations  = n_success,
    n_failed    = n_failed,
    stop_reason = stop_reason,
    final_estep = final_estep,
    maxN        = maxN,
    stop_rule   = stop_rule,
    sample_size = sample_size,
    loglik      = loglik,
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
           link         = as.integer(link),
           parent_tip_start = .pts(brts),
           parent_id    = .pid(brts))
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

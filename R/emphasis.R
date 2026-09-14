#' MCEM with thinning augmentation.
#'
#' One iteration is \code{augment_trees()} followed by \code{m_cpp()}: the
#' augmented trees are drawn by the thinning sampler, up to \code{maxN}
#' attempts per E-step, and the M-step then runs on the weights
#' \eqn{\exp(\log f - \log g - \max)} that the C++ E-step would have handed
#' it.  The two calls are what \code{em_cpp()} does internally
#' (\code{src/mcem.cpp} is \code{E_step} followed by \code{M_step} with the
#' same arguments), so the split is scoring-identical and no file under
#' \code{src/} changed; it is what lets an undecided iteration enlarge its
#' sample without redrawing, and lets this driver run the same bootstrap the
#' BDI driver runs.
#'
#' Stopping rule: \code{.ascent_controller}, the same object
#' \code{.mcem_bdi} instantiates, from the same defaults resolved once in
#' \code{.run_mcem}.  Each iteration is scored on the gain it bought in nats
#' of the EM surrogate against the Monte Carlo error of that gain.  Growth is
#' append-only: an undecided iteration does not move the iterate, so only the
#' shortfall is drawn and it pools with the draws in hand.  \code{tol},
#' \code{patience} and \code{.rel_floor} keep their meaning and are what the
#' fallback path, and \code{stop = "step"}, use.  The trace also records the
#' absolute step \code{abs_step} and \code{drift}, the same relative
#' displacement measured over the last \code{patience} iterations.
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
                      stop_opts = NULL) {
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

  # Convergence metric: max_j |theta_k,j - theta_{k-1,j}| / max(|theta_{k-1,j}|, floor).
  # The bound box does not enter; the floor keeps parameters at or near zero
  # from giving an unbounded ratio and carries their units, so the metric is
  # the same under a change of the time unit.
  floor_val <- .rel_floor(brts, link, rel_floor)
  rel_delta <- function(new, old) .rel_change(new, old, floor_val)

  fail_streak  <- 0L      # consecutive E-step failures
  n_failed     <- 0L      # E-step failures over the whole run
  n_success    <- 0L      # completed EM iterations
  prev_pars    <- pars    # last successful M-step estimate (init until the first success)
  cur_pars     <- pars    # point the next E-step samples at and the M-step starts from
  par_hist     <- list()  # successful iterates, for the windowed drift
  mcem         <- NULL
  pool         <- NULL    # draws held across a growth step
  had_success  <- FALSE   # at least one iteration completed, so `pars` is an estimate
  stop_reason  <- "max_iter"
  t0_mcem      <- proc.time()[3]
  center_pars  <- (lower_bound + upper_bound) / 2
  maxN_cap     <- 50000L

  # ── E-step: the thinning sampler alone ──
  run_e <- function(at, maxN_now, n_want) {
    e <- tryCatch(
      augment_trees(brts = brts,
                    pars = at,
                    sample_size = max(1L, as.integer(n_want)),
                    maxN = as.integer(maxN_now),
                    max_missing = max_missing,
                    max_lambda = 1e6,
                    num_threads = num_threads,
                    model = as.integer(model),
                    link = as.integer(link),
                    rho = as.numeric(rho),
                    parent_tip_start = parent_tip_start),
      error = function(e) NULL)
    if (is.null(e) || length(e$trees) == 0L) return(NULL)
    e
  }

  # Pooling two E-steps drawn at the same theta: the rejection counters add,
  # and fhat is rebuilt over the pool rather than averaged.
  e_pool <- function(a, b) {
    if (is.null(a)) return(b)
    if (is.null(b)) return(a)
    add <- function(k) .n0(a[[k]]) + .n0(b[[k]])
    out <- list(trees = c(a$trees, b$trees),
                logf = c(a$logf, b$logf), logg = c(a$logg, b$logg),
                rejected = add("rejected"),
                rejected_overruns = add("rejected_overruns"),
                rejected_lambda = add("rejected_lambda"),
                rejected_zero_weights = add("rejected_zero_weights"),
                rejected_nonfinite = add("rejected_nonfinite"),
                envelope_violations = add("envelope_violations"),
                time = add("time"))
    out$num_trees <- length(out$trees)
    out
  }

  # ── M-step: the same optimisation em_cpp() runs, on the weights the C++
  # E-step would have handed it (exp(log w - max log w)).  Verified against
  # em_cpp: identical() on the estimates, on logf, and on fhat.
  run_m <- function(e, at) {
    lw <- e$logf - e$logg
    w  <- exp(lw - max(lw[is.finite(lw)]))
    w[!is.finite(w)] <- 0
    es <- list(trees = e$trees, weights = w,
               rejected = .n0(e$rejected),
               rejected_overruns = .n0(e$rejected_overruns),
               rejected_lambda = .n0(e$rejected_lambda),
               rejected_zero_weights = .n0(e$rejected_zero_weights),
               time = .n0(e$time),
               fhat = .is_fhat(e$logf, e$logg,
                               n_zero_weight = .n0(e$rejected_zero_weights)))
    m <- tryCatch(
      m_cpp(e_step = es, init_pars = at, plugin = "rpd1",
            lower_bound = lower_bound, upper_bound = upper_bound,
            xtol_rel = xtol, num_threads = num_threads,
            model = as.integer(model), link = as.integer(link),
            rho = as.numeric(rho), rconditional = conditional),
      error = function(e) NULL)
    if (is.null(m)) return(NULL)
    list(estimates = m$estimates, weights = w, fhat = es$fhat, m_time = m$time)
  }

  # ── the stopping rule, shared with .mcem_bdi ──
  opts <- if (is.null(stop_opts))
    .ascent_opts(list(tol = tol, patience = patience, sample_size = sample_size,
                      max_iter = max_iter), lower_bound, upper_bound)
    else stop_opts
  opts$tol <- tol; opts$patience <- patience
  opts$sample_size <- as.integer(sample_size); opts$max_iter <- as.integer(max_iter)

  ctl_trees <- NULL
  refit <- function(w, init) {
    m <- tryCatch(
      m_cpp(e_step = list(trees = ctl_trees, weights = w, rejected = 0L,
                          rejected_overruns = 0L, rejected_lambda = 0L,
                          rejected_zero_weights = 0L, time = 0, fhat = 0),
            init_pars = init, plugin = "rpd1",
            lower_bound = lower_bound, upper_bound = upper_bound,
            xtol_rel = xtol, num_threads = num_threads,
            model = as.integer(model), link = as.integer(link),
            rho = as.numeric(rho), rconditional = conditional),
      error = function(e) NULL)
    if (is.null(m)) NULL else as.numeric(m$estimates)
  }
  score <- function(theta, trees) eval_logf(theta, trees,
                                            model = as.integer(model),
                                            link = as.integer(link),
                                            rho = as.numeric(rho))$logf
  ctl <- .ascent_controller(opts, refit, score)

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
                        m_step = TRUE, m_moved = NA, stop_cols = NULL) {
    par_df <- as.data.frame(as.list(stats::setNames(p, paste0("par", seq_along(p)))))
    n_rej <- n_rejected_of(results)
    lw <- results$logf - results$logg
    cbind(par_df, .ascent_trace_cols(stop_cols), data.frame(
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
      time                  = results$time
    ))
  }

  for (i in seq_len(max_iter)) {
    maxN_used <- maxN
    # On a growth step the iterate did not move, so only the shortfall is
    # drawn and it pools with the draws already held.
    have <- if (is.null(pool)) 0L else length(pool$trees)
    want <- max(0L, as.integer(ctl$S) - have)
    e_new <- if (want > 0L) run_e(cur_pars, maxN, want) else NULL
    if (!is.null(e_new)) ctl$add_draws(length(e_new$trees))
    e_raw <- if (want > 0L) {
      if (is.null(e_new)) NULL else e_pool(pool, e_new)
    } else pool
    ctl_trees <- if (is.null(e_raw)) NULL else e_raw$trees
    m_out <- if (is.null(e_raw)) NULL else run_m(e_raw, cur_pars)
    results <- if (is.null(e_raw) || is.null(m_out)) NULL else
      utils::modifyList(e_raw, list(estimates = m_out$estimates,
                                    fhat = m_out$fhat,
                                    time = .n0(e_raw$time) + .n0(m_out$m_time)))

    if (is.null(results)) {
      pool <- NULL
      ctl$reset()
      fail_streak <- fail_streak + 1L
      n_failed    <- n_failed + 1L

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

      # ── stopping rule ──
      v <- ctl$decide(theta_old = cur_pars, theta_new = new_pars,
                      m_moved = m_moved, trees = e_raw$trees,
                      w_raw = m_out$weights, logf_old = e_raw$logf,
                      delta_max = delta_max, abs_step = abs_step)
      v$draws <- ctl$draws

      if (isTRUE(v$accept)) {
        prev_pars <- new_pars
        cur_pars  <- new_pars
        pool      <- NULL
      } else {
        # Undecided on precision: the step is not taken and the draws are
        # kept, so the enlarged sample pools with them at the same theta.
        pool <- e_raw
      }
      par_hist[[n_success]] <- cur_pars
      # Relative displacement over the last `patience` completed iterations
      drift <- if (n_success > patience)
        rel_delta(cur_pars, par_hist[[n_success - patience]]) else NA_real_

      step <- trace_row(cur_pars, results, delta_max, abs_step, drift,
                        maxN_used, m_step = TRUE, m_moved = m_moved,
                        stop_cols = v)
      mcem <- rbind(mcem, step)

      if (verbose) {
        rej_str <- if (step$n_rejected > 0L) sprintf("  rej=%d", step$n_rejected) else ""
        message(sprintf(
          "Iteration %d: fhat=%.4f  dQ=%.3e  U=%.3e  case=%s  S=%d  ESS=%.0f%s%s",
          i, results$fhat, v$dQ, v$U, format(v$case), as.integer(ctl$S),
          v$ESS, rej_str,
          if (isTRUE(v$accept)) "" else "  (held: growing the sample)"))
      }

      if (isTRUE(v$stop)) { stop_reason <- v$stop_reason; break }
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
  # discarded. Recorded as the last trace row with m_step = FALSE. It is
  # attempted whenever an iteration succeeded, whatever the stop reason: the
  # alternative is to report the previous iterate's E-step as the likelihood
  # of `pars`. When it fails, loglik, loglik_var and final_IS report nothing.
  final_estep <- FALSE
  fin <- NULL
  if (had_success) {
    fin <- run_e(prev_pars, maxN, sample_size)
    if (!is.null(fin)) {
      fin$fhat <- .is_fhat(fin$logf, fin$logg,
                           n_zero_weight = .n0(fin$rejected_zero_weights))
      final_estep <- TRUE
      ctl$add_draws(length(fin$trees))
      mcem <- rbind(mcem, trace_row(prev_pars, fin, NA_real_, NA_real_, NA_real_,
                                    maxN, m_step = FALSE, m_moved = NA,
                                    stop_cols = list(S = as.numeric(ctl$S),
                                                     draws = ctl$draws)))
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
    loglik      = loglik,
    loglik_var  = loglik_var,
    final_IS    = final_IS,
    total_draws = ctl$draws,
    sample_size_final = as.integer(ctl$S)
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
           parent_tip_start = .pts(brts))
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

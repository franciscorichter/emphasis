# ---------------------------------------------------------------------------
# 03-worker.R — the function executed in each child process.  One fit per
# process (callr::r_bg), so a C++ abort or a hang costs one job, not the run.
#
# val_run_job() never throws: every path returns a row with
#   outcome in {ok, error, timeout, crash}
# and the full schema below.  "timeout" and "crash" are written by the driver.
# ---------------------------------------------------------------------------

val_job_env <- function() {
  Sys.setenv(OMP_NUM_THREADS = "1", RCPP_PARALLEL_NUM_THREADS = "1",
             TBB_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
             OPENBLAS_NUM_THREADS = "1")
  invisible(TRUE)
}

# summary of a final_IS block, tolerant of fields a pre-wave-1 build lacks
.val_is_summary <- function(fi) {
  g <- function(k) if (!is.null(fi[[k]])) as.numeric(fi[[k]])[1] else NA_real_
  lw <- fi$lw
  ess <- if (!is.null(fi$ESS)) as.numeric(fi$ESS)[1] else
    if (!is.null(lw) && any(is.finite(lw))) {
      w <- exp(lw[is.finite(lw)] - max(lw[is.finite(lw)]))
      sum(w)^2 / sum(w^2)
    } else NA_real_
  list(ESS = ess, n_rejected = g("n_rejected"),
       rejected_zero_weights = g("rejected_zero_weights"),
       acc = g("acc"),
       n_nonfinite = if (!is.null(lw)) sum(!is.finite(lw)) else NA_real_,
       n_lw = if (!is.null(lw)) length(lw) else NA_real_,
       sd_lw = if (!is.null(lw) && sum(is.finite(lw)) > 1)
         stats::sd(lw[is.finite(lw)]) else NA_real_)
}

val_run_job <- function(job, lib, trees_file, refs_file, out_file) {
  val_job_env()
  t0 <- Sys.time()
  row <- c(job, list(outcome = "error", error = NA_character_,
                     elapsed = NA_real_, host = Sys.info()[["nodename"]]))
  res <- tryCatch({
    val_load_emphasis(lib)
    val_load_refs()
    row$build <- val_build_fingerprint()
    TR  <- readRDS(trees_file)$trees
    REF <- readRDS(refs_file)$refs
    tt  <- TR[[job$tree_id]]
    rf  <- REF[[job$tree_id]]
    out <- switch(job$kind,
      cr_fit   = .val_fit(job, tt, rf, "cr"),
      dd_fit   = .val_fit(job, tt, rf, "dd"),
      fhat     = .val_fhat(job, tt, rf),
      pipeline = .val_pipeline(job, tt, rf),
      init     = .val_init_arm(job, tt, rf),
      stop("unknown job kind: ", job$kind))
    row <- c(row, out)
    row$outcome <- "ok"
    row
  }, error = function(e) {
    row$error <- conditionMessage(e)
    row
  })
  res$elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  saveRDS(res, out_file)
  invisible(NULL)
}

# --- mcem fits -------------------------------------------------------------
.val_fit <- function(job, tt, rf, model) {
  brts <- tt$brts
  lam_gen <- if (model == "cr") tt$lambda_gen else tt$lambda0_gen
  box <- if (model == "cr") val_box_cr(lam_gen, job$box_scale)
         else val_box_dd(tt$lambda0_gen, tt$mu0_gen, tt$K_gen)
  if (model == "dd" && job$box_scale != 1) {
    box$upper[1] <- box$upper[1] * job$box_scale
    box$lower[2] <- box$lower[2] * job$box_scale
    box$upper[3] <- box$upper[3] * job$box_scale
  }
  init <- if (model == "cr") val_init_cr(job$init, lam_gen, rf$mle)
          else val_init_dd(job$init, tt$lambda0_gen, tt$mu0_gen, tt$K_gen, rf$mle)
  init <- pmin(pmax(init, box$lower), box$upper)

  ctrl <- list(lower_bound = box$lower, upper_bound = box$upper,
               sampling = job$sampler, sample_size = as.integer(job$N),
               maxN = as.integer(20 * job$N), max_iter = 400L,
               max_missing = 1e4, num_threads = 1L,
               max_time = max(30, job$timeout_s - 30), verbose = FALSE)

  set.seed(job$seed)
  fit <- emphasis::estimate_rates(brts, method = "mcem", model = model,
                                  init_pars = init, control = ctrl,
                                  link = "linear", cond = NULL)

  pars <- as.numeric(fit$pars)
  lx <- if (model == "dd") rf$lx else NULL
  ll_hat  <- ll_exact(pars, brts, model = model, lx = lx)
  ll_mle  <- rf$mle_loglik
  trace   <- fit$details$mcem
  prev    <- NULL
  if (!is.null(trace) && nrow(trace) >= 2L) {
    k <- nrow(trace)
    p8 <- as.numeric(trace[k - 1L, paste0("par", 1:8)])
    prev <- if (model == "cr") p8[c(1, 5)] else p8[c(1, 2, 5, 6)]
  }
  ll_prev <- if (!is.null(prev)) ll_exact(prev, brts, model = model, lx = lx) else NA_real_

  se <- rf$se
  npar <- length(pars)
  z <- rep(NA_real_, npar)
  for (j in seq_len(npar))
    if (is.finite(se[j]) && se[j] > 0) z[j] <- (pars[j] - rf$mle[j]) / se[j]
  if (model == "cr" && isTRUE(rf$boundary_mu0)) z[2] <- NA_real_

  fi <- .val_is_summary(fit$details$final_IS)
  dd_map <- if (model == "dd") val_emphasis_to_dd(pars) else NULL

  list(
    pars = pars, init_pars = init,
    lower_bound = box$lower, upper_bound = box$upper,
    at_bound = any(abs(pars - box$lower) < 1e-9 | abs(pars - box$upper) < 1e-9),
    loglik = fit$loglik, loglik_var = fit$loglik_var,
    AIC = fit$AIC, n_pars = fit$n_pars,
    stop_reason = fit$stop_reason %||% fit$details$stop_reason %||% NA_character_,
    iterations = fit$iterations %||% fit$details$iterations %||% NA_integer_,
    ll_exact_hat = ll_hat, ll_exact_mle = ll_mle,
    delta_ell = ll_hat - ll_mle,                 # <= 0 up to reference tolerance
    ll_exact_prev = ll_prev,
    e_at_hat  = fit$loglik - ll_hat,             # reported vs exact at theta_K
    e_at_prev = fit$loglik - ll_prev,            # ... at theta_{K-1} (H20)
    z = z, mle = rf$mle, se = se,
    dist_init_se = if (all(is.finite(se[seq_len(npar)]) & se[seq_len(npar)] > 0))
      max(abs((init - rf$mle) / se[seq_len(npar)])) else NA_real_,
    dd_map_status = if (!is.null(dd_map)) dd_map$status else NA_character_,
    K_hat = if (!is.null(dd_map)) dd_map$K else NA_real_,
    final_IS = fi, trace = trace)
}

# --- fixed-theta importance-sampling grid ---------------------------------
.val_fhat <- function(job, tt, rf) {
  brts <- tt$brts
  model <- if (tt$kind == "cr") "cr" else "dd"
  mb <- if (model == "cr") c(0L, 0L, 0L) else c(1L, 0L, 0L)
  grid <- rf$theta_grid
  lx <- if (model == "dd") rf$lx else NULL
  N <- as.integer(job$N)
  reps <- as.integer(job$reps %||% 1L)
  rows <- list()
  for (i in seq_len(nrow(grid))) {
    th <- as.numeric(grid[i, ])
    pars <- if (model == "cr") c(th[1], th[2]) else th
    p8 <- if (model == "cr") c(th[1], 0, 0, 0, th[2], 0, 0, 0)
          else c(th[1], th[2], 0, 0, th[3], th[4], 0, 0)
    ll0 <- ll_exact(pars, brts, model = model, lx = lx)
    if (job$sampler == "bdi") {
      for (r in seq_len(reps)) {
        set.seed(job$seed + 1000L * i + r)
        a <- tryCatch(emphasis:::.augment_tree_bdi(brts, p8, model_bin = mb,
                                                   sample_size = N, link = 0L,
                                                   rho = 1),
                      error = function(e) e)
        if (inherits(a, "error")) {
          rows[[length(rows) + 1L]] <- data.frame(
            point = i, rep = r, sampler = "bdi", ll_exact = ll0,
            fhat = NA_real_, g = NA_real_, ess = NA_real_, sd_lw = NA_real_,
            n_valid = NA_real_, acc = NA_real_, n_nonfinite = NA_real_,
            err = conditionMessage(a), stringsAsFactors = FALSE)
        } else {
          lw <- a$weights
          ess <- if (any(is.finite(lw))) {
            w <- exp(lw[is.finite(lw)] - max(lw[is.finite(lw)]))
            sum(w)^2 / sum(w^2) } else NA_real_
          rows[[length(rows) + 1L]] <- data.frame(
            point = i, rep = r, sampler = "bdi", ll_exact = ll0,
            fhat = a$fhat, g = a$fhat - ll0, ess = ess,
            sd_lw = if (sum(is.finite(lw)) > 1) stats::sd(lw[is.finite(lw)]) else NA_real_,
            n_valid = a$n_valid %||% length(lw),
            acc = a$acc %||% NA_real_,
            n_nonfinite = a$n_nonfinite %||% NA_real_,
            err = NA_character_, stringsAsFactors = FALSE)
        }
      }
    } else {
      for (r in seq_len(reps)) {
        set.seed(job$seed + 1000L * i + r)
        a <- tryCatch(emphasis:::augment_trees(brts, p8, N, 20L * N, 1e4, 1e6,
                                               1L, model = mb, link = 0L, rho = 1),
                      error = function(e) e)
        if (inherits(a, "error")) {
          rows[[length(rows) + 1L]] <- data.frame(
            point = i, rep = r, sampler = "dynamic_fresh", ll_exact = ll0,
            fhat = NA_real_, g = NA_real_, ess = NA_real_, sd_lw = NA_real_,
            n_valid = NA_real_, acc = NA_real_, n_nonfinite = NA_real_,
            err = conditionMessage(a), stringsAsFactors = FALSE)
          next
        }
        ev <- emphasis:::eval_logf(p8, a$trees, model = mb, link = 0L, rho = 1)
        lw <- ev$logf - a$logg
        fh <- emphasis:::.is_fhat(ev$logf, a$logg,
                                 n_zero_weight = a$rejected_zero_weights)
        fin <- is.finite(lw)
        ess <- if (any(fin)) { w <- exp(lw[fin] - max(lw[fin]))
                               sum(w)^2 / sum(w^2) } else NA_real_
        rows[[length(rows) + 1L]] <- data.frame(
          point = i, rep = r, sampler = "dynamic_fresh", ll_exact = ll0,
          fhat = fh, g = fh - ll0, ess = ess,
          sd_lw = if (sum(fin) > 1) stats::sd(lw[fin]) else NA_real_,
          n_valid = length(lw), acc = NA_real_,
          n_nonfinite = sum(!fin), err = NA_character_, stringsAsFactors = FALSE)
      }
    }
  }
  list(grid = grid, fhat_rows = do.call(rbind, rows),
       se_scale = rf$se, mle = rf$mle)
}

# --- initialiser arm (I1/I2/I3) --------------------------------------------
# Does the cross-entropy search earn its place?  The three arms share ONE
# auto_bounds box, computed once per tree and cached, so the box is not a
# confounder: they differ only in what they hand to MCEM.
#
#   I1  cem    the pipeline's own route: GAM surface, then CEM inside the box
#   I2  naive  the box midpoint, which is what a user with no search would use
#   I3  gam    the GAM stage's optimum, with no CEM -- isolates whether CEM
#              adds anything over the surface it is initialised from
#
# Two deficits are recorded for each: the STARTING point's own deficit against
# the exact MLE (how good is the initialiser before MCEM touches it) and the
# final deficit after MCEM.  Cost is the E-step draws each initialiser spends.
.val_init_arm <- function(job, tt, rf) {
  brts  <- tt$brts
  model <- if (tt$kind == "cr") "cr" else "dd"
  lx    <- if (model == "dd") rf$lx else NULL
  dl    <- function(th) ll_exact(as.numeric(th), brts, model = model, lx = lx) - rf$mle_loglik

  # ONE box for all three arms.  The seed is the tree's, not the job's, so the
  # three configs compute the same box independently -- if each seeded from its
  # own job seed they would search inside different boxes and the comparison
  # would be confounded by the thing it is meant to hold fixed.
  set.seed(val_seed(paste0(job$tree_id, "-box")))
  ab <- emphasis::auto_bounds(brts, model = model, link = "linear",
                              num_threads = 1L, verbose = FALSE)
  lb <- ab$lower_bound; ub <- ab$upper_bound
  base <- list(lower_bound = lb, upper_bound = ub, num_threads = 1L,
               max_missing = 1e4, verbose = FALSE)

  t0 <- proc.time()[3]
  set.seed(job$seed)
  start <- switch(job$config,
    I2 = (lb + ub) / 2,
    I3 = {
      g <- emphasis::estimate_rates(brts, method = "gam", model = model,
             control = c(base, list(sample_size = 200L, grid_points = 12L)))
      as.numeric(g$pars)
    },
    I1 = {
      g <- emphasis::estimate_rates(brts, method = "gam", model = model,
             control = c(base, list(sample_size = 200L, grid_points = 12L)))
      c2 <- emphasis::estimate_rates(brts, method = "cem", model = model,
              init_pars = as.numeric(g$pars),
              control = c(base, list(num_particles = 50L, num_trees = 5L,
                                     max_iter = 20L)))
      as.numeric(c2$pars)
    },
    stop("unknown init config: ", job$config))
  t_init <- proc.time()[3] - t0
  start <- pmin(pmax(start, lb), ub)

  set.seed(job$seed + 1L)
  fit <- emphasis::estimate_rates(brts, method = "mcem", model = model,
           init_pars = start,
           control = c(base, list(sampling = "bdi", sample_size = 200L,
                                  max_iter = 400L,
                                  max_time = max(30, job$timeout_s - t_init - 30))))
  list(start_pars      = start,
       start_delta_ell = dl(start),
       pars            = as.numeric(fit$pars),
       delta_ell       = dl(fit$pars),
       ll_exact_mle    = rf$mle_loglik,
       stop_reason     = fit$stop_reason %||% NA_character_,
       iterations      = fit$iterations %||% NA_integer_,
       init_seconds    = t_init,
       box_contains_mle0 = all(rf$mle >= lb) && all(rf$mle <= ub),
       auto_lower = lb, auto_upper = ub)
}

# --- full pipeline ---------------------------------------------------------
.val_pipeline <- function(job, tt, rf) {
  brts <- tt$brts
  model <- if (tt$kind == "cr") "cr" else "dd"
  set.seed(job$seed)
  p <- emphasis::emphasis_pipeline(
    brts, model = model, link = "linear",
    stages = c("bounds", "gam", "cem", "mcem"),
    control = list(num_threads = 1L, rho = 1,
                   max_time = max(60, job$timeout_s - 60)),
    verbose = FALSE)
  lx <- if (model == "dd") rf$lx else NULL
  stage_pars <- lapply(p$fits, function(f) if (is.null(f)) NULL else as.numeric(f$pars))
  stage_ll <- vapply(stage_pars, function(pp)
    if (is.null(pp)) NA_real_ else ll_exact(pp, brts, model = model, lx = lx), 1)
  fin <- as.numeric(p$pars)
  # emphasis_pipeline returns the auto_bounds object in $bounds and the stage
  # table in $log (not $run_log)
  lb <- p$bounds$lower_bound %||% NA
  ub <- p$bounds$upper_bound %||% NA
  inbox <- function(th) if (any(is.na(lb)) || any(is.na(ub))) NA
                        else all(th >= lb) && all(th <= ub)
  list(pars = fin,
       ll_exact_hat = ll_exact(fin, brts, model = model, lx = lx),
       ll_exact_mle = rf$mle_loglik,
       delta_ell = ll_exact(fin, brts, model = model, lx = lx) - rf$mle_loglik,
       loglik = p$loglik %||% NA_real_, AIC = p$AIC %||% NA_real_,
       run_log = p$log, mcem_trace = p$mcem_trace,
       cond_used = p$cond, best_stage = p$best_stage %||% NA_character_,
       stage_pars = stage_pars, stage_ll_exact = stage_ll,
       auto_lower = lb, auto_upper = ub,
       box_contains_mle0 = inbox(rf$mle),
       box_contains_mle1 = if (model == "cr") inbox(rf$mle_cond1) else NA,
       mle_cond1 = if (model == "cr") rf$mle_cond1 else NA,
       mle = rf$mle, se = rf$se)
}

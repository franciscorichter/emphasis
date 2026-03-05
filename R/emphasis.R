.mcem_dynamic_fresh <- function(brts,
                      pars,
                      sample_size,
                      max_missing,
                      max_lambda,
                      lower_bound,
                      upper_bound,
                      xtol,
                      tol,
                      burnin,
                      num_threads,
                      return_trees = FALSE,
                      verbose = FALSE,
                      conditional = NULL,
                      model = c(0L, 0L, 0L),
                      link = 0L,
                      ...) {
  if (inherits(brts, "phylo")) {
    brts <- sort(ape::branching.times(brts), decreasing = TRUE)
  }
  if (!is.numeric(brts)) {
    stop("`brts` must be numeric or a `phylo` object.")
  }
  if (!is.numeric(pars)) {
    stop("`pars` must be numeric.")
  }
  if (!is.null(conditional) && !is.function(conditional)) {
    stop("`conditional` must be a function or NULL.")
  }

  copy_trees <- isTRUE(return_trees)
  maxN <- max(1L, as.integer(10 * sample_size))

  sde <- Inf
  i <- 0
  mcem <- NULL

  while (is.infinite(sde) || sde > tol) {
    i <- i + 1
    results <- tryCatch(
      em_cpp(brts = brts,
             init_pars = pars,
             sample_size = sample_size,
             maxN = maxN,
             max_missing = max_missing,
             max_lambda = max_lambda,
             lower_bound = lower_bound,
             upper_bound = upper_bound,
             xtol_rel = xtol,
             num_threads = num_threads,
             copy_trees = copy_trees,
             model = as.integer(model),
             link = as.integer(link),
             rconditional = conditional),
      error = function(e) {
        warning(".mcem_dynamic_fresh: E-step failed (", conditionMessage(e),
                "); stopping at iteration ", i, ".")
        NULL
      }
    )
    if (is.null(results)) break

    pars <- results$estimates
    par_df <- as.data.frame(as.list(stats::setNames(pars, paste0("par", seq_along(pars)))))
    step <- cbind(par_df, data.frame(
      fhat = results$fhat,
      sample_size = sample_size,
      time = results$time
    ))
    mcem <- rbind(mcem, step)

    if (verbose) {
      message(sprintf("Iteration %d: log-likelihood estimate = %f", i, results$fhat))
    }

    if (i > burnin) {
      mcem_est <- mcem[max(1, floor(nrow(mcem) / 2)):nrow(mcem), , drop = FALSE]
      mcem_est <- mcem_est[is.finite(mcem_est$fhat), , drop = FALSE]
      if (nrow(mcem_est) == 0) {
        warning("No finite log-likelihood estimates after burn-in; stopping .mcem_dynamic_fresh().")
        break
      }
      sde <- stats::sd(mcem_est$fhat) / sqrt(nrow(mcem_est))
      if (verbose) {
        message(sprintf("Iteration %d: log-likelihood SE = %f", i, sde))
      }
    } else if (verbose) {
      message(sprintf("Burn-in iteration %d/%d", i, burnin))
    }

    if (i >= 1000) {
      warning(".mcem_dynamic_fresh reached 1000 iterations without meeting tolerance.")
      break
    }
  }

  list(
    mcem = mcem,
    pars = pars,
    iterations = i,
    se = sde
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

mcEM_step <- function(brts,
                      pars,
                      sample_size,
                      soc,
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
                      ...) {
  if (inherits(brts, "phylo")) {
    brts <- ape::branching.times(brts)
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
    results <- em_cpp(brts = brts,
                      init_pars = pars,
                      sample_size = sample_size,
                      maxN = maxN,
                      soc = soc,
                      max_missing = max_missing,
                      max_lambda = max_lambda,
                      lower_bound = lower_bound,
                      upper_bound = upper_bound,
                      xtol_rel = xtol,
                      num_threads = num_threads,
                      copy_trees = copy_trees,
                      rconditional = conditional)

    pars <- results$estimates
    step <- data.frame(
      par1 = pars[1],
      par2 = pars[2],
      par3 = pars[3],
      par4 = pars[4],
      fhat = results$fhat,
      sample_size = sample_size,
      time = results$time
    )
    mcem <- rbind(mcem, step)

    if (verbose) {
      message(sprintf("Iteration %d: log-likelihood estimate = %f", i, results$fhat))
    }

    if (i > burnin) {
      mcem_est <- mcem[max(1, floor(nrow(mcem) / 2)):nrow(mcem), , drop = FALSE]
      mcem_est <- mcem_est[is.finite(mcem_est$fhat), , drop = FALSE]
      if (nrow(mcem_est) == 0) {
        warning("No finite log-likelihood estimates after burn-in; stopping mcEM_step().")
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
      warning("mcEM_step reached 1000 iterations without meeting tolerance.")
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

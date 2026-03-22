# --------------------------------------------------------------------------- #
#  Notation used throughout this file                                         #
#                                                                             #
#  z_i   : augmented tree simulated from proposal q(z | obs, theta_i)        #
#  logf(z, theta) = log p(obs, z | theta)  -- model log-likelihood           #
#  logg(z, theta) = log q(z | obs, theta)  -- proposal log-probability       #
#  lw(z_i, theta) = logf(z_i, theta) - logg(z_i, theta)                     #
#                   IS log-weight (both terms evaluated at same theta)        #
#  fhat(theta) = log_mean_exp over a set of lw values                        #
#              = IS estimate of log p(obs | theta)                            #
# --------------------------------------------------------------------------- #


# Null-safe integer coercion (for rejection counts that may be NULL)
.n0 <- function(x) if (is.null(x)) 0L else as.integer(x)

# Total rejection count across all failure types
.total_rejected <- function(raw) {
  .n0(raw$rejected) + .n0(raw$rejected_zero_weights) +
    .n0(raw$rejected_overruns) + .n0(raw$rejected_lambda)
}

# --------------------------------------------------------------------------- #
#  IS log-likelihood estimation                                               #
# --------------------------------------------------------------------------- #

#' IS log-likelihood estimate (with optional bias correction)
#'
#' Two estimators are available:
#'
#' \strong{Default} (\code{bias_correct = FALSE}): the standard log-sum-exp
#' estimator \eqn{\log\bar{L} = \log\!\left(\frac{1}{n}\sum_i L_i\right)},
#' computed in numerically stable form.  \eqn{\bar{L}} is an unbiased estimator
#' of \eqn{p(\mathrm{obs}\mid\theta)}, but \eqn{\log\bar{L}} slightly
#' underestimates \eqn{\log p} (Jensen's inequality).
#'
#' \strong{Bias-corrected} (\code{bias_correct = TRUE}): uses the exact identity
#' \deqn{\log\bar{L} = \hat\ell + \log\!\left(1 + \sum_{k=2}^{\infty}
#'   \frac{m_k}{k!}\right), \qquad
#'   m_k = \frac{1}{n}\sum_i(\log L_i - \hat\ell)^k}
#' truncated at order \code{K} (default 2).  The K=2 approximation is
#' \eqn{\hat\ell + \log(1 + m_2/2)}, equivalent to the log-normal correction
#' \eqn{\mu + \sigma^2/2}.  With the full series (\eqn{K\to\infty}) the two
#' estimators are identical; the truncation is the approximation.
#'
#' Trees that hit computational limits (\code{max_missing}, \code{max_lambda})
#' are discarded entirely -- they are simulator failures, not IS samples.
#' Trees that completed augmentation but have zero IS weight (\code{w = 0}
#' because \code{lambda = 0} under the current parameters) are legitimate IS
#' samples and are included in the denominator via \code{n_zero_weight}.
#'
#' @param logf Numeric vector of \code{log p(obs, z | theta)} for valid trees.
#' @param log_q Numeric vector of \code{log q(z | obs, theta)} for valid trees.
#' @param bias_correct Logical; use moment-based bias correction (default \code{FALSE}).
#' @param K Integer truncation order for the moment series (default \code{2}).
#' @param n_zero_weight Integer; number of completed augmentations that had
#'   zero IS weight.  These contribute \code{w = 0} to the IS sum and inflate
#'   the denominator accordingly (default \code{0}).
#' @keywords internal
.is_fhat <- function(logf, log_q, bias_correct = FALSE, K = 2L,
                     n_zero_weight = 0L) {
  lw      <- logf - log_q
  n_valid <- length(lw)
  if (n_valid == 0L) return(NA_real_)

  # Denominator includes zero-weight trees (completed augmentations where
  # w=0 because lambda=0 under the model).  Overrun/lambda rejections are
  # excluded: those are computational failures independent of theta.
  S_completed <- n_valid + as.integer(n_zero_weight)

  if (!bias_correct || n_valid == 1L) {
    # Standard log-mean-exp
    max_lw <- max(lw)
    if (!is.finite(max_lw)) return(NA_real_)
    return(log(sum(exp(lw - max_lw))) + max_lw - log(S_completed))
  }

  # Moment-based bias-corrected estimator truncated at order K.
  # When n_zero_weight > 0, the zero-weight trees contribute w=0 to the sum
  # but inflate the denominator, so we fall back to the standard estimator
  # which handles S_completed correctly.
  if (n_zero_weight > 0L) {
    max_lw <- max(lw)
    if (!is.finite(max_lw)) return(NA_real_)
    return(log(sum(exp(lw - max_lw))) + max_lw - log(S_completed))
  }

  ell_hat <- mean(lw)
  if (!is.finite(ell_hat)) return(NA_real_)

  d <- lw - ell_hat
  moment_sum <- 0
  for (k in seq(2L, K)) {
    moment_sum <- moment_sum + mean(d^k) / factorial(k)
  }

  inner <- 1 + moment_sum
  if (inner <= 0) return(NA_real_)

  ell_hat + log(inner)
}


#' Effective Sample Size from IS log-weights
#'
#' \eqn{\text{ESS} = (\sum w_i)^2 / \sum w_i^2}, computed in log-space for
#' numerical stability.  Measures IS efficiency: ESS close to \eqn{n} means
#' all trees contribute roughly equally; ESS close to 1 means one tree
#' dominates.
#' @keywords internal
.ess_from_lw <- function(lw) {
  lw <- lw[is.finite(lw)]
  if (length(lw) == 0L) return(NA_real_)
  m <- max(lw)
  w <- exp(lw - m)
  (sum(w))^2 / sum(w^2)
}


# --------------------------------------------------------------------------- #
#  Particle simulation                                                         #
# --------------------------------------------------------------------------- #

#' Draw augmented trees at a single particle (internal inference use only)
#'
#' Calls \code{\link{augment_trees}} to sample \code{sample_size} augmented
#' trees from the proposal \eqn{q(z \mid \text{obs}, \theta)}.
#'
#' Returns the raw \code{augment_trees} output: \code{logf}, \code{logg},
#' \code{trees}, and rejection counts for adaptive limit escalation.
#'
#' @keywords internal
.simulate_particle <- function(brts, pars8, model_bin, link,
                               sample_size, maxN, max_missing, max_lambda,
                               num_threads) {
  maxN_int <- if (is.null(maxN) || is.na(maxN))
    max(2000L, 200L * as.integer(sample_size)) else as.integer(maxN)
  tryCatch(
    augment_trees(
      brts        = as.numeric(brts),
      pars        = as.numeric(pars8),
      sample_size = as.integer(sample_size),
      maxN        = maxN_int,
      max_missing = as.integer(max_missing),
      max_lambda  = as.numeric(max_lambda),
      num_threads = as.integer(num_threads),
      model       = as.integer(model_bin),
      link        = as.integer(link)
    ),
    error = function(e) NULL
  )
}


# --------------------------------------------------------------------------- #
#  Population state                                                            #
# --------------------------------------------------------------------------- #

#' Initialise the particle population
#'
#' Returns a list with four parallel components of length \code{num_points}:
#' \describe{
#'   \item{\code{pars}}{Data frame of parameter vectors.}
#'   \item{\code{fhat}}{Numeric; \code{NA} = needs evaluation.}
#'   \item{\code{log_q}}{List; per-particle vector of \code{log q(z | obs, theta)}
#'     for cached trees. \code{NULL} = no trees cached.}
#'   \item{\code{trees}}{List; per-particle list of raw augmented-tree data
#'     frames (C++ format for \code{loglikelihood()}). \code{NULL} = not cached.}
#' }
#' @keywords internal
.init_population <- function(num_points, lower_bound, upper_bound) {
  pars <- .init_particle_grid(num_points, lower_bound, upper_bound)
  list(
    pars  = pars,
    fhat  = rep(NA_real_, num_points),
    log_q = vector("list", num_points),
    trees = vector("list", num_points)
  )
}

#' @keywords internal
.init_particle_grid <- function(num_points, lower_bound, upper_bound) {
  num_dim <- length(lower_bound)
  pars <- matrix(
    stats::runif(num_dim * num_points, min = lower_bound, max = upper_bound),
    nrow = num_points, ncol = num_dim, byrow = TRUE
  )
  df <- as.data.frame(pars)
  names(df) <- paste0("par", seq_len(num_dim))
  df
}


# --------------------------------------------------------------------------- #
#  Particle evaluation                                                         #
# --------------------------------------------------------------------------- #

#' Evaluate particles independently (Mode 1)
#'
#' Only evaluates particles with \code{NA} fhat.  For each such particle
#' \eqn{\theta_i}, simulates \code{sample_size} trees from
#' \eqn{q(z \mid \text{obs}, \theta_i)} and computes:
#' \deqn{\hat{f}(\theta_i) = \log\text{mean}\exp(\text{logf}(z, \theta_i) -
#'   \log q(z, \theta_i))}
#' Elite particles with cached \code{fhat} are skipped entirely.  Their
#' cached trees and \code{log_q} values are also stored for potential use
#' by the shared evaluator if the mode is later switched.
#'
#' @return Named list: updated \code{pop}, \code{rej_lambda} and
#'   \code{rej_overruns} vectors (length = number of particles evaluated).
#' @keywords internal
.eval_independent <- function(pop, input, num_threads) {
  to_eval  <- which(is.na(pop$fhat))
  if (length(to_eval) == 0L) {
    return(list(pop = pop, rej_lambda = 0L, rej_overruns = 0L,
                n_simulated = 0L))
  }

  # Adaptive: parallelize at R level only when there's enough work per
  # particle.  When R-parallel, each worker uses 1 C++ thread.
  # When sequential, pass all threads to C++ TBB.
  use_r_parallel <- (num_threads > 1L &&
                     .Platform$OS.type == "unix" &&
                     input$sample_size >= 10L &&
                     length(to_eval) >= num_threads)
  cpp_threads <- if (use_r_parallel) 1L else num_threads

  # Worker: evaluate one particle
  eval_one <- function(idx) {
    i   <- to_eval[idx]
    raw <- .simulate_particle(
      input$brts, as.numeric(pop$pars[i, ]), input$model, input$link,
      input$sample_size, input$maxN,
      input$max_missing, input$max_lambda, cpp_threads
    )
    if (!is.null(raw) && length(raw$logf) > 0L && !anyNA(raw$logf)) {
      fhat_i  <- .is_fhat(raw$logf, raw$logg,
                          bias_correct = isTRUE(input$bias_correct),
                          n_zero_weight = .n0(raw$rejected_zero_weights))
      log_q_i <- raw$logg
      trees_i <- raw$trees
    } else {
      fhat_i  <- NA_real_
      log_q_i <- NULL
      trees_i <- NULL
    }
    list(
      i        = i,
      fhat     = fhat_i,
      log_q    = log_q_i,
      trees    = trees_i,
      rej_lam  = if (!is.null(raw)) raw$rejected_lambda   else 0L,
      rej_over = if (!is.null(raw)) raw$rejected_overruns else 0L
    )
  }

  if (use_r_parallel) {
    results <- parallel::mclapply(
      seq_along(to_eval), eval_one,
      mc.cores = num_threads
    )
  } else {
    results <- lapply(seq_along(to_eval), eval_one)
  }

  # Collect results back into pop
  rej_lam  <- 0L
  rej_over <- 0L
  for (res in results) {
    if (!inherits(res, "try-error") && !is.null(res)) {
      i <- res$i
      if (!is.na(res$fhat)) {
        pop$fhat[i]    <- res$fhat
        pop$log_q[[i]] <- res$log_q
        pop$trees[[i]] <- res$trees
      }
      rej_lam  <- rej_lam  + res$rej_lam
      rej_over <- rej_over + res$rej_over
    }
  }

  list(pop          = pop,
       rej_lambda   = rej_lam,
       rej_overruns = rej_over,
       n_simulated  = length(to_eval))
}


#' Evaluate particles using the shared (pooled) tree estimator (Mode 2)
#'
#' Simulates fresh trees for every particle, pools them into a single set,
#' and cross-evaluates: for each particle \eqn{\theta_j}, computes
#' \deqn{\hat{f}(\theta_j) = \log\text{mean}\exp\bigl(
#'   \text{logf}(z_i, \theta_j) - \log q(z_i, \theta_j)\bigr)}
#' where both numerator and denominator are evaluated at the target particle
#' \eqn{\theta_j}.  The trees \eqn{z_i} may have been simulated at any
#' particle; what matters for IS correctness is that the proposal density
#' \eqn{q} is evaluated at the same \eqn{\theta_j} used for \eqn{p}.
#'
#' Tree caches are cleared between iterations (see \code{.resample_particles}
#' with \code{clear_cache = TRUE}), so the pool contains only the current
#' iteration's trees.
#'
#' @return Named list: updated \code{pop}, \code{rej_lambda},
#'   \code{rej_overruns}, \code{n_simulated}.
#' @keywords internal
.eval_shared <- function(pop, input, num_threads) {
  n <- nrow(pop$pars)

  # Step 1: simulate trees only for particles that do not have a cache yet
  needs_sim <- which(sapply(pop$trees, is.null))
  rej_lam   <- 0L
  rej_over  <- 0L

  for (i in needs_sim) {
    raw <- .simulate_particle(
      input$brts, as.numeric(pop$pars[i, ]), input$model, input$link,
      input$sample_size, input$maxN,
      input$max_missing, input$max_lambda, num_threads
    )
    if (!is.null(raw) && length(raw$logf) > 0L && !anyNA(raw$logf)) {
      pop$log_q[[i]] <- raw$logg
      pop$trees[[i]] <- raw$trees
    }
    rej_lam  <- rej_lam  + (if (!is.null(raw)) raw$rejected_lambda   else 0L)
    rej_over <- rej_over + (if (!is.null(raw)) raw$rejected_overruns else 0L)
  }

  # Step 2: pool all trees
  has_trees  <- !sapply(pop$trees, is.null)
  if (!any(has_trees)) {
    return(list(pop = pop, rej_lambda = rej_lam, rej_overruns = rej_over,
                n_simulated = length(needs_sim)))
  }
  all_trees <- do.call(c, pop$trees[has_trees])

  # Step 3: for each particle theta_j, evaluate logf AND logg at theta_j
  # over the full pool.  Both are evaluated at the TARGET particle, giving
  # correct IS weights: logf(z_i, theta_j) - logg(z_i, theta_j).
  for (j in seq_len(n)) {
    ev <- tryCatch(
      eval_logf(
        pars  = as.numeric(pop$pars[j, ]),
        trees = all_trees,
        model = as.integer(input$model),
        link  = as.integer(input$link)
      ),
      error = function(e) NULL
    )
    if (!is.null(ev)) {
      pop$fhat[j] <- .is_fhat(ev$logf, ev$logg, bias_correct = FALSE)
    }
  }

  list(pop          = pop,
       rej_lambda   = rej_lam,
       rej_overruns = rej_over,
       n_simulated  = length(needs_sim))
}


#' Dispatch particle evaluation to the appropriate mode
#' @keywords internal
.eval_particles <- function(pop, input, num_threads) {
  if (isTRUE(input$shared_trees)) {
    .eval_shared(pop, input, num_threads)
  } else {
    .eval_independent(pop, input, num_threads)
  }
}


# --------------------------------------------------------------------------- #
#  Particle resampling                                                         #
# --------------------------------------------------------------------------- #

#' @keywords internal
.perturb_value <- function(par, lower, upper, sd_val, max_tries = 100L) {
  if (upper == lower) return(upper)
  for (. in seq_len(max_tries)) {
    val <- stats::rnorm(1L, mean = par, sd = sd_val)
    if (val > lower && val < upper) return(val)
  }
  min(max(stats::rnorm(1L, mean = par, sd = sd_val), lower), upper)
}

#' Resample particle population
#'
#' Selects the top \code{disc_prop} fraction of valid particles (elites) and
#' fills the remaining slots with perturbed copies of elites.  All particles
#' (including elites) have their \code{fhat} reset to \code{NA} so they are
#' re-evaluated at the next iteration, preventing "lucky particle" lock-in.
#'
#' When \code{clear_cache = FALSE} (Mode 1), elite tree/log_q caches are
#' carried forward.  When \code{clear_cache = TRUE} (Mode 2), all caches are
#' discarded so every particle simulates fresh trees at the next iteration.
#' This prevents stale trees from unconverged early iterations from polluting
#' the shared pool.
#'
#' @param clear_cache Logical; if \code{TRUE}, discard all tree/log_q caches.
#' @return Updated population list.
#' @keywords internal
.resample_particles <- function(pop, num_points, disc_prop,
                                lower_bound, upper_bound, sd_vec,
                                clear_cache = FALSE) {
  valid    <- which(!is.na(pop$fhat))
  if (length(valid) == 0L) stop("No valid particles to resample from.")

  fhat_valid <- pop$fhat[valid]
  k          <- max(1L, ceiling(disc_prop * length(valid)))
  elite_ix   <- valid[order(fhat_valid, decreasing = TRUE)[seq_len(k)]]

  elite_pars  <- pop$pars[elite_ix,  , drop = FALSE]
  elite_fhat  <- rep(NA_real_, k)     # always re-evaluate
  if (clear_cache) {
    elite_log_q <- vector("list", k)
    elite_trees <- vector("list", k)
  } else {
    elite_log_q <- pop$log_q[elite_ix]
    elite_trees <- pop$trees[elite_ix]
  }

  n_add <- num_points - k
  if (n_add > 0L) {
    src_ix   <- sample(k, n_add, replace = TRUE)
    new_pars <- elite_pars[src_ix, , drop = FALSE]
    for (j in seq_along(sd_vec)) {
      new_pars[, j] <- vapply(
        new_pars[, j], .perturb_value,
        lower   = lower_bound[j],
        upper   = upper_bound[j],
        sd_val  = sd_vec[j],
        FUN.VALUE = numeric(1L)
      )
    }
    list(
      pars  = rbind(elite_pars, new_pars),
      fhat  = c(elite_fhat, rep(NA_real_, n_add)),
      log_q = c(elite_log_q, vector("list", n_add)),
      trees = c(elite_trees, vector("list", n_add))
    )
  } else {
    list(pars = elite_pars, fhat = elite_fhat,
         log_q = elite_log_q, trees = elite_trees)
  }
}


# --------------------------------------------------------------------------- #
#  Bootstrap variance of IS log-likelihood                                     #
# --------------------------------------------------------------------------- #

#' Bootstrap variance of the bias-corrected IS log-likelihood estimate
#'
#' Draws \code{B} bootstrap replicates (sampling with replacement from the IS
#' log-weights \code{lw = logf - log_q}) and computes
#' \eqn{\hat\ell_{\rm BC,K}} for each replicate.  Returns the sample variance
#' of the \code{B} estimates, which approximates
#' \eqn{\mathrm{Var}(\hat\ell_{\rm BC,K}(L))} due to Monte Carlo noise.
#'
#' This variance is used downstream for AIC uncertainty quantification:
#' \eqn{\mathrm{Var}(\Delta\mathrm{AIC}) = 4[\mathrm{Var}(\hat\ell_1) +
#' \mathrm{Var}(\hat\ell_2)]}, enabling a Gaussian test of equal AIC.
#'
#' @param logf Numeric vector of \code{log p(obs, z | theta)} for valid trees.
#' @param log_q Numeric vector of \code{log q(z | obs, theta)} for valid trees.
#' @param K Moment truncation order (default \code{2}).
#' @param B Number of bootstrap replicates (default \code{200}).
#' @return Numeric scalar: estimated variance of the log-likelihood estimator.
#'   \code{NA} if fewer than 2 valid log-weights are available.
#' @keywords internal
.bootstrap_fhat_var <- function(logf, log_q, K = 2L, B = 200L) {
  lw <- logf - log_q
  n  <- length(lw)
  if (n < 2L) return(NA_real_)

  boots <- vapply(seq_len(B), function(.) {
    idx <- sample(n, n, replace = TRUE)
    .is_fhat(logf[idx], log_q[idx], bias_correct = TRUE, K = K)
  }, numeric(1L))

  valid <- is.finite(boots)
  if (sum(valid) < 2L) return(NA_real_)
  stats::var(boots[valid])
}


# --------------------------------------------------------------------------- #
#  CEM optimiser                                                               #
# --------------------------------------------------------------------------- #

#' Monte Carlo Cross-Entropy Method for diversification models
#'
#' Implements the Cross-Entropy Method (CEM) with an IS log-likelihood
#' objective.  A population of parameter particles is evaluated via IS
#' (\code{fhat}), the elite fraction is retained, and the remainder is
#' refilled by perturbing elite copies.
#'
#' \strong{Adaptive SD annealing:} The perturbation SD is decayed
#' multiplicatively (factor \code{sd_decay}, default 0.85) only when the
#' best log-likelihood has not improved by more than \code{tol} in the
#' current iteration.  When the likelihood is still improving, the SD is
#' held constant to maintain exploration.  A minimum SD floor
#' (\code{sd_min_frac * initial_sd}, default 1\% of initial) prevents
#' complete collapse.  This replaces the former linear annealing schedule
#' which coupled the SD decay to \code{max_iter} and could cause
#' premature collapse on easy problems or insufficient exploration on
#' hard ones.
#'
#' \strong{Stopping rules} (first triggered wins):
#' \enumerate{
#'   \item Plateau: best \code{fhat} has not improved by more than \code{tol}
#'     for \code{patience} consecutive iterations \emph{and} the SD has
#'     reached its floor -- the population has both converged in likelihood
#'     and exhausted its perturbation budget.
#'   \item Hard cap: \code{max_iter} iterations reached (safety limit).
#' }
#'
#' @param brts Numeric vector of branching times.
#' @param max_iter Maximum number of iterations (hard cap).
#' @param num_points Population size.
#' @param max_missing Max extinct lineages per augmentation.
#' @param sd_vec Initial perturbation SD per parameter (8-element).
#' @param lower_bound,upper_bound Search bounds (8-element).
#' @param maxN Max augmentation attempts per particle. Default \code{10}.
#' @param sample_size Trees simulated per particle per evaluation. Default \code{1}.
#' @param shared_trees If \code{TRUE}, use pooled-tree estimator (Mode 2).
#'   Default \code{FALSE} (independent, Mode 1).
#' @param bias_correct If \code{TRUE}, apply moment-based bias correction to
#'   the IS log-likelihood. Default \code{FALSE}.
#' @param disc_prop Elite fraction retained per iteration. Default \code{0.5}.
#' @param tol Minimum improvement in best \code{fhat} required to avoid
#'   incrementing the plateau counter.  Default \code{1e-4}.
#' @param patience Consecutive plateau iterations before early stopping.
#'   Default \code{5}.
#' @param sd_decay Multiplicative factor applied to \code{sd_vec} when the
#'   best log-likelihood has not improved.  Default \code{0.85}.  Set to
#'   \code{NULL} to use the legacy linear annealing schedule.
#' @param sd_min_frac Minimum SD as a fraction of the initial \code{sd_vec}.
#'   Default \code{0.01} (1\%).  SD will not decay below this floor.
#' @param verbose Print iteration summaries. Default \code{FALSE}.
#' @param num_threads Parallel threads. Default \code{1}.
#' @param model Length-3 binary integer vector \code{c(use_N, use_P, use_E)}.
#' @param link Integer link code: \code{0} = linear, \code{1} = exponential.
#' @return A list:
#'   \describe{
#'     \item{\code{best_loglik}}{fhat of best particle per completed iteration
#'       (length = actual iterations run, \eqn{\leq} \code{max_iter}).}
#'     \item{\code{best_pars}}{Parameter matrix, one row per iteration.}
#'     \item{\code{obtained_estim}}{Globally best parameter vector.}
#'     \item{\code{loglik_var}}{Bootstrap variance of \code{fhat} at best
#'       particle (\code{NA} if \code{sample_size = 1}).}
#'     \item{\code{converged}}{Character string naming the stopping rule that
#'       fired: \code{"annealing"}, \code{"plateau"}, or \code{"max_iter"}.}
#'     \item{\code{time}}{Elapsed \code{proc.time()}.}
#'   }
#' @keywords internal
#' @rawNamespace useDynLib(emphasis)
#' @rawNamespace import(nloptr)
#' @rawNamespace import(Rcpp)
#' @rawNamespace importFrom(RcppParallel, RcppParallelLibs)
emphasis_cem <- function(brts,
                         max_iter,
                         num_points,
                         max_missing,
                         sd_vec,
                         lower_bound,
                         upper_bound,
                         maxN         = 10L,
                         sample_size  = 1L,
                         shared_trees = FALSE,
                         bias_correct = FALSE,
                         disc_prop    = 0.5,
                         tol          = 1e-4,
                         patience     = 5L,
                         sd_decay     = 0.85,
                         sd_min_frac  = 0.01,
                         verbose      = FALSE,
                         num_threads  = 1L,
                         model        = c(0L, 0L, 0L),
                         link         = 0L,
                         cond_fun     = NULL,
                         max_time     = NULL) {

  if (!is.numeric(lower_bound) || !is.numeric(upper_bound))
    stop("lower_bound and upper_bound must be numeric vectors.")
  if (length(upper_bound) != length(lower_bound))
    stop("lower_bound and upper_bound must have the same length.")
  if (length(upper_bound) != length(sd_vec))
    stop("sd_vec must have the same length as lower_bound.")

  init_time <- proc.time()

  # Annealing schedule: adaptive (default) or legacy linear
  use_adaptive <- !is.null(sd_decay)
  if (use_adaptive) {
    sd_floor <- sd_vec * sd_min_frac   # minimum SD per parameter
  } else {
    alpha <- sd_vec / max_iter          # legacy linear annealing step
  }

  input <- list(
    brts         = brts,
    sample_size  = sample_size,
    maxN         = maxN,
    max_missing  = max_missing,
    max_lambda   = 1e6,
    lower_bound  = lower_bound,
    upper_bound  = upper_bound,
    shared_trees = shared_trees,
    bias_correct = bias_correct,
    model        = model,
    link         = link
  )

  n_pars      <- length(lower_bound)
  pop         <- .init_population(num_points, lower_bound, upper_bound)
  best_loglik <- numeric(max_iter)
  best_pars   <- matrix(NA_real_, nrow = max_iter, ncol = n_pars)

  # Per-iteration history (all particles, not just best)
  hist_fhat_all   <- vector("list", max_iter)
  hist_n_valid    <- integer(max_iter)
  hist_rej_lam    <- integer(max_iter)
  hist_rej_over   <- integer(max_iter)

  plateau_count <- 0L
  stop_reason   <- "max_iter"
  k_ran         <- 0L          # iterations actually completed
  final_pop     <- NULL        # population state at last completed iteration

  for (k in seq_len(max_iter)) {

    # --- Stop rule 1: annealing exhausted (legacy only) ---------------------
    if (!use_adaptive && all(sd_vec <= 0)) {
      stop_reason <- "annealing"
      break
    }

    result <- .eval_particles(pop, input, num_threads)
    pop    <- result$pop

    # Apply conditioning correction: fhat_cond = fhat - log(P_tree)
    if (!is.null(cond_fun)) {
      for (i in which(!is.na(pop$fhat))) {
        pop$fhat[i] <- pop$fhat[i] - cond_fun(as.numeric(pop$pars[i, ]))
      }
    }

    # Rescue: if ALL particles still have NA fhat, escalate limits
    all_failed <- FALSE
    while (all(is.na(pop$fhat))) {
      input$maxN        <- input$maxN        * 10L
      input$max_missing <- input$max_missing * 10L
      if (input$maxN > 10000L) {
        warning("emphasis_cem: all particles failed even after escalating limits. ",
                "Try increasing num_particles or widening the parameter bounds.")
        all_failed <- TRUE
        break
      }
      result <- .eval_particles(pop, input, num_threads)
      pop    <- result$pop
    }
    if (all_failed) {
      stop_reason <- "all_failed"
      break
    }

    # Adaptive limit escalation
    if (result$rej_overruns > 0L) input$max_missing <- input$max_missing * 1.1

    valid      <- !is.na(pop$fhat)
    fhat_valid <- pop$fhat[valid]
    pars_valid <- pop$pars[valid, , drop = FALSE]

    best_ix        <- which.max(fhat_valid)
    best_loglik[k] <- fhat_valid[best_ix]
    best_pars[k, ] <- as.numeric(pars_valid[best_ix, ])

    # Record per-iteration history
    hist_fhat_all[[k]] <- pop$fhat                 # all fhat (NA for invalid)
    hist_n_valid[k]    <- sum(valid)
    hist_rej_lam[k]    <- result$rej_lambda
    hist_rej_over[k]   <- result$rej_overruns

    k_ran     <- k
    final_pop <- pop   # snapshot before resample -- has fhat, log_q, trees

    # --- Stop rule 2: plateau -----------------------------------------------
    improving <- FALSE
    if (k > 1L && is.finite(best_loglik[k]) && is.finite(best_loglik[k - 1L])) {
      improvement <- best_loglik[k] - best_loglik[k - 1L]
      improving <- (improvement >= tol)
      plateau_count <- if (improving) 0L else plateau_count + 1L
    }
    if (use_adaptive) {
      # Adaptive: require plateau AND SD at floor
      sd_at_floor <- all(sd_vec <= sd_floor * (1 + 1e-8))
      if (plateau_count >= patience && sd_at_floor) {
        stop_reason <- "plateau"
        break
      }
    } else {
      # Legacy: plateau alone suffices
      if (plateau_count >= patience) {
        stop_reason <- "plateau"
        break
      }
    }

    # Time budget check
    if (!is.null(max_time)) {
      elapsed <- (proc.time() - init_time)[3]
      if (elapsed > max_time) {
        stop_reason <- "time_budget"
        if (verbose) cat(sprintf("Time budget reached (%.0fs > %ds)\n", elapsed, as.integer(max_time)))
        break
      }
    }

    pop    <- .resample_particles(pop, num_points, disc_prop,
                                  lower_bound, upper_bound, sd_vec,
                                  clear_cache = isTRUE(shared_trees))

    # SD annealing
    if (use_adaptive) {
      # Adaptive: decay only when not improving, respect floor
      if (!improving) {
        sd_vec <- pmax(sd_vec * sd_decay, sd_floor)
      }
    } else {
      # Legacy: linear decay every iteration
      sd_vec <- sd_vec - alpha
    }

    if (verbose) {
      lbl <- if (sample_size <= 1L) "best_lw" else "loglik*"
      sd_frac <- if (use_adaptive) mean(sd_vec / (sd_floor / sd_min_frac)) else
                 mean(sd_vec / (sd_vec + alpha))
      cat(sprintf(
        "Iter %3d/%d  %s=%8.4f  valid=%d/%d  plateau=%d/%d  sd=%.1f%%\n",
        k, max_iter, lbl, best_loglik[k],
        sum(valid), num_points,
        plateau_count, patience,
        sd_frac * 100))
    }
  }

  # Trim pre-allocated storage to actual iterations run
  best_loglik <- best_loglik[seq_len(k_ran)]
  best_pars   <- best_pars[seq_len(k_ran), , drop = FALSE]
  history <- list(
    fhat_all     = hist_fhat_all[seq_len(k_ran)],
    n_valid      = hist_n_valid[seq_len(k_ran)],
    rej_lambda   = hist_rej_lam[seq_len(k_ran)],
    rej_overruns = hist_rej_over[seq_len(k_ran)]
  )

  # Early exit if no iterations completed successfully
  if (k_ran == 0L) {
    if (verbose)
      cat(sprintf("Stopped: %s -- no iterations completed.\n", stop_reason))
    return(list(
      best_loglik    = numeric(0),
      best_pars      = best_pars[integer(0), , drop = FALSE],
      obtained_estim = rep(NA_real_, n_pars),
      loglik_var     = NA_real_,
      converged      = stop_reason,
      history        = list(fhat_all = list(), n_valid = integer(0),
                            rej_lambda = integer(0), rej_overruns = integer(0)),
      final_pop      = NULL,
      best_IS        = NULL,
      time           = proc.time() - init_time
    ))
  }

  # Use the final iteration's best particle, not the global max.
  # The global max over all iterations selects lucky single-tree spikes rather
  # than the converged estimate.
  best_pars_v <- best_pars[k_ran, ]

  # Final evaluation of best particle using num_trees draws.
  # Bootstrap variance is computed automatically from these draws whenever
  # num_trees > 1 (B = 200 bootstrap replicates, bias-corrected estimator).
  raw_best <- .simulate_particle(
    input$brts, as.numeric(best_pars_v), input$model, input$link,
    sample_size = as.integer(sample_size),
    maxN        = input$maxN,
    max_missing = input$max_missing,
    max_lambda  = input$max_lambda,
    num_threads = num_threads
  )

  best_IS    <- NULL
  loglik_var <- NA_real_
  if (!is.null(raw_best) && length(raw_best$logf) > 0L) {
    lw_all <- raw_best$logf - raw_best$logg
    fhat_best <- .is_fhat(raw_best$logf, raw_best$logg,
                           n_zero_weight = .n0(raw_best$rejected_zero_weights))
    # Apply conditioning correction
    if (!is.null(cond_fun) && is.finite(fhat_best)) {
      fhat_best <- fhat_best - cond_fun(as.numeric(best_pars_v))
    }
    best_IS <- list(
      logf       = raw_best$logf,
      log_q      = raw_best$logg,
      lw         = lw_all,
      trees      = raw_best$trees,
      fhat       = fhat_best,
      ESS        = .ess_from_lw(lw_all),
      n_trees    = as.integer(sample_size),
      n_rejected = .total_rejected(raw_best)
    )
    lw_fin <- lw_all[is.finite(lw_all)]
    if (length(lw_fin) >= 2L) {
      loglik_var <- .bootstrap_fhat_var(raw_best$logf, raw_best$logg,
                                        K = 2L, B = 200L)
    }
  }

  if (verbose)
    cat(sprintf("Stopped: %s after %d iterations.\n", stop_reason, k_ran))

  list(
    best_loglik    = best_loglik,
    best_pars      = best_pars,
    obtained_estim = best_pars_v,
    loglik_var     = loglik_var,
    converged      = stop_reason,
    history        = history,
    final_pop      = final_pop,
    best_IS        = best_IS,
    time           = proc.time() - init_time
  )
}

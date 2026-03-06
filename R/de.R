# --------------------------------------------------------------------------- #
#  Notation used throughout this file                                         #
#                                                                             #
#  z_i   : augmented tree simulated from proposal q(z | obs, theta_i)        #
#  logf(z, theta) = log p(obs, z | theta)  -- model log-likelihood           #
#  log_q(z, theta) = log q(z | obs, theta) -- proposal log-probability       #
#  lw(z_i, theta_j) = logf(z_i, theta_j) - log_q(z_i, theta_i)             #
#                     IS log-weight (denominator: simulation particle theta_i)#
#  fhat(theta_j) = log_mean_exp over a set of lw values                      #
#                = IS estimate of log p(obs | theta_j)                        #
# --------------------------------------------------------------------------- #


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
#' Rejected trees (\eqn{L_i = 0}) are excluded from the moment computation but
#' enter the denominator via \eqn{\log(n_{\rm valid}/n_{\rm total})}.
#'
#' @param bias_correct Logical; use moment-based bias correction (default \code{FALSE}).
#' @param K Integer truncation order for the moment series (default \code{2}).
#' @keywords internal
.is_fhat <- function(logf, log_q, n_rejected = 0L, bias_correct = FALSE, K = 2L) {
  lw      <- logf - log_q
  n_valid <- length(lw)
  if (n_valid == 0L) return(NA_real_)
  n_total <- n_valid + n_rejected

  if (!bias_correct || n_valid == 1L) {
    # Standard log-mean-exp (no useful correction with a single sample)
    max_lw <- max(lw)
    if (!is.finite(max_lw)) return(NA_real_)
    return(log(sum(exp(lw - max_lw))) + max_lw - log(n_total))
  }

  # Moment-based bias-corrected estimator truncated at order K
  # ell_hat = (1/n_valid) * sum(log L_i)   [mean of IS log-weights, valid only]
  ell_hat <- mean(lw)
  if (!is.finite(ell_hat)) return(NA_real_)

  # Central moments m_k = (1/n_valid) * sum((lw - ell_hat)^k)
  d <- lw - ell_hat
  moment_sum <- 0
  for (k in seq(2L, K)) {
    moment_sum <- moment_sum + mean(d^k) / factorial(k)
  }

  inner <- 1 + moment_sum
  if (inner <= 0) return(NA_real_)

  # Rejection adjustment: n_rejected trees have L_i = 0 → reduce by n_valid/n_total
  ell_hat + log(inner) + log(n_valid / n_total)
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
#' Calls \code{.augment_tree_internal} to sample \code{sample_size} augmented
#' trees from the proposal \code{q(z | obs, theta)}.
#'
#' Returns the raw \code{mc_loglik} output, which contains:
#' \code{logf} (log p(obs, z | theta) per tree), \code{logg} (log_q per tree),
#' \code{trees} (raw data frames for \code{loglikelihood()}), and rejection
#' counts for adaptive limit escalation.
#'
#' @keywords internal
.simulate_particle <- function(brts, pars8, model_bin, link,
                               sample_size, maxN, max_missing, max_lambda,
                               num_threads) {
  pars8_full <- as.numeric(.expand_pars(.contract_pars(pars8, model_bin), model_bin))
  maxN_int   <- if (is.null(maxN) || is.na(maxN))
    max(2000L, 200L * as.integer(sample_size)) else as.integer(maxN)
  tryCatch(
    mc_loglik(
      brts        = as.numeric(brts),
      pars        = pars8_full,
      sample_size = as.integer(sample_size),
      maxN        = maxN_int,
      max_missing = as.integer(max_missing),
      max_lambda  = as.numeric(max_lambda),
      lower_bound = rep(-1e6, 8L),
      upper_bound = rep( 1e6, 8L),
      xtol_rel    = 1e-3,
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
  rej_lam  <- rep(0L, length(to_eval))
  rej_over <- rep(0L, length(to_eval))

  for (idx in seq_along(to_eval)) {
    i   <- to_eval[idx]
    raw <- .simulate_particle(
      input$brts, as.numeric(pop$pars[i, ]), input$model, input$link,
      input$sample_size, input$maxN,
      input$max_missing, input$max_lambda, num_threads
    )
    if (!is.null(raw) && length(raw$logf) > 0L && !anyNA(raw$logf)) {
      # fhat(theta_i) = log_mean_exp(logf(z, theta_i) - log_q(z, theta_i))
      pop$fhat[i]    <- .is_fhat(raw$logf, raw$logg,
                                 n_rejected   = raw$rejected,
                                 bias_correct = isTRUE(input$bias_correct))
      pop$log_q[[i]] <- raw$logg    # log q(z | obs, theta_i) per tree
      pop$trees[[i]] <- raw$trees   # raw data frames for cross-evaluation
    }
    rej_lam[idx]  <- if (!is.null(raw)) raw$rejected_lambda   else 0L
    rej_over[idx] <- if (!is.null(raw)) raw$rejected_overruns else 0L
  }

  list(pop          = pop,
       rej_lambda   = rej_lam,
       rej_overruns = rej_over,
       n_simulated  = length(to_eval))
}


#' Evaluate particles using the shared (pooled) tree estimator (Mode 2)
#'
#' For particles without cached trees, simulates new trees first.  Then pools
#' all available trees across the entire population and, for each particle
#' \eqn{\theta_j}, computes:
#' \deqn{\hat{f}(\theta_j) = \log\text{mean}\exp\bigl(
#'   \text{logf}(z_i, \theta_j) - \log q(z_i, \theta_i)\bigr)}
#' where the sum runs over all pooled trees \eqn{z_i} regardless of which
#' particle simulated them.
#'
#' Unlike Mode 1, elite particles' \code{fhat} is \emph{always} recomputed
#' because the pool grows with each iteration (new trees enter the pool).
#' The savings versus Mode 1 is that elite trees are not re-simulated.
#'
#' \strong{Bias note:} this is a naive pooled IS estimator.  When proposals
#' differ across particles (\eqn{q(z \mid \theta_i) \neq q(z \mid \theta_j)}),
#' the estimator is biased.  Set \code{bias_correct = TRUE} to apply the
#' Multiple IS (MIS) correction once it is implemented.
#'
#' \strong{C++ limitation:} \code{loglikelihood()} currently hardcodes
#' \code{model = c(0,0,0)}; shared mode is therefore only correct for the
#' CR model until that is fixed.
#'
#' @param bias_correct Logical; if \code{TRUE}, apply MIS bias correction
#'   (not yet implemented — placeholder).
#' @return Named list: updated \code{pop}, \code{rej_lambda},
#'   \code{rej_overruns}, \code{n_simulated}.
#' @keywords internal
.eval_shared <- function(pop, input, num_threads, bias_correct = FALSE) {
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

  # Step 2: pool all trees and their log_q values
  has_trees  <- !sapply(pop$trees, is.null)
  if (!any(has_trees)) {
    return(list(pop = pop, rej_lambda = rej_lam, rej_overruns = rej_over,
                n_simulated = length(needs_sim)))
  }
  all_trees <- do.call(c, pop$trees[has_trees])
  all_log_q <- unlist(pop$log_q[has_trees])   # log q(z_i | obs, theta_i)

  # Step 3: for each particle theta_j, compute logf(z_i, theta_j) over the
  # full pool and aggregate into fhat(theta_j).
  # fhat(theta_j) = log_mean_exp(logf(z_i, theta_j) - log_q(z_i, theta_i))
  #
  # With bias_correct = TRUE: replace log_q(z_i, theta_i) with the MIS
  # denominator log((1/N) * sum_l q(z_i | obs, theta_l)) — TODO.
  for (j in seq_len(n)) {
    ev <- tryCatch(
      loglikelihood(
        pars         = as.numeric(pop$pars[j, ]),
        trees        = all_trees,
        logg         = all_log_q,
        plugin       = "",
        num_rejected = 0L
      ),
      error = function(e) NULL
    )
    if (!is.null(ev)) {
      if (bias_correct && !is.null(ev$logf)) {
        # Apply moment-based bias correction using the per-tree logf values
        # returned by loglikelihood() alongside all_log_q as IS denominators.
        # lw_j[i] = logf(z_i, theta_j) - log_q(z_i, theta_i) for all pooled trees
        pop$fhat[j] <- .is_fhat(ev$logf, all_log_q, n_rejected = 0L,
                                bias_correct = TRUE)
      } else {
        # Standard log-mean-exp (ev$fhat already computed by loglikelihood())
        pop$fhat[j] <- ev$fhat
      }
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
    .eval_shared(pop, input, num_threads,
                 bias_correct = isTRUE(input$bias_correct))
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

#' Resample particle population, preserving elite tree cache
#'
#' Selects the top \code{disc_prop} fraction of valid particles (elites) and
#' carries their tree/log_q caches forward.  All particles (including elites)
#' have their \code{fhat} reset to \code{NA} so they are re-evaluated with a
#' freshly simulated tree at the next iteration.  This prevents "lucky
#' particle" lock-in — a particle whose previous high \code{fhat} was driven
#' by a single favourable tree draw rather than by being a genuinely good
#' parameter point.
#'
#' \strong{Mode 1} benefits from elites concentrating the proposal distribution
#' (new particles are perturbed copies of elites); the re-evaluation cost is
#' the price of unbiased position estimates.
#'
#' \strong{Mode 2:} \code{.eval_shared} always recomputes \code{fhat} over the
#' full pooled tree set anyway; this reset has no additional effect there.
#'
#' @return Updated population list.
#' @keywords internal
.resample_particles <- function(pop, num_points, disc_prop,
                                lower_bound, upper_bound, sd_vec) {
  valid    <- which(!is.na(pop$fhat))
  if (length(valid) == 0L) stop("No valid particles to resample from.")

  fhat_valid <- pop$fhat[valid]
  k          <- max(1L, ceiling(disc_prop * length(valid)))
  elite_ix   <- valid[order(fhat_valid, decreasing = TRUE)[seq_len(k)]]

  elite_pars  <- pop$pars[elite_ix,  , drop = FALSE]
  elite_fhat  <- rep(NA_real_, k)     # always re-evaluate; prevents lucky-particle lock-in
  elite_log_q <- pop$log_q[elite_ix]
  elite_trees <- pop$trees[elite_ix]

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
#' @param logf Numeric vector of \code{log p(obs, z | theta)} values.
#' @param log_q Numeric vector of \code{log q(z | obs, theta)} values.
#' @param n_rejected Integer; number of rejected (zero-weight) trees.
#' @param K Moment truncation order (default \code{2}).
#' @param B Number of bootstrap replicates (default \code{200}).
#' @return Numeric scalar: estimated variance of the log-likelihood estimator.
#'   \code{NA} if fewer than 2 valid log-weights are available.
#' @keywords internal
.bootstrap_fhat_var <- function(logf, log_q, n_rejected = 0L, K = 2L, B = 200L) {
  lw <- logf - log_q
  n  <- length(lw)
  if (n < 2L) return(NA_real_)

  boots <- vapply(seq_len(B), function(.) {
    idx <- sample(n, n, replace = TRUE)
    .is_fhat(logf[idx], log_q[idx], n_rejected = 0L, bias_correct = TRUE, K = K)
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
#' refilled by perturbing elite copies.  The perturbation SD is annealed
#' linearly from \code{sd_vec} to zero over at most \code{max_iter}
#' iterations.
#'
#' \strong{Stopping rules} (first triggered wins):
#' \enumerate{
#'   \item Annealing exhausted: \code{all(sd_vec <= 0)} — the population has
#'     fully collapsed; further iterations produce identical particles.
#'   \item Plateau: best \code{fhat} has not improved by more than \code{tol}
#'     for \code{patience} consecutive iterations.
#'   \item Hard cap: \code{max_iter} iterations reached.
#' }
#'
#' @param brts Numeric vector of branching times.
#' @param max_iter Maximum number of iterations (hard cap).
#' @param num_points Population size.
#' @param max_missing Max extinct lineages per augmentation.
#' @param sd_vec Initial perturbation SD per parameter (8-element).
#' @param lower_bound,upper_bound Search bounds (8-element).
#' @param maxN Max augmentation attempts per particle. Default \code{10}.
#' @param max_lambda Max speciation rate. Default \code{500}.
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
#' @param verbose Print iteration summaries. Default \code{FALSE}.
#' @param num_threads Parallel threads. Default \code{1}.
#' @param model Length-3 binary integer vector \code{c(use_N, use_P, use_E)}.
#' @param link Integer link code: \code{0} = linear, \code{1} = exponential.
#' @param n_boot Bootstrap replicates for MC variance at best particle.
#'   \code{0} (default) skips bootstrap; \code{loglik_var} will be \code{NA}.
#' @return A list:
#'   \describe{
#'     \item{\code{best_loglik}}{fhat of best particle per completed iteration
#'       (length = actual iterations run, \eqn{\leq} \code{max_iter}).}
#'     \item{\code{best_pars}}{Parameter matrix, one row per iteration.}
#'     \item{\code{obtained_estim}}{Globally best parameter vector.}
#'     \item{\code{loglik_var}}{Bootstrap variance of \code{fhat} at best
#'       particle (\code{NA} if \code{n_boot = 0}).}
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
                         max_lambda   = 500,
                         sample_size  = 1L,
                         shared_trees = FALSE,
                         bias_correct = FALSE,
                         disc_prop    = 0.5,
                         tol          = 1e-4,
                         patience     = 5L,
                         verbose      = FALSE,
                         num_threads  = 1L,
                         model        = c(0L, 0L, 0L),
                         link         = 0L,
                         n_boot       = 0L) {

  if (!is.numeric(lower_bound) || !is.numeric(upper_bound))
    stop("lower_bound and upper_bound must be numeric vectors.")
  if (length(upper_bound) != length(lower_bound))
    stop("lower_bound and upper_bound must have the same length.")
  if (length(upper_bound) != length(sd_vec))
    stop("sd_vec must have the same length as lower_bound.")

  init_time <- proc.time()
  alpha <- sd_vec / max_iter   # linear annealing step

  input <- list(
    brts         = brts,
    sample_size  = sample_size,
    maxN         = maxN,
    max_missing  = max_missing,
    max_lambda   = max_lambda,
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

    # --- Stop rule 1: annealing exhausted -----------------------------------
    if (all(sd_vec <= 0)) {
      stop_reason <- "annealing"
      break
    }

    result <- .eval_particles(pop, input, num_threads)
    pop    <- result$pop

    # Rescue: if ALL particles still have NA fhat, escalate limits
    while (all(is.na(pop$fhat))) {
      input$maxN        <- input$maxN        * 10L
      input$max_missing <- input$max_missing * 10L
      input$max_lambda  <- input$max_lambda  * 10
      if (input$maxN > 10000L)
        stop("emphasis_cem: all particles failed; maxN limit exceeded.")
      result <- .eval_particles(pop, input, num_threads)
      pop    <- result$pop
    }

    # Adaptive limit escalation
    if (result$rej_lambda   > 0L) input$max_lambda  <- input$max_lambda  * 1.1
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
    final_pop <- pop   # snapshot before resample — has fhat, log_q, trees

    # --- Stop rule 2: plateau -----------------------------------------------
    if (k > 1L && is.finite(best_loglik[k]) && is.finite(best_loglik[k - 1L])) {
      improvement <- best_loglik[k] - best_loglik[k - 1L]
      plateau_count <- if (improvement < tol) plateau_count + 1L else 0L
    }
    if (plateau_count >= patience) {
      stop_reason <- "plateau"
      break
    }

    pop    <- .resample_particles(pop, num_points, disc_prop,
                                  lower_bound, upper_bound, sd_vec)
    sd_vec <- sd_vec - alpha

    if (verbose)
      cat(sprintf("Iter %d/%d  best fhat: %.4f  plateau: %d/%d  sim: %d/%d\n",
                  k, max_iter, best_loglik[k],
                  plateau_count, patience,
                  result$n_simulated, num_points))
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

  global_best <- which.max(best_loglik)
  best_pars_v <- best_pars[global_best, ]

  # Final evaluation of best particle: always simulate to get IS data
  # (drives both best_IS diagnostics and bootstrap variance)
  n_final  <- max(as.integer(sample_size), as.integer(n_boot), 1L)
  raw_best <- .simulate_particle(
    input$brts, as.numeric(best_pars_v), input$model, input$link,
    sample_size = n_final,
    maxN        = input$maxN,
    max_missing = input$max_missing,
    max_lambda  = input$max_lambda,
    num_threads = num_threads
  )

  best_IS    <- NULL
  loglik_var <- NA_real_
  if (!is.null(raw_best) && length(raw_best$logf) > 0L) {
    lw_all <- raw_best$logf - raw_best$logg
    lw_fin <- lw_all[is.finite(lw_all)]
    best_IS <- list(
      logf       = raw_best$logf,
      log_q      = raw_best$logg,
      lw         = lw_all,
      trees      = raw_best$trees,
      fhat       = .is_fhat(raw_best$logf, raw_best$logg,
                            n_rejected = raw_best$rejected),
      ESS        = .ess_from_lw(lw_all),
      n_trees    = n_final,
      n_rejected = raw_best$rejected
    )
    if (n_boot > 0L && length(lw_fin) >= 2L) {
      loglik_var <- .bootstrap_fhat_var(raw_best$logf, raw_best$logg,
                                        n_rejected = raw_best$rejected,
                                        K = 2L, B = as.integer(n_boot))
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

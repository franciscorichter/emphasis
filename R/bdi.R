# ═══════════════════════════════════════════════════════════════
# BDI exact augmentation sampler
#
# Draws augmented trees from the exact conditional distribution
# p(A | T, theta) via an aggregate Gillespie algorithm.
# Under constant rates: ESS = S (zero-variance proposal).
# Under DD/endogenous: self-consistent backward-forward iteration.
# ═══════════════════════════════════════════════════════════════


# --------------------------------------------------------------------------- #
#  Survival probability helpers                                                #
# --------------------------------------------------------------------------- #

# CR analytical solution
# p(t) = (lam0-mu0) / (lam0 - mu0*exp(-(lam0-mu0)*(tp-t))), valid for any
# lam0, mu0 >= 0. At lam0 == mu0 the general form is 0/0; its limit is
# 1/(1 + lam0*(tp-t)), used when |lam0-mu0| is below a relative tolerance.
#
# Written with s = |lam0-mu0| so that the exponential argument is never
# positive: exp(+s*(tp-t)) overflows for mu0 > lam0 once s*(tp-t) passes
# ~709, which is reachable through a user-supplied box on the exponential
# link.  For mu0 > lam0 the numerator and denominator are both multiplied
# by exp(-s*(tp-t)) before the division.
#
# Under Bernoulli sampling with fraction rho the quantity the conditioning
# needs is the probability of leaving at least one SAMPLED descendant.  Writing
# Z for the number of descendants alive at tp, that is 1 - E[(1-rho)^Z], and
# with the geometric law of Z it evaluates to
#
#   p_rho(t) = rho*d / (d*E + rho*lam0*(1 - E)),   E = exp(-d*(tp-t)),
#
# with d = lam0 - mu0; the critical limit is rho / (1 + rho*lam0*(tp-t)).
# Both reduce to the complete-sampling forms at rho = 1, and p_rho(tp) = rho:
# a lineage alive at the present leaves a sampled descendant (itself) exactly
# when it is sampled.
#
# Under the Stadler reparameterisation lam' = rho*lam0, mu' = mu0 - lam0*(1-rho)
# the growth rate d is unchanged and the denominator above is exactly the
# complete-sampling denominator lam' - mu'*E, so
#
#   p_rho(t; lam0, mu0) = rho * p_1(t; rho*lam0, mu0 - lam0*(1 - rho)).
#
# The factor rho is not cosmetic: p_1 of anything is 1 at tp, and the terminal
# value here is rho.  The identity is checked to 7e-16 in
# tests/testthat/test-bdi-extended.R; it is a check, not the implementation.
#
# The rho >= 1 branch is the original code, untouched, so a complete-sampling
# call returns the same bits it returned before rho entered this function.
.bdi_p_cr <- function(t, lam0, mu0, tp, rho = 1) {
  d <- lam0 - mu0
  if (rho >= 1) {
    if (abs(d) <= 1e-12 * max(abs(lam0), abs(mu0))) {
      return(1 / (1 + lam0 * (tp - t)))
    }
    s  <- abs(d)
    em <- expm1(-s * (tp - t))            # in [-1, 0)
    if (d > 0) {
      # lam0 - mu0*exp(-d*tau) = d - mu0*expm1(-d*tau), both terms positive
      return(d / (d - mu0 * em))
    } else {
      # (mu0 - lam0)*exp(-s*tau) / (mu0 - lam0*exp(-s*tau))
      return(s * exp(-s * (tp - t)) / (s - lam0 * em))
    }
  }

  if (abs(d) <= 1e-12 * max(abs(lam0), abs(mu0))) {
    return(rho / (1 + rho * lam0 * (tp - t)))
  }
  s  <- abs(d)
  em <- expm1(-s * (tp - t))              # in [-1, 0)
  if (d > 0) {
    # d*E + rho*lam0*(1 - E) = d*(1 + em) - rho*lam0*em, both terms >= 0
    rho * d / (d * (1 + em) - rho * lam0 * em)
  } else {
    # numerator and denominator both multiplied by exp(-s*(tp - t))
    rho * s * exp(-s * (tp - t)) / (s - rho * lam0 * em)
  }
}


# --------------------------------------------------------------------------- #
#  Analytical integrals for exact BDI logg (CR only)                           #
# --------------------------------------------------------------------------- #

#' Exact integral of total BDI rate over [t1, t2] for CR with constant n, k.
#'
#' total_rate(t) = (n+2k)*lam0*(1-p(t)) + n*mu0/(1-p(t))
#'
#' Uses the identities:
#'   int lam0*(1-p) dt = mu0*(t2-t1) + ln[p(t1)/p(t2)]
#'   int mu0/(1-p) dt  = lam0*(t2-t1) + ln[(1-p(t1))/(1-p(t2))]
#'
#' Valid on both sides of lam0 = mu0.  Both log-ratios are evaluated with
#' s = |lam0 - mu0| in the exponent, so the argument of exp is never
#' positive and the mu0 > lam0 side does not overflow; writing the factors
#' as \code{s - r*expm1(-s*tau)} also removes the cancellation of
#' \code{lam0 - mu0*E} near lam0 = mu0.
#'
#' Under incomplete sampling the same two identities hold with p replaced by
#' its rho-sampled counterpart -- they follow from dp/dt = lam*p^2 - (lam-mu)*p
#' alone, which carries no terminal condition -- so only the closed forms of
#' the two log-ratios change.  With
#'
#'   A(t) = d*E + rho*lam0*(1 - E)          (so p = rho*d/A)
#'   B(t) = A(t) - rho*d                    (so 1 - p = B/A)
#'
#' and E = exp(-d*(tp - t)), the integrals are
#'
#'   int lam0*(1-p) dt = min(lam0,mu0)*(t2-t1) + ln[A(t2)/A(t1)]
#'   int mu0/(1-p)  dt = max(lam0,mu0)*(t2-t1) + ln[B(t1)/B(t2)]
#'
#' At rho = 1, B(tp) = 0 and the second integral diverges as t2 -> tp: a
#' doomed lineage must die before the present.  At rho < 1, B(tp) = d*(1-rho)
#' is non-zero and the integral is finite -- a doomed lineage may survive to
#' the present as an unsampled extant tip, which is the configuration
#' \code{.bdi_augment_one} then has to emit.
#'
#' @return The integral value (Inf when rho = 1, n > 0 and t2 == tp).
#' @keywords internal
.bdi_integral_cr <- function(t1, t2, n, k, lam0, mu0, tp, rho = 1) {
  if (t2 - t1 < 1e-15) return(0)
  if (rho < 1) return(.bdi_integral_cr_rho(t1, t2, n, k, lam0, mu0, tp, rho))
  d <- lam0 - mu0
  if (abs(d) <= 1e-12 * max(abs(lam0), abs(mu0))) {
    # Critical case: lam0 = mu0 (same switch as .bdi_p_cr)
    # p(t) = 1/(1+lam0*(tp-t)), 1-p = lam0*(tp-t)/(1+lam0*(tp-t))
    #   int lam0*(1-p) dt = lam0*(t2-t1) - ln(a1/a2)
    #   int mu0/(1-p)  dt = lam0*(t2-t1) + ln[(tp-t1)/(tp-t2)]
    a1 <- 1 + lam0 * (tp - t1)
    a2 <- 1 + lam0 * (tp - t2)
    I_lam <- lam0 * (t2 - t1) - log(a1 / a2)
    if (n > 0L) {
      if (tp - t2 < 1e-300) return(Inf)   # t2 at tp: divergent
      I_mu <- lam0 * (t2 - t1) + log((tp - t1) / (tp - t2))
    } else {
      I_mu <- 0
    }
    return((n + 2L * k) * I_lam + n * I_mu)
  }

  # Both branches in terms of s = |d|, so exp is only ever evaluated at a
  # non-positive argument.  With r1 the smaller of the two rates and r2 the
  # larger, exp(-s*tau) factors out of each log-ratio and contributes
  # -s*(t2-t1) to I_lam and +s*(t2-t1) to I_mu, turning the mu0 prefactor of
  # I_lam into r1 and the lam0 prefactor of I_mu into r2:
  #   I_lam = r1*(t2-t1) + ln[(s - r1*expm1(-s*tau2)) / (s - r1*expm1(-s*tau1))]
  #   I_mu  = r2*(t2-t1) + ln[expm1(-s*tau1) / expm1(-s*tau2)]
  # For lam0 > mu0 this is the shipped form with lam0 - mu0*E written as
  # s - mu0*expm1(-s*tau), a sum of two positive terms.
  s   <- abs(d)
  r1  <- min(lam0, mu0)
  r2  <- max(lam0, mu0)
  em1 <- expm1(-s * (tp - t1))
  em2 <- expm1(-s * (tp - t2))

  I_lam <- r1 * (t2 - t1) + log((s - r1 * em2) / (s - r1 * em1))

  if (n > 0L) {
    if (abs(em2) < 1e-300) return(Inf)   # t2 at tp: divergent (either sign of d)
    I_mu <- r2 * (t2 - t1) + log(em1 / em2)
  } else {
    I_mu <- 0
  }

  (n + 2L * k) * I_lam + n * I_mu
}


#' The rho < 1 branch of \code{\link{.bdi_integral_cr}}.
#'
#' Kept in its own function so that the complete-sampling path is the code
#' the validation study measured, byte for byte.
#'
#' A and B are written so that every exponential argument is non-positive and
#' every term of every sum is non-negative, on both sides of lam0 = mu0:
#'
#'   d > 0:  A = d*(1 + em) - rho*lam0*em,  em = expm1(-d*(tp - t))
#'           B = d*(1 - rho) + em*(d - rho*lam0)
#'   d < 0:  both multiplied by exp(-s*(tp - t)), s = |d|, and negated
#'           A = s - rho*lam0*em,           em = expm1(-s*(tp - t))
#'           B = s*(1 - rho)*exp(-s*(tp-t)) + em*(s + rho*lam0)
#'   the scale factor exp(-s*(tp - t)) contributes -s*(t2 - t1) to ln[A2/A1]
#'   and +s*(t2 - t1) to ln[B1/B2], which is exactly what turns the min/max
#'   prefactors around.
#'
#' B is 1 - p times A and p is in (0, 1], so B > 0 throughout for either sign
#' of d - rho*lam0.
#' @keywords internal
.bdi_integral_cr_rho <- function(t1, t2, n, k, lam0, mu0, tp, rho) {
  d <- lam0 - mu0
  if (abs(d) <= 1e-12 * max(abs(lam0), abs(mu0))) {
    # p = rho/(1 + rho*lam0*tau), 1 - p = (1 - rho + rho*lam0*tau)/(1 + ...)
    a1 <- 1 + rho * lam0 * (tp - t1); a2 <- 1 + rho * lam0 * (tp - t2)
    b1 <- 1 - rho + rho * lam0 * (tp - t1)
    b2 <- 1 - rho + rho * lam0 * (tp - t2)
    I_lam <- lam0 * (t2 - t1) + log(a2 / a1)
    I_mu  <- if (n > 0L) lam0 * (t2 - t1) + log(b1 / b2) else 0
    return((n + 2L * k) * I_lam + n * I_mu)
  }

  s   <- abs(d)
  r1  <- min(lam0, mu0)
  r2  <- max(lam0, mu0)
  em1 <- expm1(-s * (tp - t1))
  em2 <- expm1(-s * (tp - t2))

  if (d > 0) {
    A1 <- d * (1 + em1) - rho * lam0 * em1
    A2 <- d * (1 + em2) - rho * lam0 * em2
    B1 <- d * (1 - rho) + em1 * (d - rho * lam0)
    B2 <- d * (1 - rho) + em2 * (d - rho * lam0)
  } else {
    A1 <- s - rho * lam0 * em1
    A2 <- s - rho * lam0 * em2
    B1 <- s * (1 - rho) * exp(-s * (tp - t1)) - em1 * (s + rho * lam0)
    B2 <- s * (1 - rho) * exp(-s * (tp - t2)) - em2 * (s + rho * lam0)
  }

  I_lam <- r1 * (t2 - t1) + log(A2 / A1)
  I_mu  <- if (n > 0L) r2 * (t2 - t1) + log(B1 / B2) else 0

  (n + 2L * k) * I_lam + n * I_mu
}


#' Find next event time using exact time-change for CR.
#'
#' Given cumulative hazard H(s) = .bdi_integral_cr(t_cur, s, ...),
#' find t* such that H(t*) = U where U ~ Exp(1).
#'
#' @return Event time, or value > t_max if no event before boundary.
#' @keywords internal
.bdi_find_event_time_cr <- function(t_cur, U, n, k, lam0, mu0, tp, t_max,
                                    rho = 1) {
  # Quick check: will H reach U before t_max?
  # When n>0 and t_max==tp, H→Inf so event always exists -- but only at
  # complete sampling.  At rho < 1 the extinction integral is finite at tp
  # (a doomed lineage may survive as an unsampled extant tip), so the last
  # segment takes the ordinary branch and can return "no event".
  is_last_seg <- (rho >= 1) && (tp - t_max < 1e-10) && (n > 0L)
  if (!is_last_seg) {
    H_max <- .bdi_integral_cr(t_cur, t_max, n, k, lam0, mu0, tp, rho)
    if (U >= H_max) return(t_max + 1)  # no event
  }

  # Upper bound for root search
  t_hi <- if (is_last_seg) tp - 1e-14 else t_max
  f <- function(s) .bdi_integral_cr(t_cur, s, n, k, lam0, mu0, tp, rho) - U
  stats::uniroot(f, c(t_cur + 1e-15, t_hi), tol = 1e-13)$root
}

#' Can the BDI sampler be used for this model/link/rho?
#'
#' The sampler conditions on a rate that is a function of the lineage count
#' alone: \code{.bdi_p_cr} and \code{.bdi_solve_p_backward} solve one survival
#' probability \code{p(t)} shared by every lineage alive at \code{t}, and the
#' aggregate Gillespie of \code{.bdi_augment_one} moves \code{n_alive} lineages
#' with one pair of rates.  That is what an N-only model gives it.
#'
#' A D-dependent model (\code{model_bin[3] == 1}) does not: its rate depends on
#' each lineage's own pendant age \code{E_s}, so the survival probability is a
#' functional of the lineage's own history and no single \code{p(t)} exists.
#' The construction this file is built on does not extend there, and
#' D-dependent models stay on the thinning proposal.  The M covariate
#' (\code{model_bin[2]}) is a clade-level mean and is not in scope here either:
#' the sampler's own mean-field state would have to be closed on \code{P} as
#' well as \code{N}, which is written but not gated.
#'
#' Under the \code{gaussian} link the rate is
#' \code{beta_0 * exp(-(beta_N*N - 1)^2 / 2)}, still a function of \code{N}
#' alone.  For \code{cr} the covariate part is zero and the rate is the
#' constant \code{beta_0 * exp(-1/2)}, so constant rates on the gaussian link
#' are the constant-rates case under a reparameterisation and are sampled
#' exactly, as they are on the other two links.  For \code{dd} they are not:
#' \code{lambda(N)} rises to a peak at \code{N = 1/beta_N} and falls after it,
#' and the Picard iteration of \code{\link{.bdi_iterate}} does not contract on
#' the far side.  Measured on a 20-tip tree over a 4 x 8 x 4 grid in
#' \code{(beta_0, beta_N, rho)} with a 200-sweep budget -- more than three times
#' what this function allows -- it failed to reach \code{tol} in 36 of 128
#' cells, with \code{delta} running to 104.  That is divergence, not slow
#' convergence.  The failing region is not an interval that could be excluded:
#' at \code{rho = 0.5, beta_0 = 1} the iteration converges at
#' \code{beta_N = 0.08}, fails at 0.1, and converges again at 0.2, 0.3 and 0.5;
#' and the whole pattern moves with \code{rho}.  Where it does converge it can
#' take 164 sweeps.  \code{dd} on the gaussian link therefore stays on the
#' thinning proposal.  The linear and exponential \code{dd} rates are monotone
#' in \code{N} and reached the fixed point in every cell of the same sweep, in
#' at most 27 iterations at \code{rho = 1}.  They are slower below it: an
#' exponential \code{dd} cell needs 63 sweeps at \code{rho = 0.2}, and
#' converging linear \code{dd} cells have needed 126 and 134, which is what
#' sets the \code{rho < 1} budget.
#'
#' Incomplete sampling is in scope for every model and link the gate otherwise
#' accepts.  \code{p(t)} becomes the probability of leaving a \emph{sampled}
#' descendant (\code{.bdi_p_cr}, \code{.bdi_solve_p_backward}, both of which
#' take \code{rho}), and \code{.bdi_to_tree_df} emits the \code{5e10} sentinel
#' for a missing lineage still alive at the present, so the draws carry
#' unsampled extant lineages and \code{N(t)} counts them.
#'
#' @param model_bin Length-3 binary model vector \code{c(use_N, use_M, use_D)}.
#' @param link Integer link code: 0 linear, 1 exponential, 2 gaussian.
#' @param rho Sampling fraction in \code{(0, 1]}. Default \code{1}.
#' @keywords internal
.bdi_supported <- function(model_bin, link, rho = 1.0) {
  model_bin <- as.integer(model_bin)
  link      <- as.integer(link)
  # A rho above 1 is out of range and the callers reject it (.check_rho); it
  # reaches the C++ layer as 1, so it is read here as complete sampling too
  # rather than silently routed elsewhere.
  rho_ok    <- is.numeric(rho) && length(rho) == 1L && is.finite(rho) &&
    rho > 0 && rho <= 1 + 1e-12
  if (!rho_ok || !(length(model_bin) %in% c(3L, 4L))) return(FALSE)
  model_bin <- .pad_model_bin(model_bin)
  # No M covariate and no D covariate: each would close the mean-field state
  # on something the iteration does not carry.  ED is admitted on the linear
  # and exponential links: its clade mean is P-hat/N-hat, which the iteration
  # already solves for, so the proposal can read the mean-field ED rate and
  # leave the per-lineage departure from it to the importance weights.
  if (model_bin[2L] != 0L || model_bin[3L] != 0L) return(FALSE)
  if (model_bin[4L] != 0L) return(link %in% c(0L, 1L))
  if (link %in% c(0L, 1L)) return(TRUE)
  # gaussian: constant rates only (see above).
  link == 2L && model_bin[1L] == 0L
}

#' Why \code{\link{.bdi_supported}} refused, as a phrase for a message.
#'
#' Returns \code{NULL} when the sampler is available.  The reasons are ordered
#' so that the first thing wrong with the call is the one reported.
#' @inheritParams .bdi_supported
#' @keywords internal
.bdi_unsupported_reason <- function(model_bin, link, rho = 1.0) {
  model_bin <- as.integer(model_bin)
  link      <- as.integer(link)
  if (!(is.numeric(rho) && length(rho) == 1L && is.finite(rho) &&
        rho > 0 && rho <= 1 + 1e-12))
    return(sprintf("rho = %s (outside (0, 1])", format(rho)))
  if (!(length(model_bin) %in% c(3L, 4L)))
    return("a model vector that is not length 3 or 4")
  model_bin <- .pad_model_bin(model_bin)
  if (model_bin[4L] != 0L && !(link %in% c(0L, 1L)))
    return(paste0("an ED-dependent model on the gaussian link: the mean-field ",
                  "ED rate the proposal is built on is not available there"))
  if (model_bin[3L] != 0L)
    return(paste0("a D-dependent model: the BDI conditional distribution is ",
                  "built on one survival probability p(t) shared by every ",
                  "lineage alive at t, and a D-model's rate depends on each ",
                  "lineage's own pendant age, so no single p(t) exists"))
  if (model_bin[2L] != 0L)
    return(paste0("an M-dependent model: the sampler's mean-field state would ",
                  "have to be closed on the clade-mean pendant age as well as ",
                  "on N"))
  if (!(link %in% c(0L, 1L, 2L)))
    return(sprintf("link code %d", link))
  if (link == 2L && model_bin[1L] != 0L)
    return(paste0("a diversity-dependent model on the gaussian link: ",
                  "lambda(N) = beta_0*exp(-(beta_N*N - 1)^2/2) is not monotone ",
                  "in N and the sampler's Picard iteration does not converge ",
                  "on the far side of its peak"))
  NULL
}

# Compute speciation rate from 8-param vector + model.
#
# The gaussian link is the package's quadratic-exponential rate
# (inst/include/model.hpp): beta_0 is the peak rate and the intercept is
# excluded from the quadratic, so the covariate part alone enters it.  For an
# N-only model that leaves beta_N*N, and for cr it leaves 0 -- a constant rate
# beta_0*exp(-1/2).
.bdi_lam <- function(pars8, N, P, E, model_bin, link) {
  if (link == 2L) {
    eta_cov <- pars8[2] * N + pars8[3] * P + pars8[4] * E
    return(pars8[1] * exp(-0.5 * (eta_cov - 1)^2))
  }
  eta <- pars8[1] + pars8[2] * N + pars8[3] * P + pars8[4] * E
  if (link == 0L) max(0, eta) else exp(eta)
}

.bdi_mu <- function(pars8, N, P, E, model_bin, link) {
  if (link == 2L) {
    eta_cov <- pars8[6] * N + pars8[7] * P + pars8[8] * E
    return(pars8[5] * exp(-0.5 * (eta_cov - 1)^2))
  }
  eta <- pars8[5] + pars8[6] * N + pars8[7] * P + pars8[8] * E
  if (link == 0L) max(0, eta) else exp(eta)
}


# --------------------------------------------------------------------------- #
#  Backward-forward iteration (DD / general endogenous)                        #
# --------------------------------------------------------------------------- #

#' Solve survival probability p(t) backward from tp to 0.
#' ODE: dp/dt = lam(t)*p^2 - (lam(t)-mu(t))*p, p(tp) = rho.
#' Rates use the mean-field covariates (N̂, P̂, Ê) — supports DD/PD/EP and
#' mixed models. For CR (all model_bin = 0) the covariate arguments are
#' ignored inside .bdi_lam/.bdi_mu.
#'
#' The ODE itself carries no sampling fraction; \code{rho} enters only as the
#' terminal condition \code{p(tp) = rho}, which is the probability that a
#' lineage alive at the present leaves a sampled descendant (itself).  That is
#' the whole of the incomplete-sampling change on the backward side, and it
#' reproduces the closed form of \code{.bdi_p_cr} under constant rates.
#' @keywords internal
.bdi_solve_p_backward <- function(pars8, model_bin, link, bt, tp,
                                  Nhat_fun, Phat_fun, Ehat_fun, t_grid,
                                  rho = 1) {
  n_grid <- length(t_grid)
  t_rev  <- rev(t_grid)
  p_vals <- numeric(n_grid)
  p      <- min(rho, 1.0)

  for (i in seq_along(t_rev)) {
    ti <- t_rev[i]
    p_vals[n_grid - i + 1L] <- p
    if (i < length(t_rev)) {
      dt <- t_rev[i] - t_rev[i + 1L]
      f <- function(t, p) {
        Nh <- Nhat_fun(t); Ph <- Phat_fun(t); Eh <- Ehat_fun(t)
        la <- .bdi_lam(pars8, Nh, Ph, Eh, model_bin, link)
        mu <- .bdi_mu(pars8, Nh, Ph, Eh, model_bin, link)
        la * p^2 - (la - mu) * p
      }
      # RK4 backward step
      k1 <- f(ti, p)
      k2 <- f(ti - dt / 2, p - dt / 2 * k1)
      k3 <- f(ti - dt / 2, p - dt / 2 * k2)
      k4 <- f(ti - dt, p - dt * k3)
      p  <- p - dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
      p  <- min(max(p, 0), 1)
    }
  }
  p_vals
}


#' Solve the forward mean ODE jointly for (m, P_miss).
#'
#' State: s = (m, P_miss) with
#'   dm/dt       = (lam_bdi - mu_bdi) * m + nu      (self-consistent BDI mean)
#'   dP_miss/dt  = m                                (mean-field: each alive
#'                                                   missing lineage extends
#'                                                   its pendant edge at rate 1)
#'
#' Rates lam, mu are evaluated at (N̂, P̂, Ê) to support DD/PD/EP + mixed.
#' @return list(m_vals, pm_vals) on t_grid.
#' @keywords internal
.bdi_solve_mean_forward <- function(pars8, model_bin, link, bt, tp,
                                    p_fun, Nhat_fun, Phat_fun, Ehat_fun,
                                    t_grid) {
  bt <- sort(bt)
  m  <- 0; pm <- 0
  m_vals  <- numeric(length(t_grid))
  pm_vals <- numeric(length(t_grid))

  k_of_t <- function(t) 2L + sum(bt <= t)

  for (i in seq_along(t_grid)) {
    ti <- t_grid[i]
    m_vals[i]  <- m
    pm_vals[i] <- pm
    if (i < length(t_grid)) {
      dt <- t_grid[i + 1L] - t_grid[i]
      f <- function(t, state) {
        mi <- state[1L]; pmi <- state[2L]
        Nh <- Nhat_fun(t); Ph <- Phat_fun(t); Eh <- Ehat_fun(t)
        la <- .bdi_lam(pars8, Nh, Ph, Eh, model_bin, link)
        mu <- .bdi_mu(pars8, Nh, Ph, Eh, model_bin, link)
        p   <- p_fun(t); omp <- 1 - p
        if (omp < 1e-14) return(c(-mi * 1e6, mi))
        la_bdi <- la * omp
        mu_bdi <- mu / omp
        k      <- k_of_t(t)
        nu     <- 2 * k * la * omp
        c((la_bdi - mu_bdi) * mi + nu, mi)
      }
      s  <- c(m, pm)
      k1 <- f(ti,            s)
      k2 <- f(ti + dt / 2,   s + dt / 2 * k1)
      k3 <- f(ti + dt / 2,   s + dt / 2 * k2)
      k4 <- f(ti + dt,       s + dt     * k3)
      s  <- s + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
      m  <- max(s[1L], 0)
      pm <- max(s[2L], 0)
    }
  }
  list(m_vals = m_vals, pm_vals = pm_vals)
}


#' Solve the forward mean + second-moment ODE (Gaussian closure).
#'
#' Augments `.bdi_solve_mean_forward` with variance state needed to cancel
#' the process-level Jensen residual E[P/N] − P̂/N̂ ≈ −Cov(P,N)/N̂² + P̂ Var(N)/N̂³.
#' State: s = (m, P_miss, Var(N), Cov(N,P)) with
#'   dm/dt         = (lam_bdi − mu_bdi) m + nu
#'   dP_miss/dt    = m
#'   dVar(N)/dt    = 2(lam_bdi − mu_bdi) Var(N) + (lam_bdi + mu_bdi) m + nu
#'   dCov(N,P)/dt  = (lam_bdi − mu_bdi) Cov(N,P) + Var(N)
#'
#' The Var(N) ODE is the standard jump-diffusion variance infusion for a BDI
#' process; the Cov(N,P) ODE comes from the deterministic coupling dP = N dt
#' between events. Cross-jumps at deaths (ΔN = −1, ΔP = −age) are dropped
#' at this order — they are second-order corrections to cNP and would require
#' tracking the age distribution.
#' @return list(m_vals, pm_vals, vN_vals, cNP_vals) on t_grid.
#' @keywords internal
.bdi_solve_meanvar_forward <- function(pars8, model_bin, link, bt, tp,
                                       p_fun, Nhat_fun, Phat_fun, Ehat_fun,
                                       t_grid) {
  bt <- sort(bt)
  m <- 0; pm <- 0; vN <- 0; cNP <- 0
  m_vals   <- numeric(length(t_grid))
  pm_vals  <- numeric(length(t_grid))
  vN_vals  <- numeric(length(t_grid))
  cNP_vals <- numeric(length(t_grid))

  k_of_t <- function(t) 2L + sum(bt <= t)

  for (i in seq_along(t_grid)) {
    ti <- t_grid[i]
    m_vals[i]   <- m
    pm_vals[i]  <- pm
    vN_vals[i]  <- vN
    cNP_vals[i] <- cNP
    if (i < length(t_grid)) {
      dt <- t_grid[i + 1L] - t_grid[i]
      f <- function(t, state) {
        mi   <- state[1L]; pmi  <- state[2L]
        vNi  <- state[3L]; cNPi <- state[4L]
        Nh <- Nhat_fun(t); Ph <- Phat_fun(t); Eh <- Ehat_fun(t)
        la <- .bdi_lam(pars8, Nh, Ph, Eh, model_bin, link)
        mu <- .bdi_mu(pars8, Nh, Ph, Eh, model_bin, link)
        p   <- p_fun(t); omp <- 1 - p
        if (omp < 1e-14) return(c(-mi * 1e6, mi, 0, 0))
        la_bdi <- la * omp
        mu_bdi <- mu / omp
        k      <- k_of_t(t)
        nu     <- 2 * k * la * omp
        drift  <- la_bdi - mu_bdi
        infuse <- la_bdi + mu_bdi
        c(drift * mi + nu,                           # dm
          mi,                                        # dP_miss
          2 * drift * vNi + infuse * mi + nu,        # dVar(N)
          drift * cNPi + vNi)                        # dCov(N,P)
      }
      s  <- c(m, pm, vN, cNP)
      k1 <- f(ti,          s)
      k2 <- f(ti + dt / 2, s + dt / 2 * k1)
      k3 <- f(ti + dt / 2, s + dt / 2 * k2)
      k4 <- f(ti + dt,     s + dt     * k3)
      s   <- s + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
      m   <- max(s[1L], 0)
      pm  <- max(s[2L], 0)
      vN  <- max(s[3L], 0)
      cNP <- s[4L]
    }
  }
  list(m_vals = m_vals, pm_vals = pm_vals,
       vN_vals = vN_vals, cNP_vals = cNP_vals)
}


#' Observed-tree phylogenetic diversity under the emphasis tip_start=0
#' convention.
#'
#' The BDI tree data frame (.bdi_to_tree_df) sets tip_start = 0 for every
#' observed lineage, so their pendant PD at t is k(t) * t, with k counting the
#' two crown lineages.  Missing species carry their actual birth time.
#' P(t) = P_obs(t) + P_miss(t).
#' @keywords internal
.bdi_P_obs <- function(t, bt) {
  k <- 2L + sum(bt <= t)
  k * t
}

#' Backward-forward iteration for self-consistent BDI rates.
#'
#' Mean-field state: N̂(t) = k(t) + m̂(t), P̂(t) = P_obs(t) + P̂_miss(t),
#' Ê(t) = P̂(t) / N̂(t).  All three enter the rate function; inactive
#' covariates (pars8 slopes pinned to 0) make the extension transparent
#' for CR/DD and activate PD/EP/mixed models.
#' The sweep budget is 60 because the iteration slows as \code{rho} falls: the
#' unsampled extant lineages raise \code{N̂}, which feeds back into the rates
#' that produced them.  On the 20-tip DD tree at \code{K = 20} the fixed point
#' is reached in 6 sweeps at \code{rho = 1}, 7 at 0.9, 14 at 0.5 and 27 at 0.2.
#' The old budget of 20 therefore returned a non-converged mean field at small
#' \code{rho} -- \code{delta = 2.9e-3} against \code{tol = 1e-4} in that last
#' case -- and the loop still breaks on \code{tol}, so raising the cap cannot
#' change a run that already converged.
#'
#' @return List with p_fun, Nhat_fun, Phat_fun, Ehat_fun, and the fixed-point
#'   diagnostics converged / delta / iterations.
#' @keywords internal
.bdi_iterate <- function(pars8, model_bin, link, bt, tp,
                         max_iter = NULL, tol = 1e-4, n_grid = 500,
                         use_gaussian_closure = TRUE, rho = 1,
                         pd_mode = c("pendant", "faith"),
                         damping = "auto") {
  pd_mode <- match.arg(pd_mode)
  # Under-relaxation.  The sweep is a Picard iteration x <- F(x), and where its
  # map expands it diverges: on a 416-tip tree at rates whose equilibrium is
  # critical the residual reaches 305 after 20 sweeps and 329 after 200, so a
  # bigger budget makes it worse (audit H108).  Taking a partial step,
  # x <- (1 - w) x + w F(x), contracts a map whose expansion is mild enough,
  # at the cost of needing more sweeps when it was already contracting.
  #
  # damping = 1 is the plain iteration and is what every caller got before.
  # "auto" starts there and halves w whenever the residual grows, which costs
  # nothing where the map already contracts and is the only setting that needs
  # no knowledge of the map.
  auto_damp <- identical(damping, "auto")
  omega <- if (auto_damp) 1 else as.numeric(damping)
  if (!auto_damp && (!is.finite(omega) || omega <= 0 || omega > 1))
    stop(".bdi_iterate: damping must be \"auto\" or a number in (0, 1].",
         call. = FALSE)
  omega_min <- 1 / 64
  # The budget is 20 at complete sampling, which is what the validation study
  # measured: raising it changes the mean field of any run that had NOT reached
  # tolerance within 20, and with it the proposal and the estimate (measured
  # fhat gaps of 8.1e-4, 1.9e-5 and 11.64 nats on ordinary dd/linear cells).
  # Below rho = 1 the fixed point is slower -- exponential dd needs 63 sweeps
  # at rho = 0.2 and converging linear-dd cells 126 and 134 -- so the budget
  # there is set from the measurement rather than inherited.
  if (is.null(max_iter)) max_iter <- if (rho >= 1) 20L else 200L
  bt     <- sort(bt)
  t_grid <- seq(0, tp, length.out = n_grid)

  k_of_t <- function(t) 2L + sum(bt <= t)
  k_vals     <- sapply(t_grid, k_of_t)
  # Two conventions, because two covariates need two different quantities.
  #
  #   "pendant"  every observed lineage is dated from the crown (tip_start = 0
  #              in the frame), so its pendant edge at t is t and the clade's
  #              is k(t) * t.  This is what the D covariate is defined on and
  #              what the sampler has always used.
  #   "faith"    the sum of the tree's branch lengths, which is what the
  #              fair-proportion ED of the lineages sums to.  It grows at the
  #              lineage count, dPD/dt = k(t), so PD(t) = integral of k -- not
  #              k(t) * t, which counts every lineage's whole history as if it
  #              had been alive since the crown and overstates PD severalfold.
  P_obs_vals <- if (pd_mode == "faith") {
    2 * t_grid + vapply(t_grid, function(t) sum(pmax(0, t - bt)), 1)
  } else k_vals * t_grid

  # Initial guess: m = 0, P_miss = 0 ⇒ N̂ = k, P̂ = P_obs, Ê = P̂/N̂.
  Nhat_vals <- k_vals
  Phat_vals <- P_obs_vals
  Ehat_vals <- ifelse(Nhat_vals > 0, Phat_vals / Nhat_vals, 0)
  vN_vals   <- numeric(length(t_grid))
  cNP_vals  <- numeric(length(t_grid))

  # Scale for convergence check (normalize P̂ by a stable magnitude).
  Pscale <- max(max(Phat_vals), 1)

  p_vals <- rep(0, length(t_grid))
  delta  <- NA_real_
  prev_delta <- Inf
  n_iter_used <- 0L

  for (iter in seq_len(max_iter)) {
    n_iter_used <- iter
    Nhat_fun <- stats::approxfun(t_grid, Nhat_vals, rule = 2)
    Phat_fun <- stats::approxfun(t_grid, Phat_vals, rule = 2)
    Ehat_fun <- stats::approxfun(t_grid, Ehat_vals, rule = 2)

    p_vals <- .bdi_solve_p_backward(pars8, model_bin, link, bt, tp,
                                    Nhat_fun, Phat_fun, Ehat_fun, t_grid,
                                    rho = rho)
    p_fun  <- stats::approxfun(t_grid, p_vals, rule = 2)

    if (use_gaussian_closure) {
      mean_res <- .bdi_solve_meanvar_forward(pars8, model_bin, link, bt, tp,
                                             p_fun, Nhat_fun, Phat_fun, Ehat_fun,
                                             t_grid)
      m_vals   <- mean_res$m_vals
      pm_vals  <- mean_res$pm_vals
      vN_vals  <- mean_res$vN_vals
      cNP_vals <- mean_res$cNP_vals
    } else {
      mean_res <- .bdi_solve_mean_forward(pars8, model_bin, link, bt, tp,
                                          p_fun, Nhat_fun, Phat_fun, Ehat_fun,
                                          t_grid)
      m_vals   <- mean_res$m_vals
      pm_vals  <- mean_res$pm_vals
      vN_vals  <- rep(0, length(t_grid))
      cNP_vals <- rep(0, length(t_grid))
    }

    Nhat_new <- k_vals + m_vals
    Phat_new <- P_obs_vals + pm_vals
    # Second-order Taylor correction: E[P/N] ≈ P̂/N̂ − Cov(P,N)/N̂² + P̂·Var(N)/N̂³
    # (zero under plain mean-field closure ⇒ backwards-compatible).
    Ehat_raw <- ifelse(Nhat_new > 0,
                       Phat_new / Nhat_new
                       - cNP_vals / pmax(Nhat_new, 1e-12)^2
                       + Phat_new * vN_vals / pmax(Nhat_new, 1e-12)^3,
                       0)
    # Clip to [0, tp]: E_s is a pendant age.
    Ehat_new <- pmin(pmax(Ehat_raw, 0), tp)

    # The residual of the map, not the step that is taken: it is zero at a
    # fixed point whatever the relaxation is, so `tol` means the same thing at
    # every omega.
    delta_N <- max(abs(Nhat_new - Nhat_vals))
    delta_P <- max(abs(Phat_new - Phat_vals)) / Pscale
    delta   <- max(delta_N, delta_P)

    if (auto_damp && is.finite(prev_delta) && delta > prev_delta &&
        omega > omega_min) {
      omega <- max(omega / 2, omega_min)
    }
    prev_delta <- delta

    if (omega >= 1) {
      Nhat_vals <- Nhat_new
      Phat_vals <- Phat_new
      Ehat_vals <- Ehat_new
    } else {
      Nhat_vals <- (1 - omega) * Nhat_vals + omega * Nhat_new
      Phat_vals <- (1 - omega) * Phat_vals + omega * Phat_new
      # E-hat is a function of the pair, so it is recomputed from the relaxed
      # pair rather than relaxed itself, which would leave the three mutually
      # inconsistent.
      Eh <- ifelse(Nhat_vals > 0,
                   Phat_vals / Nhat_vals
                   - cNP_vals / pmax(Nhat_vals, 1e-12)^2
                   + Phat_vals * vN_vals / pmax(Nhat_vals, 1e-12)^3,
                   0)
      Ehat_vals <- pmin(pmax(Eh, 0), tp)
    }
    Pscale    <- max(max(Phat_vals), 1)

    if (delta < tol) break
  }

  Nhat_fun <- stats::approxfun(t_grid, Nhat_vals, rule = 2)
  Phat_fun <- stats::approxfun(t_grid, Phat_vals, rule = 2)
  Ehat_fun <- stats::approxfun(t_grid, Ehat_vals, rule = 2)
  p_fun    <- stats::approxfun(t_grid, p_vals, rule = 2)

  list(p_fun    = p_fun,
       Nhat_fun = Nhat_fun,
       Phat_fun = Phat_fun,
       Ehat_fun = Ehat_fun,
       vN_vals  = vN_vals,
       cNP_vals = cNP_vals,
       t_grid   = t_grid,
       # The Picard iteration is not proved to contract, so a caller reads
       # these three rather than assuming a fixed point was reached.  Measured
       # on a 20-tip tree: the linear and exponential dd rates, which are
       # monotone in N, converged in every cell of a sweep over slope and rho
       # (rho down to 0.2), in 3 to 15 sweeps.  The gaussian dd rate
       # lambda(N) = beta_0*exp(-(beta_N*N - 1)^2/2) peaks at N = 1/beta_N and
       # falls after it; the map N -> rate -> N then oscillates and delta grew
       # to 20 in 25 of 84 cells.  That is why .bdi_supported refuses dd on the
       # gaussian link outright rather than trusting this flag per theta: the
       # failing region moves with rho, and an MCEM run walks through it.
       converged = isTRUE(delta < tol),
       delta     = delta,
       iterations = n_iter_used)
}


# --------------------------------------------------------------------------- #
#  BDI Gillespie augmentation (one tree)                                       #
# --------------------------------------------------------------------------- #

#' Draw one augmented tree via the BDI process.
#'
#' Under CR: uses exact time-change method with analytical cumulative hazard.
#'   The event times are drawn from the correct inhomogeneous BDI process,
#'   giving zero-variance IS weights (logf - logg = constant).
#' Under DD: uses approximate Gillespie with piecewise-constant rates.
#'   IS weights have nonzero variance (importance sampling, not exact).
#'
#' At \code{rho < 1} a missing lineage still alive at tp is not a rejection but
#' an \emph{unsampled extant} lineage: the conditioning event is "leaves no
#' sampled descendant", and at the present that is satisfied by being alive and
#' unsampled, an event of probability \code{1 - rho}.  Those lineages are
#' returned in \code{$unsampled} (birth times) and written with the \code{5e10}
#' sentinel by \code{\link{.bdi_to_tree_df}}, so \code{N(t)} counts them.
#'
#' @return List with \code{$reason}: \code{"accepted"} (then also
#'   \code{$species}, \code{$unsampled}, \code{$n_alive_at_tp}, \code{$logg}),
#'   \code{"max_missing"} (more than \code{max_missing} missing lineages
#'   drawn) or \code{"survivor"} (a missing lineage still alive at tp; only
#'   reachable at \code{rho = 1}, where such a tree has f = 0).
#' @keywords internal
#' Attachment weights: each lineage\'s own rate, around the aggregate
#'
#' The proposal\'s aggregate rate already carries the clade mean of the
#' covariate, so a lineage\'s own rate is that aggregate plus its own
#' departure from the mean, \code{lam + beta * (c_s - mean(c))} -- not
#' \code{lam + beta * c_s}, which shifts every lineage instead of spreading
#' them and leaves the mean in twice.
#'
#' Weights are floored away from zero rather than at it.  A lineage the
#' model gives a positive rate must stay reachable, or the proposal no
#' longer covers the target and the estimate is biased rather than merely
#' noisy.  The floor is a thousandth of the aggregate, which costs a
#' bounded amount of weight variance and keeps the support.
#'
#' @param cov Per-lineage covariate values (here the pendant age).
#' @param lam The aggregate rate the mean field supplies.
#' @param beta The covariate\'s coefficient.
#' @return A positive weight per lineage, averaging \code{lam}.
#' @keywords internal
.attach_w <- function(cov, lam, beta) {
  n <- length(cov)
  if (n == 0L) return(numeric(0))
  if (!is.finite(lam) || lam <= 0) return(rep(1, n))
  v <- lam + beta * (cov - mean(cov))
  pmax(v, lam * 1e-3)
}


.bdi_augment_one <- function(bt, pars8, model_bin, link, tp,
                             p_fun = NULL, Nhat_fun = NULL,
                             Phat_fun = NULL, Ehat_fun = NULL,
                             max_missing = 1e4L, rho = 1,
                             track_parents = FALSE, first_aug_id = 0L,
                             step_max = Inf, obs_pid = NULL, beta_ed = 0) {
  bt <- sort(bt)
  is_cr <- all(model_bin == 0L)
  complete <- (rho >= 1)

  lam0 <- .bdi_lam(pars8, 0, 0, 0, model_bin, link)
  mu0  <- .bdi_mu(pars8, 0, 0, 0, model_bin, link)

  boundaries <- c(0, bt, tp)
  cur_la_raw <- NA_real_
  alive   <- numeric(0)
  n_alive <- 0L
  species <- list()
  n_total <- 0L
  logg    <- 0

  # Parent recording.  The lineage a birth attaches to is read off the same
  # uniform that chose the event, so the random stream is untouched and every
  # draw is the draw the sampler made before (dev/bdi_invariance.R).
  alive_id  <- integer(0)     # id of each missing lineage alive now
  alive_par <- integer(0)     # the lineage it was born from
  sp_id     <- integer(0)     # ids, in the order species[] records deaths
  sp_par    <- integer(0)
  next_id   <- as.integer(first_aug_id)
  # The observed lineages alive when k of them are: the two crown lineages,
  # then the daughter of each observed branching event.  This is the naming
  # the C++ layer uses and .observed_parent_id() reports in.
  obs_ids <- function(k) c(-2L, -3L, if (k > 2L) seq_len(k - 2L) - 1L else integer(0))

  # Rate-proportional attachment.  The proposal's *timing* has to be
  # mean-field -- the conditioned waiting time is built on one total rate and
  # one shared survival probability p(t), and there is no per-lineage p(t).
  # The *attachment* does not: which lineage a birth joins is a discrete
  # choice at a known instant, so it can be drawn in proportion to that
  # lineage's own rate.  The two decisions are separable, and only the second
  # one is cheap to do exactly.
  #
  # The per-lineage quantity used here is the pendant age t - tip_start, which
  # is what ED is mostly made of and is free to carry.  It is a proxy for ED,
  # not ED, and the weight prices the difference like every other one.
  #
  # Measured paired against the uniform draw -- same tree, same seed, 300
  # draws, 4 trees x 3 seeds per cell -- the gain tracks the effect size,
  # which is what the parent-identity term in the weight does too:
  #
  #             25 tips   50    100    200
  #   -0.35 l0    1.66   1.63   1.35   1.58     uniform ESS  9,  10,  8,  2
  #   -0.15 l0    0.76   1.08   0.62   0.98     uniform ESS 105, 101, 76, 70
  #
  # It helps where the sampler is starving and costs a little where it is
  # already rich, so it is the default.  At a weak effect the lineages barely
  # differ, the tilt is mostly the proxy's own error, and uniform is as good
  # or better -- pass attach = "uniform" there if the draw is cheap anyway.
  attach_ed <- isTRUE(track_parents) && is.finite(beta_ed) && beta_ed != 0
  # Tip start of each observed lineage, carried forward: a lineage is pendant
  # from its last split, so this is updated as each observed event is passed,
  # never precomputed over events that have not happened yet.
  ts_obs <- new.env(parent = emptyenv())
  assign("-2", 0, envir = ts_obs); assign("-3", 0, envir = ts_obs)
  note_event <- function(e) {              # observed event e (1-based into bt)
    if (!attach_ed || is.null(obs_pid) || e < 1L || e > length(obs_pid)) return(invisible())
    assign(as.character(obs_pid[e]), bt[e], envir = ts_obs)   # the splitter
    assign(as.character(e - 1L),     bt[e], envir = ts_obs)   # its daughter
    invisible()
  }
  ts_at <- function(k) {
    vapply(obs_ids(k), function(i) {
      nm <- as.character(i)
      if (exists(nm, envir = ts_obs, inherits = FALSE)) get(nm, envir = ts_obs) else 0
    }, 0)
  }
  # draw an index in proportion to v, using the uniform already in hand
  pick <- function(u, v) {
    cs <- cumsum(v)
    tot <- cs[length(cs)]
    if (!is.finite(tot) || tot <= 0) return(min(length(v), as.integer(floor(u * length(v))) + 1L))
    min(length(v), sum(cs < u * tot) + 1L)
  }

  k_of_t <- function(t) 2L + sum(bt <= t)

  for (seg_idx in seq_len(length(boundaries) - 1L)) {
    t0 <- boundaries[seg_idx]
    t1 <- boundaries[seg_idx + 1L]
    k  <- k_of_t(t0 + 1e-12)
    t  <- t0
    if (seg_idx >= 2L) note_event(seg_idx - 1L)   # the split that opened this segment

    while (t < t1) {
      if (is_cr) {
        # ── Exact time-change method for CR ──
        U <- stats::rexp(1)
        t_star <- .bdi_find_event_time_cr(t, U, n_alive, k, lam0, mu0, tp, t1,
                                          rho)

        if (t_star >= t1) {
          # No event before boundary
          logg <- logg - .bdi_integral_cr(t, t1, n_alive, k, lam0, mu0, tp, rho)
          break
        }

        # Event at t_star: survival integral contribution = -U
        logg <- logg - U
        t    <- t_star

        # Compute exact rates at event time for type selection
        p    <- .bdi_p_cr(t, lam0, mu0, tp, rho)
        omp  <- 1 - p
        la   <- lam0 * omp
        mu_r <- mu0 / max(omp, 1e-15)
        nu   <- 2 * k * lam0 * omp
        total <- n_alive * (la + mu_r) + nu

      } else {
        # ── Approximate Gillespie for DD / PD / EP / mixed ──
        # Rates use the mean-field covariates from the backward-forward
        # iteration. Under DD only N̂ is active; PD activates P̂;
        # EP activates Ê = P̂/N̂ (mean isolation time).
        p   <- p_fun(t)
        omp <- 1 - p
        if (omp < 1e-12) break
        Nh <- Nhat_fun(t); Ph <- Phat_fun(t); Eh <- Ehat_fun(t)
        la_raw <- .bdi_lam(pars8, Nh, Ph, Eh, model_bin, link)
        mu_raw <- .bdi_mu(pars8, Nh, Ph, Eh, model_bin, link)
        cur_la_raw <- la_raw      # the aggregate the attachment tilts around
        la   <- la_raw * omp
        mu_r <- mu_raw / max(omp, 1e-15)
        nu   <- 2 * k * la_raw * omp
        total <- n_alive * (la + mu_r) + nu

        if (total < 1e-15) break
        dt <- stats::rexp(1, total)
        # The rates are held at their value at t for the whole step, so a step
        # long enough for the lineage count to move is charged at rates that
        # no longer hold.  step_max re-evaluates them on a mesh: a step cut
        # short carries no event, only its own survival term.  At step_max =
        # Inf the cut is the segment boundary and this is the plain Gillespie
        # step it was.
        cap <- min(step_max, t1 - t)
        if (dt > cap) {
          logg <- logg - total * cap
          t <- t + cap
          if (t >= t1 - 1e-13) break
          next
        }
        t <- t + dt
        logg <- logg - total * dt
      }

      # ── Event type selection (shared) ──
      # logg uses per-lineage rates (labeled density) to match C++ logf:
      #   logf has log(λ) per speciation = per-lineage rate
      #   so logg must use per-lineage BDI rates, not total rates.
      r <- stats::runif(1) * total
      if (r < n_alive * la) {
        logg    <- logg + log(la)           # per-species birth rate
        if (track_parents) {
          if (attach_ed) {
            # in proportion to each lineage's own rate, through its pendant
            # age; the uniform that chose the event is reused as the variate,
            # so the timing stream is untouched
            v  <- .attach_w(t - alive, cur_la_raw, beta_ed)
            u  <- r / (n_alive * la)
            j  <- pick(u, v)
            # the attachment is no longer 1/n_alive, and logg says so
            sv <- sum(v)
            if (is.finite(sv) && sv > 0)
              logg <- logg + log(n_alive * v[j] / sv)
          } else {
            # every alive missing lineage carries the same rate la, so r/la is
            # already a uniform index into them
            j <- min(n_alive, as.integer(floor(r / la)) + 1L)
          }
          alive_par <- c(alive_par, alive_id[j])
          alive_id  <- c(alive_id, next_id); next_id <- next_id + 1L
        }
        alive   <- c(alive, t)
        n_alive <- n_alive + 1L
        n_total <- n_total + 1L
      } else if (r < n_alive * la + nu) {
        logg    <- logg + log(la)           # per-lineage immigration rate
        if (track_parents) {
          if (attach_ed) {
            v  <- .attach_w(t - ts_at(k), cur_la_raw, beta_ed)
            u  <- (r - n_alive * la) / nu
            j  <- pick(u, v)
            sv <- sum(v)
            if (is.finite(sv) && sv > 0)
              logg <- logg + log(k * v[j] / sv)
          } else {
            # the k observed lineages share nu equally
            j <- min(k, as.integer(floor((r - n_alive * la) / (nu / k))) + 1L)
          }
          alive_par <- c(alive_par, obs_ids(k)[j])
          alive_id  <- c(alive_id, next_id); next_id <- next_id + 1L
        }
        alive   <- c(alive, t)
        n_alive <- n_alive + 1L
        n_total <- n_total + 1L
      } else if (n_alive > 0L) {
        logg <- logg + log(mu_r)            # per-species death rate
        idx  <- sample.int(n_alive, 1L)
        species[[length(species) + 1L]] <- c(alive[idx], t)
        if (track_parents) {
          sp_id  <- c(sp_id, alive_id[idx]); sp_par <- c(sp_par, alive_par[idx])
          alive_id <- alive_id[-idx]; alive_par <- alive_par[-idx]
        }
        alive   <- alive[-idx]
        n_alive <- n_alive - 1L
      }

      if (n_total > max_missing) return(list(reason = "max_missing"))
    }
  }

  # At complete sampling the exact CR process leaves no survivor (the
  # extinction hazard diverges at tp) and the approximate DD Gillespie
  # occasionally does; such a tree has f = 0, so it is rejected.
  # At rho < 1 a survivor is an unsampled extant lineage and is kept.
  if (complete && n_alive > 0L) return(list(reason = "survivor"))

  list(reason = "accepted", species = species,
       unsampled = alive, n_alive_at_tp = n_alive, logg = logg,
       sp_id = sp_id, sp_par = sp_par,
       uns_id = alive_id, uns_par = alive_par)
}


# --------------------------------------------------------------------------- #
#  Convert BDI output to emphasis tree data frame                              #
# --------------------------------------------------------------------------- #

#' Convert BDI species list to the tree data frame format used by emphasis.
#' Columns: brts, n, t_ext, pd, tip_start, id, parent_id
#'
#' The parent recorded for a lineage born at \code{t} is the last observed
#' branching at or before \code{t}.  A lineage born before the first observed
#' branching has none, and \code{max(0L, ...)} makes it observed node 0, which
#' is not yet born then; \code{\link{.aug_to_Ltable}} refuses that attachment
#' rather than building a \code{tas} with a negative edge.  See the
#' \code{.augment_tree_bdi} documentation for what that costs and what fixing
#' it would take.
#' @keywords internal
.bdi_to_tree_df <- function(species, bt, tp, unsampled = numeric(0),
                            obs_parent_id = NULL, obs_focal_ts = NULL,
                            aug_id = NULL, aug_par = NULL,
                            uns_id = NULL, uns_par = NULL) {
  # With the observed topology named, the frame carries the real parent of
  # every node -- observed ones from the caller, augmented ones from the draw
  # that made them -- instead of the last-branching guess below, and gains the
  # two columns the covariate sweep reads (audit finding H45).
  topology <- !is.null(obs_parent_id)
  bt_sorted <- sort(bt)
  n_obs     <- length(bt_sorted)
  n_aug     <- length(species)
  n_uns     <- length(unsampled)
  # obs + (spec + ext) per extinct missing + spec per unsampled extant
  # + closing node at tp
  n_total   <- n_obs + 2L * n_aug + n_uns + 1L

  # Pre-allocate vectors
  v_brts      <- numeric(n_total)
  v_t_ext     <- numeric(n_total)
  v_tip_start <- numeric(n_total)
  v_id        <- integer(n_total)
  v_parent_id <- integer(n_total)
  v_focal_ts  <- rep(-1, n_total)      # ts_unknown
  v_clade     <- rep(if (topology) 1L else 0L, n_total)

  # Observed speciation nodes
  idx <- seq_len(n_obs)
  v_brts[idx]      <- bt_sorted
  v_t_ext[idx]     <- 1e11   # t_ext_tip
  v_tip_start[idx] <- 0
  v_id[idx]        <- seq(0L, n_obs - 1L)
  v_parent_id[idx] <- if (topology) as.integer(obs_parent_id) else -1L
  if (topology && !is.null(obs_focal_ts)) v_focal_ts[idx] <- as.numeric(obs_focal_ts)

  # The last observed branching at or before `birth`; see the note above and
  # in .augment_tree_bdi on what this assignment costs (audit finding H45).
  parent_of <- function(birth) {
    parent <- which(bt_sorted <= birth) - 1L
    if (length(parent) == 0L) 0L else max(0L, max(parent))
  }

  # Augmented species: speciation + extinction nodes
  if (n_aug > 0L) {
    off <- n_obs
    for (i in seq_along(species)) {
      sp    <- species[[i]]
      birth <- sp[1]; death <- sp[2]
      sid   <- if (topology) as.integer(aug_id[i]) else as.integer(n_obs + i - 1L)
      parent <- if (topology) as.integer(aug_par[i]) else parent_of(birth)
      j <- off + 2L * (i - 1L)
      # Speciation node
      v_brts[j + 1L]      <- birth
      v_t_ext[j + 1L]     <- death
      v_tip_start[j + 1L] <- birth
      v_id[j + 1L]        <- sid
      v_parent_id[j + 1L] <- parent
      # Extinction node
      v_brts[j + 2L]      <- death
      v_t_ext[j + 2L]     <- 0     # t_ext_extinct
      v_tip_start[j + 2L] <- birth
      v_id[j + 2L]        <- sid
      v_parent_id[j + 2L] <- parent
    }
  }

  # Unsampled extant lineages (rho < 1): a speciation node carrying the
  # 5e10 sentinel and no extinction node -- they are alive at the present, so
  # they raise N(t) from their birth onward and never lower it again.  The
  # C++ scorer reads 5e10 as is_unsampled (inst/include/model_helpers.hpp):
  # the birth contributes log(lambda) like any other, the lineage contributes
  # no log(mu), and Model::loglik charges it log(1 - rho).
  if (n_uns > 0L) {
    off <- n_obs + 2L * n_aug
    for (i in seq_along(unsampled)) {
      birth <- unsampled[i]
      sid   <- if (topology) as.integer(uns_id[i]) else as.integer(n_obs + n_aug + i - 1L)
      j <- off + i
      v_brts[j]      <- birth
      v_t_ext[j]     <- 5e10       # t_ext_unsampled
      v_tip_start[j] <- birth
      v_id[j]        <- sid
      v_parent_id[j] <- if (topology) as.integer(uns_par[i]) else parent_of(birth)
    }
  }

  # Closing node at tp (present-day marker, matches C++ convention)
  v_brts[n_total]      <- tp
  v_t_ext[n_total]     <- 1e11
  v_tip_start[n_total] <- 0
  v_id[n_total]        <- -1L
  v_parent_id[n_total] <- -1L

  # Sort by brts
  ord <- order(v_brts)
  v_brts      <- v_brts[ord]
  v_t_ext     <- v_t_ext[ord]
  v_tip_start <- v_tip_start[ord]
  v_id        <- v_id[ord]
  v_parent_id <- v_parent_id[ord]
  v_focal_ts  <- v_focal_ts[ord]
  v_clade     <- v_clade[ord]

  # Compute n: n[i] = lineage count during [brts[i-1], brts[i]).
  # C++ convention: n[i] = n_after(event i-1) = n[i-1] + type(event i-1).
  v_n    <- numeric(n_total)
  n_cur  <- 2
  v_n[1] <- n_cur
  if (n_total >= 2L) {
    for (i in 2L:n_total) {
      n_cur  <- n_cur + if (v_t_ext[i - 1L] == 0) -1 else 1
      v_n[i] <- n_cur
    }
  }

  out <- data.frame(
    brts      = v_brts,
    n         = v_n,
    t_ext     = v_t_ext,
    pd        = rep(0, n_total),
    tip_start = v_tip_start,
    id        = v_id,
    parent_id = v_parent_id,
    stringsAsFactors = FALSE
  )
  if (!topology) return(out)
  # The sweep reads the convention off the last node, so every row carries it.
  out$focal_tip_start <- v_focal_ts
  out$clade           <- v_clade
  out[c("brts", "n", "t_ext", "pd", "tip_start", "focal_tip_start",
        "clade", "id", "parent_id")]
}


# --------------------------------------------------------------------------- #
#  Main BDI augmentation wrapper (replaces .augment_tree_internal)             #
# --------------------------------------------------------------------------- #

#' BDI augmentation: draw sample_size augmented trees.
#'
#' Returns the same structure as the augment_trees() C++ function,
#' list(trees, logf, logg, weights, fhat), plus the draw counts.
#'
#' Draws are attempted until \code{sample_size} have completed or the budget
#' of \code{5 * sample_size} attempts is spent.  A draw is rejected when its
#' missing-lineage count exceeds \code{max_missing} (counted in
#' \code{n_rejected_max_missing}) or, under DD, when a missing lineage is
#' still alive at the present (counted in \code{n_rejected}; such a tree has
#' f = 0 at rho = 1).  Completed draws with a non-finite log-weight are
#' dropped from the returned tree set (the M-step set) but stay in the fhat
#' denominator with weight zero.
#'
#' fhat = log(sum_w / n_valid) + max_lw + log(acc), where n_valid counts every
#' completed draw and acc = n_valid / (n_valid + n_rejected) is the acceptance
#' rate of the survivor channel.  max_missing overflows are excluded from the
#' denominator, as in the thinning E-step (S_completed in src/E_step.cpp).
#'
#' The log(acc) correction is derived at \code{rho = 1}, where a surviving
#' missing lineage has f = 0 and the rejection is a property of the target.
#' At \code{rho < 1} an unsampled survivor is a legitimate configuration and is
#' emitted rather than rejected, so \code{n_rejected} is 0, \code{acc} is 1 and
#' the factor is inert.  \code{acc} is the
#' stop-at-\code{n_valid} plug-in estimate of the acceptance rate, which is
#' biased at small \code{sample_size}: at the function's own default
#' \code{sample_size = 1L} a draw with three survivor rejections moves fhat
#' by log(1/4) = -1.39.  At the MCEM default of a few hundred draws the bias
#' is of order 1e-3.
#'
#' @section Every pre-first-split lineage is attached to observed node 0:
#' \code{.bdi_to_tree_df} records, as the parent of an augmented lineage born
#' at \code{t}, the last observed branching at or before \code{t}
#' (\code{max(which(bt_sorted <= birth) - 1L)}).  For a lineage born before the
#' \strong{first} observed branching there is none, and the expression falls
#' back to observed node 0 -- a lineage that is not yet born at \code{t}.  The
#' sampler has no notion of the two crown lineages, which is where such a birth
#' really attaches, and it never records the reserved crown ids the thinning
#' sampler uses.
#'
#' The consequence used to be silent and wrong: \code{tas} carried an edge of
#' negative length wherever it happened.  On seven trees measured here, the
#' parent build turned all 200 draws into a \code{phylo} and 5 to 200 of them
#' had a negative edge, the shortest \code{-7.15}.
#'
#' \code{\link{.aug_to_Ltable}} now refuses an attachment older than its
#' recorded parent and returns \code{NULL}, so those draws produce no tree at
#' all.  The failure is loud where it used to be silent, and the throughput is
#' the price: on the same seven trees, \code{simulate_tree(method = "bdi")}
#' built 0, 106, 124, 129, 153, 156 and 168 trees out of 200, with the refused
#' draws in every case exactly those carrying a lineage born before the first
#' observed branching.  A tree whose first observed branching is late loses
#' most of its draws; one whose crown splits early loses few.
#'
#' The parent assignment itself is unchanged and still wrong.  Fixing it means
#' giving the BDI sampler the two crown lineages and the same uniform draw over
#' labelled attachments the thinning sampler makes (audit finding H45).
#'
#' @return List: \code{trees}, \code{logf}, \code{logg}, \code{weights}
#'   (finite-weight draws only), \code{fhat}, \code{n_valid} (completed
#'   draws), \code{n_nonfinite} (completed draws dropped for a non-finite
#'   log-weight), \code{n_attempts}, \code{n_rejected} (survivors at tp),
#'   \code{n_rejected_max_missing}, \code{acc}.  \code{trees} is empty when
#'   no draw completed or every completed draw has a non-finite log-weight;
#'   a warning is issued when fewer than \code{sample_size} draws completed.
#' @keywords internal
.augment_tree_bdi <- function(tree,
                              pars,
                              model_bin   = c(0L, 0L, 0L),
                              sample_size = 1L,
                              max_missing = 1e4,
                              link        = 0L,
                              rho         = 1.0,
                              use_gaussian_closure = TRUE,
                              topology    = NULL,
                              mesh        = NULL,
                              attach      = c("rate", "uniform"),
                              damping     = "auto") {
  attach <- match.arg(attach)
  brts  <- .extract_brts(tree)
  mb4   <- .pad_model_bin(model_bin)
  # Accept compact, 8-element (no ED) or 10-element (with ED) pars
  pars_full <- if (length(pars) %in% c(8L, 10L)) pars else .expand_pars(pars, mb4)
  use_ed    <- mb4[4L] != 0L
  if (use_ed && length(pars_full) < 10L)
    stop(".augment_tree_bdi: an ED model needs the 10-element parameter vector",
         call. = FALSE)
  # The proposal is the mean-field one: every lineage alive at t is given the
  # rate the clade's average lineage has there.  The clade mean of
  # fair-proportion ED is the pendant PD per lineage, P-hat/N-hat, which the
  # iteration already carries as E-hat -- so the ED coefficient enters the
  # proposal through the same slot the mean isolation time does.  What each
  # lineage's own ED does to its rate is left to the importance weights.
  pars8 <- pars_full[1:8]
  mb_prop <- mb4[1:3]
  if (use_ed) {
    pars8[4L] <- pars_full[9L]
    pars8[8L] <- pars_full[10L]
    mb_prop[3L] <- 1L
  }
  tp    <- brts[1L]

  # Convert emphasis brts (crown-age first, decreasing) to forward-time bt
  # emphasis brts: tp, t_{n-1}, ..., t_1 (decreasing, present = 0)
  # forward time bt: tp - brts (increasing from 0)
  bt <- sort(tp - brts[-1L])

  # The observed topology, when the branching times carry it: which lineage
  # splits at each event and when that lineage became a pendant tip.  With it
  # the draw records the real parent of every augmented birth, and the tree it
  # returns can be scored by a per-lineage covariate.  Without it the frame is
  # what it has always been -- see .bdi_to_tree_df.
  n_obs   <- length(bt)
  # The attributes carry one entry per branching time, the last a sentinel for
  # the present; the events are the first n_obs, in forward-time order, which
  # is the order bt is in.
  ok_len  <- function(x) length(x) == n_obs || length(x) == n_obs + 1L
  obs_pid <- utils::head(.pid(brts), n_obs)
  obs_pts <- utils::head(.pts(brts), n_obs)
  use_top <- if (is.null(topology)) {
    ok_len(.pid(brts)) && ok_len(.pts(brts))
  } else isTRUE(topology)
  if (use_top && (length(obs_pid) != n_obs || length(obs_pts) != n_obs))
    stop(".augment_tree_bdi: topology requested but the branching times carry ",
         "no parent_id/parent_tip_start of the right length", call. = FALSE)

  is_cr <- all(mb_prop == 0L)

  # How finely the approximate Gillespie re-evaluates its frozen rates, as a
  # number of steps per crown age.  NULL leaves the step uncapped.
  step_max <- if (is.null(mesh) || !is.finite(mesh) || mesh <= 0) Inf else tp / mesh

  # Solve BDI rates
  p_fun <- Nhat_fun <- Phat_fun <- Ehat_fun <- NULL
  mf_converged <- NA
  mf_delta     <- NA_real_
  if (!is_cr) {
    sol <- .bdi_iterate(pars8, mb_prop, link, bt, tp,
                        use_gaussian_closure = use_gaussian_closure,
                        rho = rho,
                        pd_mode = if (use_ed) "faith" else "pendant",
                        damping = damping)
    p_fun    <- sol$p_fun
    Nhat_fun <- sol$Nhat_fun
    Phat_fun <- sol$Phat_fun
    Ehat_fun <- sol$Ehat_fun
    mf_converged <- isTRUE(sol$converged)
    mf_delta     <- sol$delta
    # The mean-field covariates the whole DD proposal is built on are the
    # fixed point of .bdi_iterate.  When the iteration stops on max_iter
    # instead, the proposal is built on whatever the last sweep produced.  It
    # is still a valid proposal -- the weights carry the error -- but the
    # caller is told rather than left to assume a fixed point was reached.
    if (!mf_converged) warning(sprintf(paste0(
      "BDI augmentation: the mean-field iteration did not converge ",
      "(delta = %.2e after %d sweeps); the proposal uses the last sweep."),
      sol$delta, sol$iterations), call. = FALSE)
  }

  # Draw augmented trees.  A draw is rejected on max_missing overflow or,
  # under DD (approximate Gillespie), when a missing lineage survives to tp.
  # Attempts continue until sample_size draws have completed or the budget
  # is spent; the two rejection channels are counted separately.
  trees      <- vector("list", sample_size)
  logg       <- numeric(sample_size)
  n_valid    <- 0L
  n_attempts <- 0L
  n_rej_surv <- 0L
  n_rej_mm   <- 0L
  max_tries  <- 5L * sample_size

  for (attempt in seq_len(max_tries)) {
    if (n_valid >= sample_size) break
    n_attempts <- n_attempts + 1L
    aug <- .bdi_augment_one(bt, pars8, mb_prop, link, tp,
                            p_fun, Nhat_fun, Phat_fun, Ehat_fun,
                            max_missing, rho = rho,
                            track_parents = use_top, first_aug_id = n_obs,
                            step_max = step_max, obs_pid = obs_pid,
                            beta_ed = if (use_ed && attach == "rate") pars_full[9L] else 0)
    if (aug$reason == "survivor")    { n_rej_surv <- n_rej_surv + 1L; next }
    if (aug$reason == "max_missing") { n_rej_mm   <- n_rej_mm   + 1L; next }

    n_valid <- n_valid + 1L
    df <- if (use_top) {
      .bdi_to_tree_df(aug$species, bt, tp, unsampled = aug$unsampled,
                      obs_parent_id = obs_pid, obs_focal_ts = obs_pts,
                      aug_id = aug$sp_id, aug_par = aug$sp_par,
                      uns_id = aug$uns_id, uns_par = aug$uns_par)
    } else {
      .bdi_to_tree_df(aug$species, bt, tp, unsampled = aug$unsampled)
    }
    # The frame names the events; the sweep turns them into the state the
    # covariates are read from (tip_start, focal_tip_start and pendant PD).
    if (use_top) df <- eval_pendant_sweep(df)
    trees[[n_valid]] <- df
    logg[n_valid]     <- aug$logg
  }

  # Trim to actual count
  trees <- trees[seq_len(n_valid)]
  logg  <- logg[seq_len(n_valid)]

  if (n_valid < sample_size) {
    warning(sprintf(paste0(
      "BDI augmentation: %d of %d requested draws completed in %d attempts ",
      "(%d survivors at tp, %d over max_missing = %d)."),
      n_valid, sample_size, n_attempts, n_rej_surv, n_rej_mm,
      as.integer(max_missing)), call. = FALSE)
  }

  # Compute logf (model log-likelihood) via C++ eval_logf.
  # eval_logf also returns a thinning-based logg — we discard it
  # and use the Gillespie-accumulated logg from above instead.
  if (n_valid > 0L) {
    ev   <- eval_logf(pars_full, trees,
                      model = as.integer(mb4),
                      link  = as.integer(link),
                      rho   = as.numeric(rho))
    logf <- ev$logf
  } else {
    logf <- numeric(0)
  }

  weights <- logf - logg

  # A completed draw with a non-finite log-weight (logf = -Inf when a rate
  # is zero on the augmented tree) has weight zero: it stays in the fhat
  # denominator and leaves the M-step set.
  finite      <- is.finite(weights)
  n_nonfinite <- sum(!finite)
  trees   <- trees[finite]
  logf    <- logf[finite]
  logg    <- logg[finite]
  weights <- weights[finite]

  # Survivor rejections have f = 0 at rho = 1, so the completed draws are a
  # sample from g / P(accept | theta) and log(acc) restores that factor.
  # max_missing overflows are left out of the denominator (E_step.cpp).
  acc <- if (n_valid + n_rej_surv > 0L) n_valid / (n_valid + n_rej_surv) else NA_real_
  # one denominator for both proposals -- see .is_summary
  summ <- .is_summary(weights, n_zero_weight = n_nonfinite, n_rejected = n_rej_surv)
  fhat <- if (is.na(summ$fhat)) -Inf else summ$fhat

  list(trees   = trees,
       logf    = logf,
       logg    = logg,
       weights = weights,
       fhat    = fhat,
       n_valid = n_valid,
       n_nonfinite = n_nonfinite,
       n_attempts  = n_attempts,
       n_rejected  = n_rej_surv,
       n_rejected_max_missing = n_rej_mm,
       acc     = acc,
       # NA under cr (no mean-field system is solved); TRUE/FALSE under dd.
       mf_converged = mf_converged,
       mf_delta     = mf_delta)
}


# --------------------------------------------------------------------------- #
#  BDI MCEM loop: E-step (BDI) + M-step (C++ nlopt)                           #
# --------------------------------------------------------------------------- #

#' MCEM with BDI augmentation.
#'
#' Drop-in replacement for .mcem_dynamic_fresh(). Uses BDI exact sampling
#' for the E-step and m_cpp() for the M-step.
#'
#' Under CR: zero-variance IS, so every E-step gives the same fhat regardless
#' of the augmented trees.  Convergence depends only on the M-step.
#'
#' Stopping rule.  With theta_k the M-step output of iteration k,
#' \deqn{\delta_k = \max_j |\theta_k[j] - \theta_{k-1}[j]| / \max(|\theta_{k-1}[j]|, \epsilon)}
#' (\code{.rel_change}) with \eqn{\epsilon} from \code{.rel_floor}, and the
#' run stops as "converged" once \eqn{\delta_k < tol} for \code{patience}
#' consecutive iterations.  The scale is the parameter itself, not the
#' search box, so the same fit stops at the same point whatever bounds the
#' user passes; the floor carries the units of the parameters, so it also
#' stops at the same point when the tree is measured in another unit of
#' time.  An iteration whose
#' M-step returned its starting point unchanged carries no information about
#' stability and resets the streak.  The trace also records the absolute
#' step \code{abs_step} and \code{drift}, the same relative displacement
#' measured over the last \code{patience} iterations: under pure Monte Carlo
#' noise \code{drift} is of the order of \code{sqrt(patience)} steps, under a
#' deterministic EM drift it is \code{patience} steps.
#'
#' The reported likelihood is a separate E-step at the final iterate, so
#' \code{fhat}, ESS, \code{final_IS} and \code{loglik_var} all describe
#' \code{pars}; that E-step is the last row of \code{mcem} (\code{m_step =
#' FALSE}, no step columns) and is not counted in \code{iterations}.  When
#' that E-step fails and the final iterate is not the last one whose E-step
#' succeeded, \code{pars} falls back to that iterate and the E-step is run
#' there: a theta the sampler cannot evaluate is not returned as an
#' estimate.
#'
#' Draw counts come from the E-step, which drops non-finite draws before
#' returning: \code{n_valid} is the number of completed draws and the fhat
#' denominator, \code{num_trees} the finite subset handed to the M-step, and
#' \code{n_nonfinite} the difference.  \code{rejected} and
#' \code{rejected_max_missing} are the two rejection channels.
#'
#' Eight E-step failures end the run with \code{stop_reason =
#' "e_step_failure"} whether or not they are consecutive, so an M-step that
#' keeps proposing an iterate the sampler cannot evaluate does not alternate
#' to \code{max_iter}.
#'
#' @return A list: \code{mcem} (trace, one row per iteration plus the final
#'   E-step), \code{pars}, \code{iterations} (number of completed E+M
#'   iterations), \code{stop_reason}, \code{loglik} (fhat at \code{pars},
#'   \code{NA} when the final E-step failed), \code{loglik_var},
#'   \code{final_IS}, \code{n_failed} (E- or M-step failures).
#' @keywords internal
.mcem_bdi <- function(brts,
                      pars,
                      sample_size,
                      max_missing,
                      lower_bound,
                      upper_bound,
                      max_iter,
                      xtol,
                      tol = 1e-2,
                      patience,
                      num_threads,
                      verbose    = FALSE,
                      conditional = NULL,
                      model      = c(0L, 0L, 0L),
                      link       = 0L,
                      max_time   = NULL,
                      rho        = 1.0,
                      stop_rule  = c("rel_change", "mc_error"),
                      mc_batches = 5L,
                      mc_z       = 1.0,
                      mc_grow    = 1.5,
                      max_draws  = NULL,
                      mesh       = NULL,
                      damping    = "auto") {
  stop_rule <- match.arg(stop_rule)
  if (is.null(max_draws) || !is.finite(max_draws)) max_draws <- as.integer(8L * sample_size)
  max_draws <- as.integer(max(max_draws, sample_size))
  prev_iterate <- pars

  if (inherits(brts, "phylo"))
    brts <- sort(ape::branching.times(brts), decreasing = TRUE)
  if (!is.numeric(brts)) stop("`brts` must be numeric or a `phylo` object.")
  if (!is.numeric(pars)) stop("`pars` must be numeric.")

  model_bin <- as.integer(model)
  link_int  <- as.integer(link)
  patience  <- max(1L, as.integer(patience))

  # Floor of the relative-change denominator: a parameter below the floor in
  # magnitude is measured against the floor, so a coordinate sitting at zero
  # does not turn every step into an infinite relative change.  The floor
  # carries the units of the parameters (.rel_floor), so the statistic is
  # the same when the tree and the rates are expressed in another unit of
  # time.  .mcem_dynamic_fresh uses the same two helpers.
  floor_val  <- .rel_floor(brts, link_int)
  rel_change <- function(new, old) .rel_change(new, old, floor_val)

  # One BDI E-step at `theta`; NULL on error or when no tree was drawn.
  e_step_at <- function(theta) {
    e <- tryCatch(
      .augment_tree_bdi(
        tree        = brts,
        pars        = theta,
        model_bin   = model_bin,
        sample_size = as.integer(sample_size),
        max_missing = as.integer(max_missing),
        link        = link_int,
        rho         = as.numeric(rho),
        mesh        = mesh,
        damping     = damping
      ),
      error = function(e) NULL
    )
    if (is.null(e) || length(e$trees) == 0L) return(NULL)
    e
  }

  # Trace row.  Step columns are NA for the final E-step (no M-step).
  trace_row <- function(theta, e, elapsed_e, m_time = 0, m_step = TRUE,
                        delta_max = NA_real_, abs_step = NA_real_,
                        drift = NA_real_, m_moved = NA) {
    par_df <- as.data.frame(as.list(
      stats::setNames(theta, paste0("par", seq_along(theta)))))
    cbind(par_df, data.frame(
      fhat        = e$fhat,
      delta_max   = delta_max,
      abs_step    = abs_step,
      drift       = drift,
      m_step      = m_step,
      m_moved     = m_moved,
      rejected    = .n0(e$n_rejected),
      rejected_max_missing = .n0(e$n_rejected_max_missing),
      # .augment_tree_bdi drops non-finite draws before returning and
      # reports the count; the second term covers a caller that does not.
      n_nonfinite = .n0(e$n_nonfinite) + sum(!is.finite(e$logf)),
      n_valid     = .n0(e$n_valid),
      num_trees   = length(e$trees),
      ESS         = .ess_from_lw(e$weights),
      # "mc_error" only: the step in units of its own Monte Carlo standard
      # error, and the draws the E-step used at that iteration.
      mc_z        = NA_real_,
      sample_size = sample_size,
      time        = elapsed_e * 1000 + m_time
    ))
  }

  streak      <- 0L
  fail_streak <- 0L
  n_efail     <- 0L            # E-step failures, consecutive or not
  n_failed    <- 0L
  prev_pars   <- pars          # last iterate whose E-step succeeded
  history     <- list(pars)    # theta_0, theta_1, ... for the drift column
  mcem        <- NULL
  n_iter      <- 0L
  stop_reason <- "max_iter"
  t0_mcem     <- proc.time()[3]

  for (i in seq_len(max_iter)) {
    # ── E-step: BDI augmentation ──
    t0_e  <- proc.time()[3]
    e_raw <- e_step_at(pars)
    elapsed_e <- as.numeric(proc.time()[3] - t0_e)

    if (is.null(e_raw)) {
      fail_streak <- fail_streak + 1L
      n_efail     <- n_efail + 1L
      n_failed    <- n_failed + 1L
      streak      <- 0L
      # Restart from the last iterate whose E-step succeeded; the box centre
      # is not a known-good point and may itself be where the sampler fails.
      pars <- prev_pars
      if (verbose) message(sprintf(
        "Iteration %d: E-step failed (%d consecutive, %d in total) - restarting from the last successful iterate",
        i, fail_streak, n_efail))
      # Eight failures end the run whether or not they are consecutive: an
      # M-step that keeps proposing the same unusable iterate alternates
      # success and failure forever and the streak never builds.
      if (fail_streak >= 8L || n_efail >= 8L) {
        stop_reason <- "e_step_failure"; break
      }
      next
    }

    # M-step set: every draw .augment_tree_bdi returns.  It has already
    # dropped the non-finite ones (they stay in its fhat denominator and
    # are counted in the trace); an E-step with nothing left is mapped to
    # NULL by e_step_at, so this guard is defensive.
    lw <- e_raw$weights
    if (length(lw) == 0L) {
      fail_streak <- fail_streak + 1L
      n_efail     <- n_efail + 1L
      n_failed    <- n_failed + 1L
      streak      <- 0L
      pars <- prev_pars
      if (verbose) message(sprintf(
        "Iteration %d: E-step returned no finite log-weight (%d consecutive)",
        i, fail_streak))
      if (fail_streak >= 8L || n_efail >= 8L) {
        stop_reason <- "e_step_failure"; break
      }
      next
    }
    fail_streak <- 0L
    prev_pars   <- pars

    # m_cpp objective: sum loglik(theta, tree_i) * w[i] with w as direct
    # multipliers.  BDI log-weights are constant and negative under CR, so
    # convert to self-normalised IS weights (mean 1, all positive).
    w_norm <- exp(lw - max(lw))
    w_norm <- w_norm / sum(w_norm) * length(w_norm)

    e_step <- list(
      trees                 = e_raw$trees,
      weights               = w_norm,
      rejected              = .n0(e_raw$n_rejected),
      rejected_overruns     = 0L,
      rejected_lambda       = 0L,
      rejected_zero_weights = .n0(e_raw$n_nonfinite) + sum(!is.finite(e_raw$logf)),
      time                  = elapsed_e * 1000,
      fhat                  = e_raw$fhat
    )

    # ── M-step: C++ optimisation ──
    m_result <- tryCatch(
      m_cpp(e_step     = e_step,
            init_pars  = pars,
            plugin     = "rpd1",
            lower_bound = lower_bound,
            upper_bound = upper_bound,
            xtol_rel   = xtol,
            num_threads = as.integer(num_threads),
            model      = model_bin,
            link       = link_int,
            rho        = as.numeric(rho),
            rconditional = conditional),
      error = function(e) NULL
    )

    if (is.null(m_result)) {
      fail_streak <- fail_streak + 1L
      n_failed    <- n_failed + 1L
      streak      <- 0L
      if (verbose) message(sprintf("Iteration %d: M-step failed", i))
      if (fail_streak >= 8L) { stop_reason <- "m_step_failure"; break }
      next
    }

    new_pars  <- as.numeric(m_result$estimates)
    m_moved   <- any(new_pars != pars)
    abs_step  <- max(abs(new_pars - pars))
    delta_max <- rel_change(new_pars, pars)
    n_iter    <- n_iter + 1L
    history[[n_iter + 1L]] <- new_pars
    drift <- if (n_iter > patience)
      rel_change(new_pars, history[[n_iter + 1L - patience]]) else NA_real_
    pars      <- new_pars

    mcem <- rbind(mcem, trace_row(
      pars, e_raw, elapsed_e, m_time = m_result$time,
      delta_max = delta_max, abs_step = abs_step, drift = drift,
      m_moved = m_moved))

    if (verbose) {
      message(sprintf(
        "Iteration %d: fhat=%.4f  delta=%.2e  step=%.2e  streak=%d/%d  ESS=%.0f/%d%s",
        i, e_raw$fhat, delta_max, abs_step, streak, patience,
        .ess_from_lw(e_raw$weights), length(e_raw$trees),
        if (m_moved) "" else "  (M-step returned its start)"))
    }

    # Check convergence.  An M-step that returned its start unchanged says
    # nothing about stability (the objective may have been undefined there),
    # so it does not count toward patience under either rule.
    #
    # "rel_change" is the package's original rule: `patience` consecutive
    # iterations whose relative parameter change is below `tol`, at a fixed
    # number of draws.  "mc_error" compares the step with its own Monte Carlo
    # standard error instead, estimated by re-running the M-step on batches of
    # the same draws; when the step is inside that noise the draws are grown,
    # and the run stops only once they have reached max_draws.  The estimate
    # returned is then the mean of the last `patience` iterates.  See
    # .mcem_dynamic_fresh for the references.
    if (stop_rule == "rel_change") {
      if (!m_moved) {
        streak <- 0L
      } else if (delta_max < tol) {
        streak <- streak + 1L
        if (streak >= patience) { stop_reason <- "converged"; break }
      } else {
        streak <- 0L
      }
    } else {
      se <- local({
        B <- as.integer(mc_batches); nn <- length(e_step$trees)
        if (B < 2L || nn < 2L * B) return(NULL)
        idx <- split(seq_len(nn), rep(seq_len(B), length.out = nn))
        est <- lapply(idx, function(ii) {
          if (!any(e_step$weights[ii] > 0)) return(NULL)
          eb <- e_step; eb$trees <- e_step$trees[ii]; eb$weights <- e_step$weights[ii]
          m <- tryCatch(m_cpp(eb, pars, "rpd1", lower_bound, upper_bound, xtol,
                              as.integer(num_threads), model = model_bin, link = link_int,
                              rho = as.numeric(rho), rconditional = conditional),
                        error = function(e) NULL)
          if (is.null(m)) NULL else as.numeric(m$estimates)[seq_along(pars)]
        })
        est <- do.call(rbind, Filter(Negate(is.null), est))
        if (is.null(est) || nrow(est) < 2L) return(NULL)
        apply(est, 2, stats::sd) / sqrt(nrow(est))
      })
      z <- if (is.null(se)) NA_real_ else max(abs(new_pars - prev_iterate) / pmax(se, .Machine$double.eps))
      mcem$mc_z[nrow(mcem)] <- z
      mcem$sample_size[nrow(mcem)] <- sample_size
      if (!m_moved || is.na(z)) {
        streak <- 0L
      } else if (z < mc_z) {
        streak <- streak + 1L
        if (sample_size < max_draws) {
          sample_size <- min(as.integer(ceiling(mc_grow * sample_size)), max_draws)
          streak <- 0L
          if (verbose) message(sprintf("  step within Monte Carlo error (z = %.2f): draws -> %d", z, sample_size))
        } else if (streak >= patience) { stop_reason <- "converged"; break }
      } else {
        streak <- 0L
      }
      prev_iterate <- new_pars
    }

    # Time budget
    if (!is.null(max_time)) {
      elapsed <- proc.time()[3] - t0_mcem
      if (elapsed > max_time) {
        stop_reason <- "time_budget"
        if (verbose) message(sprintf("Time budget reached (%.0fs)", elapsed))
        break
      }
    }
  }

  # Under "mc_error" the returned estimate is the mean of the last `patience`
  # iterates: at the fixed point they fluctuate around the maximiser with the
  # Monte Carlo error the rule measured, and averaging removes part of it.
  if (stop_rule == "mc_error" && n_iter >= 2L) {
    take <- history[seq(max(2L, n_iter + 1L - patience + 1L), n_iter + 1L)]
    take <- Filter(function(x) is.numeric(x) && length(x) == length(pars), take)
    if (length(take) >= 2L) pars <- colMeans(do.call(rbind, take))
  }

  # ── Final E-step at the returned iterate ──
  # The trace rows pair theta_k with fhat(theta_{k-1}); this E-step gives
  # fhat, ESS and the IS diagnostics at pars itself.
  t0_e    <- proc.time()[3]
  e_final <- e_step_at(pars)
  elapsed_e <- as.numeric(proc.time()[3] - t0_e)

  # An iterate the sampler cannot evaluate is not an estimate: fall back to
  # the last iterate whose E-step succeeded and report that one instead.
  if (is.null(e_final) && !identical(pars, prev_pars)) {
    n_efail  <- n_efail + 1L
    n_failed <- n_failed + 1L
    pars     <- prev_pars
    if (verbose) message(
      "Final E-step failed; returning the last iterate whose E-step succeeded.")
    t0_e      <- proc.time()[3]
    e_final   <- e_step_at(pars)
    elapsed_e <- as.numeric(proc.time()[3] - t0_e)
  }

  if (!is.null(e_final)) {
    mcem <- rbind(mcem, trace_row(pars, e_final, elapsed_e, m_step = FALSE))
  } else if (verbose) {
    message("Final E-step at the returned parameters failed; loglik is NA.")
  }

  loglik     <- NA_real_
  loglik_var <- NA_real_
  final_IS   <- NULL
  if (!is.null(e_final)) {
    loglik <- e_final$fhat
    lw     <- e_final$logf - e_final$logg
    # .augment_tree_bdi returns only draws with a finite log-weight.
    if (length(e_final$logf) >= 2L) {
      loglik_var <- .bootstrap_fhat_var(e_final$logf, e_final$logg,
                                        K = 2L, B = 200L)
    }
    final_IS <- list(
      logf  = e_final$logf,
      logg  = e_final$logg,
      lw    = lw,
      fhat  = e_final$fhat,
      ESS   = .ess_from_lw(lw),
      gap   = .jensen_gap(.ess_from_lw(lw),
                          length(lw) + .n0(e_final$n_nonfinite) +
                            .n0(e_final$n_rejected))$gap,
      gap_heavy = .jensen_gap(.ess_from_lw(lw),
                              length(lw) + .n0(e_final$n_nonfinite) +
                                .n0(e_final$n_rejected))$heavy,
      n_rejected = .n0(e_final$n_rejected),
      rejected_max_missing  = .n0(e_final$n_rejected_max_missing),
      rejected_zero_weights = .n0(e_final$n_nonfinite) +
                              sum(!is.finite(e_final$logf))
    )
  }

  list(
    mcem        = mcem,
    pars        = pars,
    iterations  = n_iter,
    stop_reason = stop_reason,
    loglik      = loglik,
    loglik_var  = loglik_var,
    final_IS    = final_IS,
    n_failed    = n_failed
  )
}

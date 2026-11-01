#' Constant-rate maximum-likelihood estimate of a reconstructed tree
#'
#' The likelihood of a crown tree under constant speciation \eqn{\lambda},
#' extinction \eqn{\mu} and sampling fraction \eqn{\rho} is available in
#' closed form (Nee et al. 1994; Stadler 2009).  With \eqn{r = \lambda - \mu}
#' and \eqn{x_1 > \dots > x_{n-1}} the branching times (crown first),
#' \deqn{p_1(t) = \frac{\rho r^2 e^{-rt}}{(\rho\lambda + (\lambda(1-\rho) - \mu) e^{-rt})^2},
#'       \qquad 1 - p_0(t) = \frac{\rho r}{\rho\lambda + (\lambda(1-\rho) - \mu) e^{-rt}},}
#' and, up to a constant,
#' \deqn{\log L = 2 \log p_1(x_1) + \sum_{i \ge 2} \{\log\lambda + \log p_1(x_i)\}
#'       - c \cdot 2 \log\{1 - p_0(x_1)\},}
#' where \eqn{c = 1} conditions on the survival of both crown lineages and
#' \eqn{c = 0} does not.  This is the likelihood the MCEM maximises for the
#' \code{cr} model, so its maximiser is the fixed point of that fit.
#'
#' @param brts Branching times (any order, crown age the largest).
#' @param cond Logical: condition on survival of the crown lineages.
#' @param rho Sampling fraction in (0, 1].
#' @return \code{c(lambda, mu)}, or \code{NULL} when the optimiser fails.
#' @keywords internal
.cr_mle <- function(brts, cond = FALSE, rho = 1) {
  x <- sort(as.numeric(brts), decreasing = TRUE)
  x <- x[x > 0]
  if (length(x) < 2L) return(NULL)
  nll <- function(th) {
    l <- exp(th[1L]); m <- exp(th[2L]); r <- l - m
    if (abs(r) < 1e-10) r <- 1e-10
    den <- rho * l + (l * (1 - rho) - m) * exp(-r * x)
    lp1 <- log(rho) + 2 * log(abs(r)) - r * x - 2 * log(abs(den))
    ll <- 2 * lp1[1L] + sum(log(l) + lp1[-1L])
    if (cond) ll <- ll - 2 * (log(rho) + log(abs(r)) - log(abs(den[1L])))
    if (!is.finite(ll)) return(1e100)
    -ll
  }
  n <- length(x) + 1L
  l0 <- max(log(n / 2) / max(x), 1e-3)
  best <- NULL
  for (f in c(0.1, 0.5, 0.9)) {
    o <- tryCatch(stats::optim(log(c(l0 / (1 - f), l0 * f / (1 - f))), nll,
                               control = list(reltol = 1e-12, maxit = 2000)),
                  error = function(e) NULL)
    if (!is.null(o) && (is.null(best) || o$value < best$value)) best <- o
  }
  if (is.null(best) || best$value >= 1e99) return(NULL)
  exp(best$par)
}

#' Default starting point of an MCEM fit
#'
#' The constant-rate fit of the tree (\code{.cr_mle}) placed in the model's
#' parameter vector: the two intercepts carry its \eqn{\lambda} and
#' \eqn{\mu} (log rates under the exponential link, rates under the linear
#' link) and every covariate coefficient is zero, so the starting model is
#' the constant-rate model.  For the \code{cr} model this is the maximiser of
#' the observed likelihood itself.  Each coordinate is then moved inside the
#' search box.  Returns \code{NULL} when the constant-rate fit fails; the
#' caller then starts at the middle of the box.
#'
#' @param brts Branching times.
#' @param model_bin Compact model vector.
#' @param link_int 0 linear, 1 exponential.
#' @param lower_bound,upper_bound Compact bounds.
#' @param cond Logical: the fit is conditioned on survival.
#' @param rho Sampling fraction.
#' @keywords internal
.default_init <- function(brts, model_bin, link_int, lower_bound, upper_bound,
                          cond = FALSE, rho = 1) {
  est <- .cr_mle(brts, cond = cond, rho = rho)
  if (is.null(est) || any(!is.finite(est)) || any(est <= 0)) return(NULL)
  n_pars <- length(lower_bound)
  init <- numeric(n_pars)
  icpt <- if (link_int == 1L) log(est) else est
  init[1L] <- icpt[1L]
  init[n_pars / 2L + 1L] <- icpt[2L]
  span <- upper_bound - lower_bound
  pmin(pmax(init, lower_bound + 1e-3 * span), upper_bound - 1e-3 * span)
}

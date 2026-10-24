# --------------------------------------------------------------------------- #
#  Survival of the clade by direct simulation, and the local surface of        #
#  log P(S) that estimate_rates(cond = "surface") conditions on.               #
# --------------------------------------------------------------------------- #
#
# Survival is crown survival, as everywhere in the package (train_GAM's labels,
# the cr formula of ?estimate_rates): both crown lineages leave a sampled
# descendant at the present.  A clade is one forward draw of
# simulate_div_tree_cpp() (one try), seeded by `seed + i` for the i-th clade,
# so two calls with the same seed draw the same clades at different
# parameters (common random numbers).  A clade that passes the lineage cap is
# counted as surviving, as train_GAM() counts it.  Under rho < 1 a clade that
# reached the present with n1 and n2 extant descendants of the two crown
# lineages contributes (1 - (1 - rho)^n1) (1 - (1 - rho)^n2), the probability
# that independent sampling of each tip with probability rho keeps one on both
# sides; the mean over clades estimates P(S) without a second random draw.

# Room above a seed for the clade indices of every call that offsets it.
.surv_seed_max <- .Machine$integer.max - 300000L
.surv_seed <- function(seed) {
  if (is.null(seed)) return(sample.int(.surv_seed_max, 1L))
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) || seed < 1)
    stop("'seed' must be a positive integer or NULL.", call. = FALSE)
  as.integer((as.numeric(seed) - 1) %% .surv_seed_max + 1)
}
.now <- function() proc.time()[[3L]]

# Run `expr` under set.seed(seed) and leave the caller's random stream as it was.
.with_seed <- function(seed, expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(if (had) assign(".Random.seed", old, envir = globalenv())
          else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
            rm(".Random.seed", envir = globalenv()))
  set.seed(seed)
  expr
}

# Survival values of clades seed + idx; stops early (returning the clades
# drawn) when `deadline` (a .now() value) passes.
.surv_values <- function(pars_full, model_bin, link_int, crown_age, rho, max_lin,
                         seed, idx, deadline = Inf) {
  v <- numeric(length(idx))
  done <- 0L
  for (j in seq_along(idx)) {
    if (j %% 64L == 1L && .now() > deadline) break
    s <- simulate_div_tree_cpp(pars_full, model_bin, crown_age, max_lin, 0L,
                               link_int, seed = as.integer((seed + idx[j] - 1) %% .Machine$integer.max + 1))
    v[j] <- if (identical(s$status, "extinct")) 0
            else if (rho >= 1 || identical(s$status, "too_large")) 1
            else {
              L <- s$Ltable
              ext <- L[, 4] == -1
              (1 - (1 - rho)^sum(ext & L[, 3] < 0)) * (1 - (1 - rho)^sum(ext & L[, 3] > 0))
            }
    done <- j
  }
  v[seq_len(done)]
}

# log P with the continuity floor (k + 1/2) / (n + 1) and its standard error;
# `s2` is the per-clade variance the standard error and the clade rule read.
.surv_summary <- function(v) {
  n <- length(v)
  k <- sum(v)
  P <- (k + 0.5) / (n + 1)
  s2 <- if (n > 1L) stats::var(v) else 0
  if (!is.finite(s2) || s2 <= 0) s2 <- P * (1 - P)
  c(logP = log(P), se = sqrt(s2 / n) / P, P = P, survivors = k, clades = n, s2 = s2)
}

# Clades for a standard error `sigma` of log P at a point with P and per-clade
# variance s2 (at rho = 1, s2 = P (1 - P) and this is (1 - P) / (P sigma^2)).
.surv_n_for <- function(P, s2, sigma, lo = 500L, hi = 40000L)
  as.integer(min(hi, max(lo, ceiling(s2 / (P^2 * sigma^2)))))

# Direct estimate at one parameter point (full layout).  With `sigma`, `n`
# clades first, then as many more as the rule asks for, up to `hi`.
.surv_direct <- function(pars_full, model_bin, link_int, crown_age, rho, max_lin,
                         seed, n = 2000L, sigma = NULL, hi = 40000L, deadline = Inf) {
  t0 <- .now()
  v <- .surv_values(pars_full, model_bin, link_int, crown_age, rho, max_lin,
                    seed, seq_len(n), deadline)
  if (!is.null(sigma) && length(v) == n) {
    a <- .surv_summary(v)
    n2 <- .surv_n_for(a[["P"]], a[["s2"]], sigma, lo = n, hi = max(hi, n))
    if (n2 > n)
      v <- c(v, .surv_values(pars_full, model_bin, link_int, crown_age, rho, max_lin,
                             seed, seq.int(n + 1L, n2), deadline))
  }
  if (!length(v)) return(c(logP = NA_real_, se = NA_real_, P = NA_real_, survivors = 0,
                           clades = 0, s2 = NA_real_, seconds = .now() - t0))
  c(.surv_summary(v), seconds = .now() - t0)
}

.surv_check_common <- function(crown_age, rho, max_lin) {
  if (!is.numeric(crown_age) || length(crown_age) != 1L || !is.finite(crown_age) || crown_age <= 0)
    stop("'crown_age' must be a positive number.", call. = FALSE)
  .check_rho(rho, arg = "rho")
  if (!is.numeric(max_lin) || length(max_lin) != 1L || !is.finite(max_lin) || max_lin < 2)
    stop("'max_lin' must be a number of at least 2.", call. = FALSE)
  invisible(NULL)
}

#' Probability that a clade survives, by direct simulation
#'
#' Estimates \eqn{\log P_\theta(S)}, the log-probability that both crown
#' lineages of a clade of crown age \code{crown_age} leave a sampled
#' descendant at the present, from forward simulations under any model and
#' link. This is the survival probability \code{estimate_rates(cond = )}
#' divides the likelihood by.
#'
#' Clade \eqn{i} is one forward draw of the simulator seeded with
#' \code{seed + i}, so calls with the same \code{seed} at different parameter
#' values draw from the same random numbers and their difference has a much
#' smaller error than either value. A clade that passes \code{max_lin}
#' lineages is counted as surviving, as \code{\link{train_GAM}} counts it.
#' Under \code{rho < 1} a clade that reaches the present with \eqn{n_1} and
#' \eqn{n_2} extant descendants of its two crown lineages contributes
#' \eqn{(1 - (1-\rho)^{n_1})(1 - (1-\rho)^{n_2})}, the probability that
#' sampling each extant tip with probability \eqn{\rho} keeps at least one
#' on each side. The estimate is \eqn{\log((k + 1/2)/(n + 1))}, \eqn{k} the
#' sum of the clades' contributions and \eqn{n} the clades drawn, and its
#' standard error is the delta-method one, \eqn{s / (\sqrt{n} P)} with
#' \eqn{s^2} the variance of the contributions (the binomial
#' \eqn{\sqrt{(1-P)/(nP)}} at \eqn{\rho = 1}).
#'
#' @param pars Parameter vector in the compact layout of \code{model} (as in
#'   \code{\link{simulate_tree}}).
#' @param model Model name, formula or binary vector, as in
#'   \code{\link{simulate_tree}}. Default \code{"cr"}.
#' @param crown_age Crown age of the clade.
#' @param rho Sampling fraction in \code{(0, 1]}. Default \code{1}.
#' @param link \code{"exponential"} (default) or \code{"linear"}.
#' @param n_clades Clades to draw; with \code{se}, the clades drawn first.
#'   Default \code{2000}.
#' @param se Target standard error of \eqn{\log P}, or \code{NULL} (default)
#'   to draw \code{n_clades} only. When set, the first \code{n_clades} give a
#'   provisional \eqn{P} and the clades are extended to the number that
#'   reaches \code{se} there, at most \code{max_clades}.
#' @param max_clades Cap on the clades drawn under \code{se}. Default
#'   \code{40000}.
#' @param max_lin Lineage cap of one clade. Default \code{1e4}.
#' @param seed Positive integer, or \code{NULL} (default) to draw one from
#'   R's generator.
#' @param max_time Seconds after which no further clade is drawn, or
#'   \code{NULL} (default) for no limit; the estimate then rests on the clades
#'   drawn.
#' @return A named numeric vector: \code{logP}, \code{se}, \code{P},
#'   \code{survivors} (the sum of the contributions), \code{clades},
#'   \code{seconds} and \code{seed}.
#' @seealso \code{\link{survival_surface}}, \code{\link{estimate_rates}}
#' @examples
#' # cr at rho = 1 against the closed form p^2 of ?estimate_rates
#' s <- survival_mc(log(c(0.5, 0.2)), model = "cr", crown_age = 3,
#'                  n_clades = 2000, seed = 1)
#' l <- 0.5; m <- 0.2; d <- l - m; T <- 3
#' c(mc = s[["logP"]], se = s[["se"]],
#'   exact = 2 * log(d / (l - m * exp(-d * T))))
#' @export
survival_mc <- function(pars, model = "cr", crown_age, rho = 1, link = "exponential",
                        n_clades = 2000L, se = NULL, max_clades = 40000L,
                        max_lin = 1e4, seed = NULL, max_time = NULL) {
  model_bin <- .resolve_model(model)
  link_int  <- .resolve_link(link)
  .surv_check_common(crown_age, rho, max_lin)
  expected_n <- 2L + 2L * sum(model_bin != 0L)
  if (!is.numeric(pars) || length(pars) != expected_n) stop(.pars_error_msg(model_bin, expected_n))
  if (!is.numeric(n_clades) || length(n_clades) != 1L || n_clades < 1)
    stop("'n_clades' must be a positive number.", call. = FALSE)
  if (!is.null(se) && (!is.numeric(se) || length(se) != 1L || !(se > 0)))
    stop("'se' must be a positive number or NULL.", call. = FALSE)
  seed <- .surv_seed(seed)
  deadline <- if (is.null(max_time)) Inf else .now() + max_time
  r <- .surv_direct(.expand_pars(as.numeric(pars), model_bin), model_bin, link_int,
                    crown_age, rho, as.integer(max_lin), seed, n = as.integer(n_clades),
                    sigma = se, hi = as.integer(max_clades), deadline = deadline)
  c(r[c("logP", "se", "P", "survivors", "clades", "seconds")], seed = seed)
}

# The quadratic basis in the scaled coordinates z and its gradient (E46-lib.R).
.qbasis <- function(z) {
  d <- length(z)
  out <- c(1, z)
  for (i in seq_len(d)) for (j in i:d) out <- c(out, z[i] * z[j])
  out
}
.qgrad_z <- function(cf, z) {
  d <- length(z)
  g <- cf[1L + seq_len(d)]
  k <- 1L + d
  for (i in seq_len(d)) for (j in i:d) {
    k <- k + 1L
    if (i == j) g[i] <- g[i] + 2 * cf[k] * z[i]
    else { g[i] <- g[i] + cf[k] * z[j]; g[j] <- g[j] + cf[k] * z[i] }
  }
  g
}
# A Latin hypercube on [-1, 1]^d.
.lhs_unit <- function(n, d) {
  m <- matrix(0, n, d)
  for (j in seq_len(d)) m[, j] <- (sample.int(n) - stats::runif(n)) / n * 2 - 1
  m
}

# The surface, from a centre in the compact layout.  Returns NULL when the
# deadline passes before every design point is drawn.
.survival_surface <- function(centre, model_bin, link_int, crown_age, rho, radius,
                              lower = NULL, upper = NULL, se = 0.02, n_clades = NULL,
                              max_clades = 40000L, max_lin = 1e4, seed, deadline = Inf) {
  t0 <- .now()
  d <- length(centre)
  p <- (d + 1L) * (d + 2L) / 2L
  npts <- as.integer(2L * p + 1L)
  clades <- 0
  # the clades per point: a provisional P at the centre (other seeds), then
  # the number that gives `se` there
  if (is.null(n_clades)) {
    pre <- .surv_direct(.expand_pars(centre, model_bin), model_bin, link_int, crown_age,
                        rho, max_lin, seed + 200000L, n = 2000L, deadline = deadline)
    if (pre[["clades"]] < 2000) return(NULL)
    clades <- clades + pre[["clades"]]
    n_clades <- .surv_n_for(pre[["P"]], pre[["s2"]], se, lo = 500L, hi = max_clades)
  }
  n_clades <- as.integer(n_clades)
  U <- .with_seed(seed + 7L, rbind(rep(0, d), .lhs_unit(npts - 1L, d)))
  X <- sweep(sweep(U, 2, radius, "*"), 2, centre, "+")
  if (!is.null(lower)) X <- pmax(X, matrix(lower, npts, d, byrow = TRUE))
  if (!is.null(upper)) X <- pmin(X, matrix(upper, npts, d, byrow = TRUE))
  Z <- sweep(sweep(X, 2, centre, "-"), 2, radius, "/")
  r <- matrix(NA_real_, npts, 3L, dimnames = list(NULL, c("logP", "se", "P")))
  for (i in seq_len(npts)) {
    if (.now() > deadline) return(NULL)
    v <- .surv_values(.expand_pars(X[i, ], model_bin), model_bin, link_int, crown_age,
                      rho, max_lin, seed, seq_len(n_clades), deadline)
    if (length(v) < n_clades) return(NULL)
    r[i, ] <- .surv_summary(v)[c("logP", "se", "P")]
    clades <- clades + n_clades
  }
  B <- t(apply(Z, 1, .qbasis))
  w <- 1 / pmax(r[, "se"], 1e-3)^2
  cf <- stats::lm.wfit(B, r[, "logP"], w)$coefficients
  cf[!is.finite(cf)] <- 0
  nm <- .par_names(model_bin)
  colnames(X) <- nm
  structure(list(centre = stats::setNames(as.numeric(centre), nm),
                 radius = stats::setNames(as.numeric(radius), nm),
                 coef = as.numeric(cf), design = cbind(as.data.frame(X), r),
                 model = model_bin, link = link_int, crown_age = crown_age, rho = rho,
                 max_lin = max_lin, se = se, n_clades = n_clades, points = npts,
                 clades = clades, seconds = .now() - t0, seed = seed),
            class = "emphasis_surface")
}

#' A local surface of the survival probability
#'
#' Fits a quadratic surface of \eqn{\log P_\theta(S)} (see
#' \code{\link{survival_mc}}) in a box around \code{pars}, from direct
#' simulation at \eqn{2p + 1} points, \eqn{p = (d+1)(d+2)/2} the
#' coefficients of a full quadratic in the \eqn{d} parameters (31 points for
#' four parameters, 57 for six): the centre and a Latin hypercube of the box
#' of half-width \code{radius}. Every point draws the same clades (common
#' random numbers), as many as give a standard error of \code{se} at the
#' centre, and the quadratic is fitted by weighted least squares, each point
#' weighted by the inverse square of its standard error. Outside the box the
#' surface continues as the tangent plane at the point where the ray from the
#' centre leaves the box; it is capped at \eqn{\log P = 0}.
#'
#' The surface is the conditioning \code{estimate_rates(cond = "surface")}
#' builds around its iterate, and it can be passed to \code{cond} itself, or
#' to \code{\link{predict_survival}}. \code{predict()} returns
#' \eqn{\log P}, \eqn{P} or the gradient of \eqn{\log P} in the parameters.
#'
#' @inheritParams survival_mc
#' @param pars The centre, in the compact layout of \code{model}.
#' @param radius Half-width of the box per parameter (recycled), or
#'   \code{NULL} (default) for \code{0.03} of the span of
#'   \code{lower_bound}--\code{upper_bound}.
#' @param lower_bound,upper_bound Optional parameter box; design points are
#'   clamped into it.
#' @param se Standard error of \eqn{\log P} at the centre the clade count is
#'   set for. Default \code{0.02}.
#' @param n_clades Clades per design point, or \code{NULL} (default) to set
#'   it from \code{se}: 2000 clades (other seeds) give a provisional \eqn{P}
#'   at the centre, and the count is the one that reaches \code{se} there,
#'   between 500 and \code{max_clades}.
#' @param max_time Seconds after which the construction stops with an
#'   error, or \code{NULL} (default) for no limit.
#' @return An object of class \code{"emphasis_surface"}: \code{centre},
#'   \code{radius}, the quadratic's \code{coef} (in the coordinates
#'   \eqn{z = (\theta - \mathrm{centre})/\mathrm{radius}}), the
#'   \code{design} (the points with their \code{logP}, \code{se} and
#'   \code{P}), \code{n_clades} per point, \code{points}, \code{clades} in
#'   all, \code{seconds}, and the model, link, crown age, \code{rho},
#'   lineage cap and seed it was built for.
#' @seealso \code{\link{survival_mc}}, \code{\link{estimate_rates}}
#' @examples
#' \donttest{
#' S <- survival_surface(log(c(0.5, 0.2)), model = "cr", crown_age = 3,
#'                       radius = 0.1, n_clades = 1000, seed = 1)
#' predict(S, log(c(0.52, 0.2)))                   # log P
#' predict(S, log(c(0.52, 0.2)), type = "gradient")
#' }
#' @export
survival_surface <- function(pars, model = "cr", crown_age, rho = 1, link = "exponential",
                             radius = NULL, lower_bound = NULL, upper_bound = NULL,
                             se = 0.02, n_clades = NULL, max_clades = 40000L,
                             max_lin = 1e4, seed = NULL, max_time = NULL) {
  model_bin <- .resolve_model(model)
  link_int  <- .resolve_link(link)
  .surv_check_common(crown_age, rho, max_lin)
  expected_n <- 2L + 2L * sum(model_bin != 0L)
  if (!is.numeric(pars) || length(pars) != expected_n) stop(.pars_error_msg(model_bin, expected_n))
  if (is.null(radius)) {
    if (is.null(lower_bound) || is.null(upper_bound))
      stop("give 'radius', or 'lower_bound' and 'upper_bound' to take 0.03 of their span.", call. = FALSE)
    radius <- 0.03 * (upper_bound - lower_bound)
  }
  radius <- rep_len(as.numeric(radius), expected_n)
  if (any(!is.finite(radius) | radius <= 0)) stop("'radius' must be positive.", call. = FALSE)
  seed <- .surv_seed(seed)
  deadline <- if (is.null(max_time)) Inf else .now() + max_time
  S <- .survival_surface(as.numeric(pars), model_bin, link_int, crown_age, rho, radius,
                         lower = lower_bound, upper = upper_bound, se = se,
                         n_clades = n_clades, max_clades = as.integer(max_clades),
                         max_lin = as.integer(max_lin), seed = seed, deadline = deadline)
  if (is.null(S)) stop("max_time ran out before the surface's design points were drawn.", call. = FALSE)
  S
}

# log P of the surface at one point (compact layout), and its gradient.
.surface_logP <- function(S, theta) {
  z <- (as.numeric(theta) - S$centre) / S$radius
  m <- max(abs(z))
  v <- if (m <= 1) sum(S$coef * .qbasis(z)) else {
    zb <- z / m
    sum(S$coef * .qbasis(zb)) + sum(.qgrad_z(S$coef, zb) * (z - zb))
  }
  min(0, v)
}
.surface_grad <- function(S, theta) {
  z <- (as.numeric(theta) - S$centre) / S$radius
  m <- max(abs(z))
  if (.surface_logP(S, theta) >= 0) return(stats::setNames(rep(0, length(z)), names(S$centre)))
  zb <- if (m <= 1) z else z / m
  stats::setNames(.qgrad_z(S$coef, zb) / S$radius, names(S$centre))
}
# Whether theta lies in the inner `frac` of the surface's box.
.surface_inside <- function(S, theta, frac = 0.5)
  max(abs((as.numeric(theta) - S$centre) / S$radius)) <= frac

# Rows of `newdata` as compact parameter vectors in the surface's order.
.surface_rows <- function(S, newdata) {
  nm <- names(S$centre)
  if (is.null(dim(newdata))) newdata <- matrix(newdata, nrow = 1L,
                                               dimnames = list(NULL, names(newdata)))
  nd <- as.data.frame(newdata)
  if (all(nm %in% names(nd))) nd <- nd[, nm, drop = FALSE]
  else if (ncol(nd) != length(nm))
    stop(sprintf("the surface has %d parameters (%s); newdata has %d columns.",
                 length(nm), paste(nm, collapse = ", "), ncol(nd)), call. = FALSE)
  as.matrix(nd)
}

#' Evaluate a survival surface
#'
#' @param object An \code{"emphasis_surface"} from \code{\link{survival_surface}}
#'   or from a fit's \code{conditioning$surface}.
#' @param newdata A parameter vector, or a matrix or data frame of them (one
#'   row a point; columns matched by name when named as the parameters,
#'   otherwise by position).
#' @param type \code{"logP"} (default), \code{"P"}, or \code{"gradient"},
#'   the gradient of \eqn{\log P} in the parameters (outside the surface's
#'   box, the slope of its tangent plane; zero where the cap at
#'   \eqn{\log P = 0} binds).
#' @param ... Ignored.
#' @return A numeric vector (one value per point), or for \code{"gradient"}
#'   a matrix with one row per point.
#' @export
predict.emphasis_surface <- function(object, newdata, type = c("logP", "P", "gradient"), ...) {
  type <- match.arg(type)
  X <- .surface_rows(object, newdata)
  if (type == "gradient")
    return(t(apply(X, 1, function(th) .surface_grad(object, th))))
  lp <- apply(X, 1, function(th) .surface_logP(object, th))
  if (type == "P") exp(lp) else lp
}

#' Print a survival surface
#' @param x An \code{"emphasis_surface"}.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
print.emphasis_surface <- function(x, ...) {
  cat(sprintf("survival surface (model = %s, crown age %g, rho = %g)\n",
              .model_label(x$model), x$crown_age, x$rho))
  cat("centre:\n"); print(round(x$centre, 6))
  cat("radius:\n"); print(signif(x$radius, 3))
  cat(sprintf("log P at the centre: %.4f (design value %.4f, se %.4f)\n",
              .surface_logP(x, x$centre), x$design$logP[1L], x$design$se[1L]))
  cat(sprintf("%d points x %d clades (%.0f clades in all, %.1f s)\n",
              x$points, x$n_clades, x$clades, x$seconds))
  invisible(x)
}

# --------------------------------------------------------------------------- #
#  estimate_rates(cond = "surface")                                            #
# --------------------------------------------------------------------------- #

# The settings of the surface conditioning, from the control list.
.surface_settings <- function(ctrl, n_tips) {
  s <- list(radius     = ctrl$surface_radius %||% 0.03,
            se         = ctrl$surface_se %||% 0.02,
            final_se   = ctrl$surface_final_se %||% 0.01,
            segment    = as.integer(ctrl$surface_segment %||% 25L),
            max_lin    = as.integer(ctrl$surface_max_lin %||% max(1000L, 20L * n_tips)),
            edge       = ctrl$surface_edge %||% 1e-3,
            restart    = ctrl$surface_restart %||% TRUE,
            max_clades = as.integer(ctrl$surface_max_clades %||% 40000L),
            final_max_clades = as.integer(ctrl$surface_final_max_clades %||% 100000L),
            seed       = ctrl$surface_seed)
  for (nm in c("radius", "se", "final_se", "edge"))
    if (!is.numeric(s[[nm]]) || length(s[[nm]]) != 1L || !(s[[nm]] > 0))
      stop(sprintf("control$surface_%s must be a positive number.", nm), call. = FALSE)
  if (!(s$segment >= 1L)) stop("control$surface_segment must be at least 1.", call. = FALSE)
  if (!is.logical(s$restart) || length(s$restart) != 1L || is.na(s$restart))
    stop("control$surface_restart must be TRUE or FALSE.", call. = FALSE)
  s
}

# Steps (i)-(iv) of the surface conditioning (README, "Conditioning on
# survival"; the experiments repository's E46).  `ctrl` is the merged control,
# with the sampler already resolved; `init_pars` the compact start of (i),
# which is also the nested start of (iv).
.fit_surface_cond <- function(tree, model, model_bin, link, link_int, init_pars,
                              control, ctrl, brts) {
  t_start <- .now()
  deadline <- if (is.null(ctrl$max_time)) Inf else t_start + ctrl$max_time
  left <- function() deadline - .now()
  lb <- as.numeric(ctrl$lower_bound); ub <- as.numeric(ctrl$upper_bound)
  span <- ub - lb
  nm <- .par_names(model_bin)
  crown <- max(brts)
  n_tips <- length(brts) + 1L
  rho <- ctrl$rho
  sc <- .surface_settings(ctrl, n_tips)
  sbase <- .surv_seed(sc$seed)
  max_iter <- as.integer(ctrl$max_iter)
  draws0 <- as.integer(ctrl$sample_size)
  verbose <- isTRUE(ctrl$verbose)
  quiet <- function(expr) if (verbose) expr else suppressWarnings(suppressMessages(expr))
  clamp <- function(th) pmin(pmax(th, lb + 1e-3 * span), ub - 1e-3 * span)
  mc <- list(surfaces = 0L, clades = 0, seconds = 0)

  # what every inner fit shares: the user's control without the surface
  # settings, the box, the sampler that was resolved, and a fixed draw cap so
  # that carrying the grown draw count forward does not move it
  base <- control[!grepl("^surface_", names(control))]
  base[c("num_trees", "sample_size")] <- NULL
  base$lower_bound <- lb; base$upper_bound <- ub
  base$sampling <- ctrl$sampling
  base$max_draws <- ctrl$max_draws %||% as.integer(8L * draws0)
  fit_one <- function(init, cond, iters, draws) {
    cc <- base
    cc$max_iter <- as.integer(iters)
    cc$num_trees <- as.integer(draws); cc$sample_size <- as.integer(draws)
    cc["max_time"] <- list(if (is.finite(deadline)) max(1, left()) else NULL)
    err <- NULL
    f <- tryCatch(quiet(estimate_rates(tree, method = "mcem", model = model, init_pars = init,
                                       control = cc, link = link, cond = cond)),
                  error = function(e) { err <<- conditionMessage(e); NULL })
    list(f = f, err = err)
  }
  draws_of <- function(f) {
    d <- tryCatch(as.integer(utils::tail(f$details$mcem$sample_size, 1L)), error = function(e) NA_integer_)
    if (length(d) == 1L && is.finite(d)) d else NA_integer_
  }
  # the unconditioned log-likelihood of a fit: the drivers report it, and
  # estimate_rates() subtracts the conditioning from it
  ll_uncond <- function(f, S) {
    v <- f$details[["loglik"]]
    if (is.numeric(v) && length(v) == 1L) return(as.numeric(v))
    if (is.null(S)) f$loglik else f$loglik + .surface_logP(S, f$pars)
  }
  direct <- function(th, seed) {
    # the final E-step's convention: the estimate at the returned point gets
    # max_time again
    dl <- if (is.null(ctrl$max_time)) Inf else .now() + ctrl$max_time
    d <- .surv_direct(.expand_pars(as.numeric(th), model_bin), model_bin, link_int, crown,
                      rho, sc$max_lin, seed, n = 2000L, sigma = sc$final_se,
                      hi = sc$final_max_clades, deadline = dl)
    mc$clades <<- mc$clades + d[["clades"]]; mc$seconds <<- mc$seconds + d[["seconds"]]
    d
  }

  # (iii): the conditioned continuation from `start`, in segments, the surface
  # rebuilt whenever the iterate leaves the inner half of its box
  seg_run <- function(start, sbase) {
    th <- clamp(start); S <- NULL; nsurf <- 0L; iters <- 0L; draws <- draws0
    nseg <- 0L; last <- NULL; last_S <- NULL; stop <- NA_character_; err <- NULL; trace <- list()
    repeat {
      if (left() <= 0) { stop <- "time_budget"; break }
      if (iters >= max_iter) { stop <- "max_iter"; break }
      if (is.null(S)) {
        S <- .survival_surface(th, model_bin, link_int, crown, rho, sc$radius * span,
                               lower = lb, upper = ub, se = sc$se,
                               max_clades = sc$max_clades, max_lin = sc$max_lin,
                               seed = .surv_seed(sbase + 1009 * nsurf), deadline = deadline)
        if (is.null(S)) { stop <- "time_budget"; break }
        nsurf <- nsurf + 1L
        mc$surfaces <<- mc$surfaces + 1L
        mc$clades <<- mc$clades + S$clades; mc$seconds <<- mc$seconds + S$seconds
        if (verbose) .vcat(sprintf("[surface] surface %d at %s: %d points x %d clades, %.1f s\n",
                                   nsurf, paste(round(th, 4), collapse = " "), S$points,
                                   S$n_clades, S$seconds))
      }
      th <- clamp(th)
      r <- fit_one(th, S, min(sc$segment, max_iter - iters), draws)
      nseg <- nseg + 1L
      if (is.null(r$f)) { stop <- "error"; err <- r$err; break }
      last <- r$f; last_S <- S
      tr <- last$details$mcem
      if (!is.null(tr) && nrow(tr)) { tr <- as.data.frame(tr); tr$segment <- nseg; trace[[nseg]] <- tr }
      iters <- iters + (if (is.finite(last$iterations)) last$iterations else 0L)
      d2 <- draws_of(last); if (is.finite(d2)) draws <- d2
      th <- as.numeric(last$pars)
      stop <- last$stop_reason
      inside <- .surface_inside(S, th)
      if (identical(stop, "converged") && inside) break
      if (stop %in% c("time_budget", "e_step_failure", "m_step_failure")) break
      if (iters >= max_iter) { stop <- "max_iter"; break }
      if (!inside) S <- NULL
    }
    out <- list(fit = last, surface = last_S, pars = th, stop_reason = stop, error = err,
                iterations = iters, segments = nseg, surfaces = nsurf, draws = draws,
                trace = if (length(trace)) do.call(rbind, trace) else NULL)
    if (!is.null(last) && all(is.finite(last$pars))) {
      d <- direct(last$pars, .surv_seed(sbase + 90000))
      out$loglik_uncond <- ll_uncond(last, last_S)
      out$logP <- d[["logP"]]; out$logP_se <- d[["se"]]; out$logP_clades <- d[["clades"]]
      out$loglik <- out$loglik_uncond - out$logP
    } else out$loglik <- NA_real_
    out
  }

  # (i) the unconditioned fit
  u <- fit_one(init_pars, NULL, max_iter, draws0)
  if (is.null(u$f)) stop(u$err, call. = FALSE)
  U <- u$f
  if (verbose) .vcat(sprintf("[surface] unconditioned fit: %s, %s iterations\n",
                             U$stop_reason, U$iterations))
  runs <- list(unconditioned = list(pars = U$pars, loglik = U$loglik, stop_reason = U$stop_reason,
                                    iterations = U$iterations))
  # (ii)-(iii) from the unconditioned estimate
  C <- seg_run(as.numeric(U$pars), sbase)
  runs$from_unconditioned <- C[setdiff(names(C), c("fit", "surface"))]
  win <- C; start <- "unconditioned"
  # (iv) from the nested start when the unconditioned gamma_0 is on an edge
  g0 <- match("gamma_0", nm)
  on_edge <- abs(U$pars[[g0]] - lb[g0]) <= sc$edge * span[g0] ||
             abs(U$pars[[g0]] - ub[g0]) <= sc$edge * span[g0]
  restart <- isTRUE(sc$restart) && isTRUE(on_edge) && left() > 0
  if (restart) {
    N <- seg_run(as.numeric(init_pars), sbase + 25000L)
    runs$from_nested <- N[setdiff(names(N), c("fit", "surface"))]
    if (is.finite(N$loglik) && (!is.finite(win$loglik) || N$loglik > win$loglik)) {
      win <- N; start <- "nested"
    }
  }

  f <- win$fit
  if (is.null(f)) {
    # no conditioned segment finished: the unconditioned fit, conditioned at
    # its estimate by a direct value
    f <- U
    d <- direct(U$pars, .surv_seed(sbase + 95000))
    win$logP <- d[["logP"]]; win$logP_se <- d[["se"]]; win$logP_clades <- d[["clades"]]
    win$loglik <- ll_uncond(U, NULL) - win$logP
    win$stop_reason <- if (identical(win$stop_reason, "error")) win$stop_reason else "time_budget"
    start <- "none"
  }
  f$loglik <- win$loglik
  if (isTRUE(is.finite(f$loglik_var)) && isTRUE(is.finite(win$logP_se))) f$loglik_var <- f$loglik_var + win$logP_se^2
  f$AIC <- if (is.finite(f$loglik)) -2 * f$loglik + 2 * f$n_pars else NA_real_
  f$stop_reason <- win$stop_reason
  f$iterations <- win$iterations %||% f$iterations
  f$cond <- TRUE
  f$conditioning <- list(
    method = "surface", start = start, restart = restart,
    surfaces = mc$surfaces, clades = mc$clades, seconds = mc$seconds,
    logP = win$logP, logP_se = win$logP_se, logP_clades = win$logP_clades,
    loglik_uncond = win$loglik_uncond %||% NA_real_,
    surface = win$surface, runs = runs,
    settings = sc[setdiff(names(sc), "seed")], seed = sbase,
    seconds_total = .now() - t_start)
  if (!is.null(win$error))
    warning("the conditioned fit stopped on an error: ", win$error, call. = FALSE)
  f$notes <- .fit_notes(f, n_tips = n_tips, link_int = link_int, emit = TRUE)
  f
}

#' The per-lineage design of an augmented tree
#'
#' One row per lineage alive on each segment between consecutive events of an
#' augmented tree: the segment's bounds, the covariates the rates read there
#' (\code{N}, the clade mean pendant age \code{M}, the lineage's pendant start
#' \code{ts} so that \code{D = (t - ts) - M}, and its \code{ED} at the segment
#' start, growing at slope one within it), and whether the lineage's own event
#' ends the segment (\code{event}: 0 none, 1 speciation, 2 extinction).  These
#' are the numbers the likelihood is built on; \code{\link{covariate_path}}
#' runs a selection path over them.
#'
#' @param tree An augmented-tree data frame with the topology (the frames
#'   \code{augment_trees} returns when given parent ids), or a \code{phylo}
#'   object, which is taken as fully observed.
#' @return A data frame with columns \code{seg}, \code{lineage}, \code{id},
#'   \code{event}, \code{t0}, \code{t1}, \code{N}, \code{M}, \code{ts},
#'   \code{ed0}.
#' @export
lineage_table <- function(tree) {
  if (inherits(tree, "phylo")) tree <- .observed_frame(tree)
  lineage_table_cpp(tree)
}

# A phylo as the estimator's own frame with no augmented lineage: a constant
# rate with no extinction inserts nothing, so the node list is the observed one.
.observed_frame <- function(phy) {
  brts <- .extract_brts(phy)
  a <- augment_trees(as.numeric(brts), c(0.3, 0, 0, 0, 0, 0, 0, 0),
                     1L, 500L, 100L, 1e6, 1L,
                     model = c(0L, 0L, 0L), link = 0L, rho = 1,
                     parent_tip_start = .pts(brts), parent_id = .pid(brts))
  a$trees[[1L]]
}

#' Selection along a path over the covariates
#'
#' Instead of comparing nested models one pair at a time, this traces the
#' differential-geometric LARS path of Augugliaro, Mineo and Wit (2013) over
#' all the covariates a rate can read --- \code{N}, \code{D} and \code{ED} ---
#' on the same augmented trees an M-step maximises over.  Under the
#' exponential link the complete-data log-likelihood is a log-linear point
#' process: one row per lineage and segment between consecutive events, with
#' the rate at the segment's midpoint, the segment's length as exposure, and
#' a count of one when that lineage's own event ends the segment.  The
#' speciation and extinction halves separate, so each rate gets its own path.
#'
#' The path is the curve on which the Rao score statistics of the active
#' covariates are equal in size and the inactive ones smaller,
#' \eqn{|r_j(\beta)| = \gamma} for \eqn{j} active, traced from the
#' \eqn{\gamma} at which the first covariate enters down to zero, with the
#' intercept always free.  It is solved here on a grid of \eqn{\gamma} by
#' Newton steps rather than by \pkg{dglars}'s predictor--corrector, because
#' the point process needs a true exposure offset, which that package does
#' not take; on a constant-exposure problem the two agree
#' (\code{tests/testthat/test-lineage-table.R}).  The order of entry says
#' which covariates the data ask for first; AIC and BIC along the path pick a
#' point, with BIC's sample size the number of events.
#'
#' Importance weights enter the likelihood as row weights, one per tree and
#' normalised to sum to one, so the objective is the expected complete-data
#' log-likelihood of one tree --- the M-step's Q function --- and the sample
#' size behind BIC is the weighted mean number of events per tree.
#'
#' @param trees A list of augmented-tree data frames with the topology, as
#'   \code{augment_trees} returns them with parent ids, or a single
#'   \code{phylo} taken as fully observed.
#' @param log_weights Optional log importance weights, one per tree; without
#'   them every tree weighs the same.
#' @param n_trees How many trees to use: the \code{n_trees} of largest weight
#'   (the first \code{n_trees} without weights), renormalised.
#' @param covariates Which covariates to offer the path; any of \code{"N"},
#'   \code{"D"}, \code{"ED"}, \code{"EDc"} (ED centred on the alive lineages).
#' @param rate \code{"speciation"}, \code{"extinction"} or both.
#' @param criterion \code{"BIC"} or \code{"AIC"} along the path.
#' @param n_gamma Number of points on the path.
#' @return A list with one element per rate: the path (a data frame with one
#'   row per \eqn{\gamma}: the coefficients, the log-likelihood, the active
#'   set size, AIC and BIC), the entry order of the covariates (the
#'   \eqn{\gamma} at which each enters), the coefficients and active set at
#'   the criterion's optimum, and the frame's size.
#' @export
covariate_path <- function(trees, log_weights = NULL, n_trees = 40L,
                           covariates = c("N", "D", "ED"),
                           rate = c("speciation", "extinction"),
                           criterion = c("BIC", "AIC"), n_gamma = 60L) {
  criterion <- match.arg(criterion)
  rate <- match.arg(rate, several.ok = TRUE)
  covariates <- match.arg(covariates, c("N", "D", "ED", "EDc"), several.ok = TRUE)
  if (inherits(trees, "phylo")) trees <- list(.observed_frame(trees))
  if (!is.list(trees) || !length(trees)) stop("'trees' must be a non-empty list of augmented trees")
  if (!is.null(log_weights)) {
    if (length(log_weights) != length(trees)) stop("one log weight per tree")
    w <- exp(log_weights - max(log_weights))
    pick <- order(w, decreasing = TRUE)[seq_len(min(n_trees, length(trees)))]
    w <- w[pick] / sum(w[pick])
  } else {
    pick <- seq_len(min(n_trees, length(trees)))
    w <- rep(1 / length(pick), length(pick))
  }
  parts <- lapply(seq_along(pick), function(k) {
    f <- .segment_rows(lineage_table_cpp(trees[[pick[k]]])); f$w <- w[k]; f })
  frame <- do.call(rbind, parts)
  ess <- 1 / sum(w^2)
  X <- cbind(N = frame$N, D = frame$D, ED = frame$ED, EDc = frame$EDc)[, covariates, drop = FALSE]
  out <- list()
  for (r in rate) {
    y <- if (r == "speciation") frame$y_spec else frame$y_ext
    path <- .dglars_pp(X, y, frame$dt, frame$w, n_gamma = n_gamma)
    crit <- if (criterion == "BIC") path$BIC else path$AIC
    best <- which.min(crit)
    cf <- unlist(path[best, c("(Intercept)", covariates)])
    entry <- vapply(covariates, function(j) {
      on <- which(path[[j]] != 0); if (length(on)) path$gamma[on[1L]] else NA_real_ }, 0)
    out[[r]] <- list(path = path, entry_gamma = entry,
                     entry_order = names(sort(entry[!is.na(entry)], decreasing = TRUE)),
                     step = best, criterion = criterion, coef = cf,
                     active = covariates[cf[covariates] != 0],
                     n_rows = nrow(frame), n_events = sum(frame$w * y), n_trees = length(pick),
                     ess = ess)
  }
  out
}

# The Poisson rows of one tree from its lineage table: one per lineage and
# segment, covariates at the segment's midpoint, its length, and the counts of
# the lineage's own events at its end.
.segment_rows <- function(tab) {
  if (!nrow(tab)) return(NULL)
  tm <- (tab$t0 + tab$t1) / 2
  ED <- tab$ed0 + (tm - tab$t0)
  data.frame(N = tab$N,
             D = (tm - tab$ts) - tab$M,
             ED = ED,
             EDc = ED - stats::ave(ED, tab$seg, FUN = mean),
             dt = tab$t1 - tab$t0,
             y_spec = as.integer(tab$event == 1L),
             y_ext = as.integer(tab$event == 2L))
}

# The dgLARS path for a log-linear point process with exposure and row weights.
#
# Rows i with covariates x_i, exposure e_i, count y_i and weight w_i;
#   l(b0, b) = sum_i w_i [ y_i (b0 + x_i b) - e_i exp(b0 + x_i b) ].
# Rao score statistic of covariate j at (b0, b):
#   r_j = U_j / sqrt(I_jj),  U_j = sum_i w_i x_ij (y_i - mu_i),  I_jj = sum_i w_i x_ij^2 mu_i.
# The path: for gamma from the first entry down to 0, the active set A and
# the point solving U_0 = 0 and r_j = s_j gamma (j in A), s_j the sign at
# entry; an inactive covariate enters when |r_j| reaches gamma.  Solved on a
# geometric grid of gamma by Newton on the (|A| + 1)-dimensional system with a
# finite-difference Jacobian, warm-started from the previous gamma.
.dglars_pp <- function(X, y, e, w = rep(1, length(y)), n_gamma = 60L, gamma_min_frac = 1e-3) {
  X <- as.matrix(X); p <- ncol(X); n <- nrow(X); ne <- sum(w * y); we <- sum(w * e)
  nm <- colnames(X)
  loglik <- function(b0, b) { eta <- b0 + drop(X %*% b); sum(w * (y * eta - e * exp(eta))) }
  scores <- function(b0, b) {
    mu <- e * exp(b0 + drop(X %*% b))
    U <- drop(crossprod(X, w * (y - mu))); I <- drop(crossprod(X^2, w * mu))
    list(U0 = sum(w * (y - mu)), r = U / sqrt(pmax(I, 1e-300)), mu = mu)
  }
  # the null point: intercept only
  b0 <- log(ne / we); b <- numeric(p); names(b) <- nm
  sc <- scores(b0, b)
  gamma_max <- max(abs(sc$r))
  gammas <- exp(seq(log(gamma_max), log(gamma_max * gamma_min_frac), length.out = n_gamma))
  active <- logical(p); sgn <- numeric(p)
  rows <- vector("list", n_gamma)
  for (k in seq_len(n_gamma)) {
    g <- gammas[k]
    # admit every inactive covariate whose score has reached the current level
    repeat {
      sc <- scores(b0, b)
      cand <- which(!active & abs(sc$r) >= g * (1 - 1e-8))
      if (!length(cand)) break
      j <- cand[which.max(abs(sc$r[cand]))]
      active[j] <- TRUE; sgn[j] <- sign(sc$r[j])
      # solve the system at this gamma with the enlarged set
      sol <- .pp_solve(b0, b, active, sgn, g, scores)
      b0 <- sol$b0; b <- sol$b
    }
    if (any(active)) { sol <- .pp_solve(b0, b, active, sgn, g, scores); b0 <- sol$b0; b <- sol$b }
    else { b0 <- log(ne / we) }
    ll <- loglik(b0, b); df <- 1L + sum(b != 0)     # nonzero coefficients: the first row is the null
    rows[[k]] <- c(gamma = g, "(Intercept)" = b0, b, loglik = ll, df = df,
                   AIC = -2 * ll + 2 * df, BIC = -2 * ll + log(max(ne, 2)) * df)
  }
  out <- as.data.frame(do.call(rbind, rows))
  names(out) <- c("gamma", "(Intercept)", nm, "loglik", "df", "AIC", "BIC")
  out
}

# Newton on F(theta) = (U_0, r_j - s_j gamma for j in A) = 0, theta = (b0, b_A),
# with a finite-difference Jacobian; a step is halved while it does not reduce
# |F|, and the solve stops at 1e-9 relative.
.pp_solve <- function(b0, b, active, sgn, g, scores, max_iter = 50L) {
  A <- which(active); q <- length(A)
  theta <- c(b0, b[A])
  Fn <- function(th) { bb <- b; bb[A] <- th[-1L]; sc <- scores(th[1L], bb)
    c(sc$U0, sc$r[A] - sgn[A] * g) }
  Fv <- Fn(theta)
  for (it in seq_len(max_iter)) {
    if (max(abs(Fv)) < 1e-9 * max(1, g)) break
    J <- matrix(0, q + 1L, q + 1L)
    for (m in seq_len(q + 1L)) {
      h <- 1e-6 * max(1, abs(theta[m])); th2 <- theta; th2[m] <- th2[m] + h
      J[, m] <- (Fn(th2) - Fv) / h
    }
    step <- tryCatch(solve(J, -Fv), error = function(err) -Fv * 1e-3)
    lam <- 1
    repeat {
      th_new <- theta + lam * step; F_new <- Fn(th_new)
      if (all(is.finite(F_new)) && (max(abs(F_new)) < max(abs(Fv)) || lam < 1e-4)) break
      lam <- lam / 2
    }
    theta <- th_new; Fv <- F_new
  }
  bb <- b; bb[A] <- theta[-1L]
  list(b0 = theta[1L], b = bb)
}

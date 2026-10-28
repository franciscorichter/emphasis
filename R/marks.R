# Latent marks (dev/latent-marks/design.md): presets that turn a process's
# parameters into the marks Model::loglik_marks reads.  Internal while the
# extension is developed on the latent-marks branch.

# The latent-rate preset (ClaDS2): a lineage's log rate offset u lives on a grid
# of S points over [-span, span]; at a split each daughter draws its offset as
# u + N(log_alpha, sigma^2), independently, discretised on the grid (the normal
# density at the grid points, renormalised per row); speciation and extinction
# are multiplied alike (constant turnover within a lineage).  The root offset
# is 0, the grid's centre, since the intercepts carry the level.  The span must
# cover the spread the walk reaches after the largest K in a tree,
# sigma * sqrt(K); the default is three times that for K = kmax.
.latent_rate_marks <- function(sigma, log_alpha = 0, S = NULL, span = NULL, kmax = 20, smax = 101L) {
  if (!is.finite(sigma) || sigma < 0) stop("sigma must be finite and non-negative")
  if (sigma == 0 && log_alpha == 0) {
    return(list(S = 1L, u = 0, log_clam = 0, log_cmu = 0, C = 1, P = 1, pi = 1))
  }
  if (is.null(span)) span <- max(3 * max(sigma, 1e-3) * sqrt(kmax), 3 * abs(log_alpha) * kmax)
  # three grid points per sigma across the span, odd, at most smax
  if (is.null(S)) S <- min(as.integer(smax), 2L * ceiling(span / (max(sigma, 1e-3) / 3)) + 1L)
  S <- as.integer(S)
  if (S %% 2L != 1L) stop("S must be odd, so that the grid holds the root offset 0")
  u <- seq(-span, span, length.out = S)
  P <- outer(u, u, function(from, to) stats::dnorm(to, from + log_alpha, max(sigma, 1e-12)))
  P <- P / rowSums(P)
  pi <- as.numeric(seq_len(S) == (S + 1L) %/% 2L)
  list(S = S, u = u, log_clam = u, log_cmu = u, C = numeric(0), P = P, pi = pi)
}

# eval_logf with the marks of a preset.
.eval_logf_marks <- function(pars, trees, model, marks, rho = 1) {
  eval_logf_marks(pars, trees, model = model, link = 1L, rho = rho, S = marks$S,
                  log_clam = marks$log_clam, log_cmu = marks$log_cmu,
                  C = if (length(marks$C)) marks$C else numeric(0),
                  P = if (marks$S > 1L) as.vector(t(marks$P)) else numeric(0), pi = marks$pi)
}

# The forest and the complete-tree frame of a simulated L-table (crown age T):
# every lineage, extinct ones included, as the node list eval_logf reads.
.ltable_forest <- function(L, T) {
  b <- T - L[, 1]; b[1:2] <- 0
  list(birth = b, death = ifelse(L[, 4] == -1, Inf, T - L[, 4]),
       parent = c(NA, NA, match(L[-(1:2), 2], L[, 3])), n = nrow(L))
}
.ltable_complete_frame <- function(L, T) {
  f <- .ltable_forest(L, T)
  ev <- rbind(data.frame(t = f$birth[-(1:2)], type = rep(1L, f$n - 2), row = seq_len(f$n - 2L) + 2L),
              data.frame(t = f$death[is.finite(f$death)], type = rep(2L, sum(is.finite(f$death))),
                         row = which(is.finite(f$death))))
  ev <- ev[order(ev$t), ]; prow <- f$parent
  ids <- c(-2L, -3L, seq_len(f$n - 2L) - 1L)
  ts <- c(0, 0, rep(NA_real_, f$n - 2L)); alive <- c(TRUE, TRUE, rep(FALSE, f$n - 2L))
  out <- vector("list", nrow(ev) + 1L)
  for (k in seq_len(nrow(ev))) {
    t <- ev$t[k]; r <- ev$row[k]; a <- which(alive); n <- length(a); pd <- sum(t - ts[a])
    if (ev$type[k] == 1L) {
      p <- prow[r]
      out[[k]] <- data.frame(brts = t, n = n, t_ext = if (is.finite(f$death[r])) f$death[r] else 10e10,
                             pd = pd, tip_start = t, focal_tip_start = ts[p], clade = 0L,
                             id = ids[r], parent_id = ids[p])
      ts[p] <- t; ts[r] <- t; alive[r] <- TRUE
    } else {
      out[[k]] <- data.frame(brts = t, n = n, t_ext = 0, pd = pd, tip_start = ts[r],
                             focal_tip_start = ts[r], clade = 0L, id = ids[r], parent_id = -1L)
      alive[r] <- FALSE
    }
  }
  a <- which(alive)
  out[[nrow(ev) + 1L]] <- data.frame(brts = T, n = length(a), t_ext = 10e10, pd = sum(T - ts[a]),
                                     tip_start = T, focal_tip_start = -1, clade = 0L,
                                     id = -1L, parent_id = -1L)
  do.call(rbind, out)
}

# One complete tree simulated with the marks of a preset: the L-table, each
# row's final mark, the status, and the complete-tree frame.  `pars` is the
# compact vector of `model`, as simulate_tree() takes it.
.simulate_marks <- function(pars, model, max_t, marks, max_N = 1e4, max_tries = 100L, seed = 0L) {
  mb <- .resolve_model(model)
  pe <- .expand_pars(pars, mb)
  for (i in seq_len(max_tries)) {
    r <- simulate_div_tree_cpp(pe, mb, max_t, as.integer(max_N), 0L, 1L,
                               seed = if (seed > 0L) as.integer(seed + i - 1L) else .draw_seed(),
                               marks_S = marks$S, log_clam = marks$log_clam, log_cmu = marks$log_cmu,
                               marks_C = if (length(marks$C)) marks$C else numeric(0),
                               marks_P = if (marks$S > 1L) as.vector(t(marks$P)) else numeric(0),
                               marks_pi = marks$pi)
    if (r$status != "extinct") break
  }
  r$frame <- if (r$status == "done") .ltable_complete_frame(r$Ltable, max_t) else NULL
  r
}

# ---- MCEM with marks: an R-level prototype (phase 3) ----
#
# The E-step is the package's: augment_trees() draws hidden lineages from the
# proposal at the current base rates (the marks ignored).  Each augmented tree
# is weighed by its mark-marginal likelihood over the proposal density,
# log w = loglik_marks(theta) - log g, which is what makes the marks enter
# without any change to the sampler.  The M-step maximises the weighted
# mark-marginal log-likelihood over the base parameters and the mark
# parameters with optim().  `marks_of(theta_m)` turns the mark parameters
# (on the optimiser's scale) into a preset.
.mcem_marks <- function(phy, model, start, start_m, marks_of, n_iter = 20L, N = 200L,
                        N_growth = 1.15, N_max = 2000L, rho = 1, max_missing = 5000L, seed = 1L,
                        verbose = FALSE) {
  mb <- .resolve_model(model)
  brts <- .extract_brts(phy)
  np <- length(start)
  th <- c(start, start_m)
  trace <- vector("list", n_iter)
  draw <- function(th, n, s) {
    a <- augment_trees(brts, .expand_pars(th[seq_len(np)], mb), sample_size = as.integer(n),
                       maxN = as.integer(50 * n), max_missing = as.integer(max_missing), max_lambda = 1e6,
                       num_threads = 1L, model = mb, link = 1L, rho = rho,
                       parent_tip_start = .pts(brts), seed = as.integer(s), parent_id = .pid(brts))
    a
  }
  wlogf <- function(th, trees) .eval_logf_marks(.expand_pars(th[seq_len(np)], mb), trees, mb,
                                                 marks_of(th[-seq_len(np)]), rho = rho)
  for (it in seq_len(n_iter)) {
    a  <- draw(th, N, seed + it)
    lw <- wlogf(th, a$trees) - a$logg
    ok <- is.finite(lw)
    m  <- max(lw[ok]); w <- exp(lw[ok] - m); w <- w / sum(w)
    ess <- 1 / sum(w^2)
    trees <- a$trees[ok]
    obj <- function(x) { v <- wlogf(x, trees); if (any(!is.finite(v))) return(1e10); -sum(w * v) }
    f <- stats::optim(th, obj, control = list(maxit = 200 + 50 * length(th)))
    th <- f$par
    trace[[it]] <- c(iter = it, N = N, ess = ess, loglik = m + log(mean(exp(lw[ok] - m))), th)
    if (verbose) message(sprintf("iter %d  N %d  ESS %.1f  logL %.3f  %s", it, N, ess,
                                 trace[[it]][["loglik"]], paste(round(th, 3), collapse = " ")))
    N <- min(as.integer(ceiling(N * N_growth)), as.integer(N_max))
  }
  list(estimate = th, trace = do.call(rbind, trace))
}

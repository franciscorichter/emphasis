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
.latent_rate_marks <- function(sigma, log_alpha = 0, S = 21L, span = NULL, kmax = 20) {
  S <- as.integer(S)
  if (S %% 2L != 1L) stop("S must be odd, so that the grid holds the root offset 0")
  if (!is.finite(sigma) || sigma < 0) stop("sigma must be finite and non-negative")
  if (sigma == 0 && log_alpha == 0) {
    return(list(S = 1L, u = 0, log_clam = 0, log_cmu = 0, C = 1, P = 1, pi = 1))
  }
  if (is.null(span)) span <- max(3 * max(sigma, 1e-3) * sqrt(kmax), 3 * abs(log_alpha) * kmax)
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

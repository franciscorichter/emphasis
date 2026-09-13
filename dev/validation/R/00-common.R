# ---------------------------------------------------------------------------
# 00-common.R — shared constants, exact references, design tables, job table.
#
# Sourced by every other script.  Contains no side effects other than defining
# objects and creating the output directories.
#
# emphasis is loaded from the library named by, in order:
#   (1) the R option `emphasis.lib`,
#   (2) the environment variable EMPHASIS_LIB,
#   (3) the default .libPaths().
# ---------------------------------------------------------------------------

VAL_ROOT <- local({
  marker <- function(d) file.exists(file.path(d, "R", "00-common.R"))
  cand <- character(0)
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) == 1L) cand <- c(cand, dirname(dirname(normalizePath(f, mustWork = FALSE))))
  cand <- c(cand, file.path(getwd(), "dev", "validation"), getwd(),
            dirname(getwd()),
            "/Users/pancho/Code/emphasis/dev/validation")
  hit <- Find(marker, cand)
  if (is.null(hit)) cand[length(cand)] else hit
})

VAL_DIR <- list(
  data    = file.path(VAL_ROOT, "data"),
  results = file.path(VAL_ROOT, "results"),
  figures = file.path(VAL_ROOT, "figures"),
  tables  = file.path(VAL_ROOT, "tables"),
  scratch = file.path(VAL_ROOT, "scratch")
)
for (d in VAL_DIR) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# --- library resolution ----------------------------------------------------

val_lib <- function() {
  lib <- getOption("emphasis.lib", NULL)
  if (is.null(lib) || !nzchar(lib)) lib <- Sys.getenv("EMPHASIS_LIB", "")
  if (nzchar(lib)) lib else NA_character_
}

val_load_emphasis <- function(lib = val_lib()) {
  if (!is.na(lib) && nzchar(lib)) .libPaths(c(lib, .libPaths()))
  suppressPackageStartupMessages(suppressWarnings(library(emphasis)))
  invisible(TRUE)
}

val_load_refs <- function() {
  suppressPackageStartupMessages(suppressWarnings({
    library(ape); library(DDD); library(TreeSim)
  }))
  invisible(TRUE)
}

# Build fingerprint: every result row carries it; 04 refuses to pool rows with
# different fingerprints, and 03 refuses to resume into a directory with one.
val_build_fingerprint <- function() {
  p <- tryCatch(system.file(package = "emphasis"), error = function(e) "")
  desc <- file.path(p, "DESCRIPTION")
  so   <- list.files(file.path(p, "libs"), pattern = "\\.(so|dylib)$",
                     full.names = TRUE, recursive = TRUE)
  md5 <- function(f) if (length(f) && file.exists(f[1]))
    unname(tools::md5sum(f[1])) else NA_character_
  list(
    version      = tryCatch(as.character(utils::packageVersion("emphasis")),
                            error = function(e) NA_character_),
    lib          = p,
    md5_desc     = md5(desc),
    md5_so       = md5(so),
    r_version    = as.character(getRversion()),
    ddd_version  = tryCatch(as.character(utils::packageVersion("DDD")),
                            error = function(e) NA_character_)
  )
}

val_fingerprint_key <- function(fp) paste(fp$version, fp$md5_desc, fp$md5_so, sep = "|")

# --- silencing DDD ---------------------------------------------------------
# bd_loglik / dd_loglik print unconditionally; every call is wrapped.
.q <- function(expr) { junk <- utils::capture.output(v <- suppressMessages(expr)); v }

# --- exact CR log-likelihood ----------------------------------------------
# Nee et al. (1994) crown likelihood in the per-lineage (labelled-history)
# convention emphasis's f uses.  brts: decreasing, crown age first, length n-1.
# Verified in 00-selfcheck.R to equal
#   DDD::bd_loglik(pars1 = c(l, m, 0, 0), pars2 = c(0, 0, 1, 0, 2), brts, 0)
# to < 1e-6, and to equal emphasis's BDI fhat exactly (constant = 0).
ll_cr_nee <- function(lambda, mu, brts) {
  if (!is.finite(lambda) || !is.finite(mu) || lambda <= 0 || mu < 0) return(-Inf)
  r  <- lambda - mu
  Tc <- brts[1L]
  if (abs(r) < 1e-8 * max(lambda, 1)) {
    # critical case lambda = mu: p1(t) -> 1 / (1 + lambda t)^2
    lp1c <- function(t) -2 * log1p(lambda * t)
    v <- 2 * lp1c(Tc) + sum(log(lambda) + lp1c(brts[-1L]))
    return(if (is.nan(v)) -Inf else v)
  }
  # p1(t) = r^2 exp(-r t) / (lambda - mu exp(-r t))^2.  The denominator is
  # squared, so it is defined for r < 0 (where it is negative for t > 0) and
  # undefined only at its zero, t = log(lambda/mu)/(mu - lambda), which lies
  # at t < 0 whenever mu > lambda.
  # For r < 0, exp(-r t) overflows at moderate |r| t (n = 200, high turnover),
  # so factor it out: log|den| = log(mu) - r t + log|1 - (lambda/mu) e^{r t}|,
  # and the -r t cancels against the -r t in the numerator.
  lp1 <- if (r > 0) function(t) {
    den <- lambda - mu * exp(-r * t)
    if (any(abs(den) < 1e-300)) return(rep(-Inf, length(t)))
    2 * log(abs(r)) - r * t - 2 * log(abs(den))
  } else function(t) {
    if (mu <= 0) return(rep(-Inf, length(t)))
    u <- 1 - (lambda / mu) * exp(r * t)
    if (any(abs(u) < 1e-300)) return(rep(-Inf, length(t)))
    2 * log(abs(r)) + r * t - 2 * log(mu) - 2 * log(abs(u))
  }
  v <- lp1(Tc) * 2 + sum(log(lambda) + lp1(brts[-1L]))
  if (is.nan(v)) -Inf else v
}

ll_cr_ddd <- function(lambda, mu, brts, cond = 0, btorph = 1, soc = 2) {
  # pars2 layout (DDD 5.2.4): c(tdmodel, cond, btorph, verbose, soc)
  .q(DDD::bd_loglik(pars1 = c(lambda, mu, 0, 0),
                    pars2 = c(0, cond, btorph, 0, soc),
                    brts = brts, missnumspec = 0))
}

ll_dd_ddd <- function(lambda, mu, K, brts, lx = NULL, cond = 0, btorph = 1,
                      soc = 2, methode = "analytical") {
  if (is.null(lx)) lx <- val_lx(length(brts) + 1L)
  # pars2 layout (DDD 5.2.4): c(lx, ddmodel, cond, btorph, verbose, soc)
  .q(DDD::dd_loglik(pars1 = c(lambda, mu, K),
                    pars2 = c(lx, 1, cond, btorph, 0, soc),
                    brts = brts, missnumspec = 0, methode = methode))
}

val_lx <- function(n) max(300, 10 * (n + 1))

# --- parameter mapping emphasis <-> DDD ------------------------------------
# dd/linear with gamma_N fixed at 0 is DDD ddmodel = 1:
#   lambda(N) = max(0, beta_0 + beta_N N),  mu = gamma_0
#   beta_0 = lambda0,  beta_N = -(lambda0 - mu0)/K,  gamma_0 = mu0
val_dd_to_emphasis <- function(lambda0, mu0, K)
  c(beta_0 = lambda0, beta_N = -(lambda0 - mu0) / K, gamma_0 = mu0, gamma_N = 0)

val_emphasis_to_dd <- function(pars) {
  b0 <- pars[1L]; bN <- pars[2L]; g0 <- pars[3L]
  if (!is.finite(b0) || !is.finite(bN) || !is.finite(g0))
    return(list(status = "non_finite", lambda0 = NA, mu0 = NA, K = NA))
  if (bN >= 0)                 # CR limit: no diversity dependence
    return(list(status = "cr_limit", lambda0 = b0, mu0 = g0, K = Inf))
  if (b0 <= g0)                # lambda0 <= mu0 has no DDD counterpart
    return(list(status = "outside_ddd", lambda0 = b0, mu0 = g0, K = NA))
  list(status = "ok", lambda0 = b0, mu0 = g0, K = -(b0 - g0) / bN)
}

# Exact log-likelihood at an emphasis estimate, either model.
ll_exact <- function(pars, brts, model = c("cr", "dd"), lx = NULL) {
  model <- match.arg(model)
  if (model == "cr") return(ll_cr_nee(pars[1L], pars[2L], brts))
  m <- val_emphasis_to_dd(pars)
  if (m$status == "non_finite" || m$status == "outside_ddd") return(NA_real_)
  if (m$status == "cr_limit")   return(ll_cr_ddd(m$lambda0, m$mu0, brts))
  ll_dd_ddd(m$lambda0, m$mu0, m$K, brts, lx = lx)
}

# --- exact CR MLE ----------------------------------------------------------
val_cr_mle <- function(brts, starts = NULL, lam_gen = 1, mu_gen = 0) {
  if (is.null(starts))
    starts <- rbind(c(lam_gen, max(mu_gen, 1e-3)),
                    c(2 * lam_gen, 0.25 * lam_gen),
                    c(0.5 * lam_gen, 1e-6))
  # +Inf far from the data would make L-BFGS-B error out, so the objective is
  # capped rather than infinite; the cap is never attained at an optimum.
  nll <- function(th) {
    v <- -ll_cr_nee(th[1L], th[2L], brts)
    if (!is.finite(v)) 1e12 else v
  }
  best <- NULL
  try_start <- function(st, method) {
    st <- pmin(pmax(st, c(1e-6, 0)), c(1e3, 1e3))
    o <- tryCatch(
      if (method == "L-BFGS-B")
        stats::optim(st, nll, method = "L-BFGS-B", lower = c(1e-6, 0),
                     upper = c(1e3, 1e3),
                     control = list(factr = 1e2, pgtol = 1e-10))
      else
        stats::optim(st, nll, method = "Nelder-Mead",
                     control = list(reltol = 1e-12, maxit = 2000)),
      error = function(e) NULL)
    if (!is.null(o) && is.finite(o$value) && o$value < 1e11 &&
        o$par[1] > 0 && o$par[2] >= 0 &&
        (is.null(best) || o$value < best$value)) best <<- o
  }
  for (i in seq_len(nrow(starts))) try_start(starts[i, ], "L-BFGS-B")
  for (i in seq_len(nrow(starts))) try_start(starts[i, ], "Nelder-Mead")
  if (!is.null(best)) try_start(best$par, "L-BFGS-B")   # polish
  if (is.null(best)) return(NULL)
  th <- best$par
  boundary_mu0 <- th[2L] < 1e-8
  # Observed information from the closed form.
  H <- tryCatch(numDeriv::hessian(function(p) ll_cr_nee(p[1L], p[2L], brts), th),
                error = function(e) NULL)
  se <- c(NA_real_, NA_real_)
  if (!is.null(H) && all(is.finite(H))) {
    V <- tryCatch(solve(-H), error = function(e) NULL)
    if (!is.null(V) && all(is.finite(diag(V))) && all(diag(V) > 0))
      se <- sqrt(diag(V))
  }
  if (boundary_mu0 || !is.finite(se[1L])) {
    # profile SE for lambda at mu = mu_hat
    f <- function(l) ll_cr_nee(l, th[2L], brts)
    h2 <- tryCatch(numDeriv::hessian(function(l) f(l), th[1L]),
                   error = function(e) NULL)
    if (!is.null(h2) && is.finite(h2[1, 1]) && h2[1, 1] < 0)
      se[1L] <- sqrt(-1 / h2[1, 1])
  }
  list(pars = c(lambda = unname(th[1L]), mu = unname(th[2L])),
       loglik = -best$value, se = unname(se),
       boundary_mu0 = boundary_mu0,
       mle_r_negative = th[1L] <= th[2L],
       convergence = best$convergence)
}

# Central-difference Hessian with a fixed, known evaluation count
# (2p^2 + 1 calls).  numDeriv's Richardson default costs an order of magnitude
# more, which matters only where one evaluation is a dd_loglik ODE solve.
val_hess_fd <- function(f, x, h = NULL) {
  p <- length(x)
  if (is.null(h)) h <- pmax(abs(x) * 1e-3, 1e-5)
  f0 <- f(x)
  if (!is.finite(f0)) return(NULL)
  H <- matrix(NA_real_, p, p)
  for (i in seq_len(p)) for (j in i:p) {
    ei <- ej <- numeric(p); ei[i] <- h[i]; ej[j] <- h[j]
    v <- (f(x + ei + ej) - f(x + ei - ej) - f(x - ei + ej) + f(x - ei - ej)) /
      (4 * h[i] * h[j])
    H[i, j] <- H[j, i] <- v
  }
  if (!all(is.finite(H))) NULL else H
}

# --- fixed-theta grid ------------------------------------------------------
# MLE plus +/- 1 and +/- 2 SE along each axis, clipped to the feasible region.
val_theta_grid_cr <- function(mle) {
  th <- mle$pars; se <- mle$se
  if (!is.finite(se[1L])) se[1L] <- 0.1 * th[1L]
  if (!is.finite(se[2L])) se[2L] <- max(0.1 * th[1L], 0.05)
  g <- list(c(th[1L], th[2L]))
  for (d in c(-2, -1, 1, 2)) {
    g[[length(g) + 1L]] <- c(max(1e-3, th[1L] + d * se[1L]), th[2L])
    g[[length(g) + 1L]] <- c(th[1L], max(0, th[2L] + d * se[2L]))
  }
  m <- do.call(rbind, g)
  m <- m[!duplicated(round(m, 10)), , drop = FALSE]
  colnames(m) <- c("lambda", "mu")
  m
}

val_theta_grid_dd <- function(mle_pars, se) {
  # mle_pars: compact emphasis (beta_0, beta_N, gamma_0, gamma_N)
  g <- list(mle_pars)
  for (j in c(1L, 2L)) for (d in c(-1, 1)) {
    p <- mle_pars
    p[j] <- p[j] + d * se[j]
    if (j == 1L) p[1L] <- max(1e-3, p[1L])
    if (j == 2L) p[2L] <- min(-1e-8, p[2L])
    g[[length(g) + 1L]] <- p
  }
  m <- do.call(rbind, g)
  colnames(m) <- c("beta_0", "beta_N", "gamma_0", "gamma_N")
  m[!duplicated(round(m, 12)), , drop = FALSE]
}

# ---------------------------------------------------------------------------
# DESIGN TABLES
# ---------------------------------------------------------------------------

# CR cells: lambda_gen = 1 fixes the time unit (the CR likelihood is scale
# invariant; 00-selfcheck.R asserts it).  eps = mu/lambda.
VAL_CR_EPS <- c(0, 0.3, 0.6, 0.9)

val_cr_cells <- function(tier) {
  switch(tier,
    smoke = data.frame(
      n = c(rep(c(20L, 50L), each = 4L), 100L),
      eps = c(rep(VAL_CR_EPS, 2L), 0.6),
      lam = 1, trees = 1L, stringsAsFactors = FALSE),
    main = rbind(
      data.frame(n = rep(c(20L, 50L, 100L), each = 4L),
                 eps = rep(VAL_CR_EPS, 3L), lam = 1,
                 trees = rep(c(20L, 20L, 12L), each = 4L),
                 stringsAsFactors = FALSE),
      # time-scale cell: same process, rates / 10 (probes H1 and the scale
      # dependence of the stopping rule)
      data.frame(n = 50L, eps = 0.5, lam = 0.1, trees = 10L,
                 stringsAsFactors = FALSE)),
    ext = data.frame(n = 200L, eps = c(0.3, 0.9), lam = 1, trees = 5L,
                     stringsAsFactors = FALSE),
    stop("unknown tier: ", tier))
}

# DD regimes.  Trees are kept only when n in [0.5 K, 1.5 K], which keeps K
# identifiable; a regime with K -> infinity makes dd_ML a 100-500 s boundary
# search and a weak reference.
val_dd_cells <- function(tier) {
  base <- data.frame(
    regime  = c("A", "B", "C"),
    lambda0 = c(0.8, 0.8, 0.8),
    mu0     = c(0.1, 0.2, 0.4),
    K       = c(40,  60,  40),
    age     = c(10,  12,  12),
    stress  = c(FALSE, FALSE, TRUE),   # C: high turnover, H58-exposed
    stringsAsFactors = FALSE)
  base$trees <- switch(tier, smoke = 1L, main = 10L, ext = 0L,
                       stop("unknown tier: ", tier))
  base[base$trees > 0L, , drop = FALSE]
}

# DD trees are accepted only inside this band around K (keeps K identifiable)
# and never above VAL_DD_NMAX, above which dd_ML(3 free) costs minutes: the
# reference, not the fit, is what limits the dd arm.
VAL_DD_BAND  <- c(0.5, 1.2)
VAL_DD_NMAX  <- 80L

val_tree_id <- function(kind, ...) paste(c(kind, ...), collapse = "-")

val_seed <- function(id) {
  s <- digest::digest2int(id)
  if (is.na(s) || s == 0L) s <- 1L
  abs(s) %% .Machine$integer.max
}

# ---------------------------------------------------------------------------
# CONFIGURATIONS
# ---------------------------------------------------------------------------
# Every mcem configuration: cond = NULL, rho = 1, link = "linear",
# num_threads = 1, max_missing = 1e4, max_iter = 400, maxN = 20 N,
# hand-supplied bounds containing the exact MLE.
#
#  C1 bdi   N=200 init_far  3 reps      C2 thin  N=200 init_far  3 reps
#  C3 bdi   N=200 init_mle  1 rep       C4 thin  N=200 init_mle  1 rep
#  C5 bdi   N=1000 init_far (subset)    C6 thin  N=1000 init_far (subset)
#  C7 bdi   N=50   init_far (subset)    C8 thin  N=50   init_far (subset)
#  C9 fixed-theta IS grid (subset)      C10 pipeline (subset)
#  C12 box x 10 (n = 50, eps = 0.6)
#  D1..D4, D5/D6, D7 the dd analogues.

val_cr_configs <- function(tier) {
  base <- data.frame(
    config  = c("C1", "C2", "C3", "C4"),
    sampler = c("bdi", "dynamic_fresh", "bdi", "dynamic_fresh"),
    N       = c(200L, 200L, 200L, 200L),
    init    = c("far", "far", "mle", "mle"),
    reps    = c(3L, 3L, 1L, 1L),
    subset  = FALSE, box_scale = 1, stringsAsFactors = FALSE)
  if (tier == "smoke") { base$reps <- 1L; return(base) }
  if (tier == "ext") {
    # n = 200: thinning x 3 replicates plus one BDI fit per tree under a cap.
    # BDI's per-iteration cost at n = 200 (2.5 s) makes replicates unaffordable,
    # so the BDI arm reports cost scaling and whether the cap is reached.
    return(data.frame(
      config  = c("C2", "C1"), sampler = c("dynamic_fresh", "bdi"),
      N = 200L, init = "far", reps = c(3L, 1L), subset = FALSE, box_scale = 1,
      stringsAsFactors = FALSE))
  }
  sub <- data.frame(
    config  = c("C5", "C6", "C7", "C8"),
    sampler = c("bdi", "dynamic_fresh", "bdi", "dynamic_fresh"),
    N       = c(1000L, 1000L, 50L, 50L),
    init    = "far", reps = 1L, subset = TRUE, box_scale = 1,
    stringsAsFactors = FALSE)
  boxarm <- data.frame(
    config = c("C12a", "C12b"), sampler = c("bdi", "dynamic_fresh"),
    N = 200L, init = "far", reps = 2L, subset = TRUE, box_scale = 10,
    stringsAsFactors = FALSE)
  rbind(base, sub, boxarm)
}

val_dd_configs <- function(tier) {
  base <- data.frame(
    config  = c("D1", "D2", "D3", "D4"),
    sampler = c("bdi", "dynamic_fresh", "bdi", "dynamic_fresh"),
    N       = 200L,
    init    = c("far", "far", "mle", "mle"),
    reps    = c(2L, 2L, 1L, 1L),
    subset  = FALSE, box_scale = 1, stringsAsFactors = FALSE)
  if (tier == "smoke") { base$reps <- 1L; return(base) }
  sub <- data.frame(
    config = c("D5", "D6"), sampler = c("bdi", "dynamic_fresh"),
    N = 1000L, init = "far", reps = 1L, subset = TRUE, box_scale = 1,
    stringsAsFactors = FALSE)
  rbind(base, sub)
}

# Subset S: the first trees by index in each cell (4 at n >= 100, 5 otherwise).
val_subset_size <- function(n, tier) {
  if (tier == "smoke") return(1L)
  if (tier == "ext")   return(2L)   # BDI-only fixed-theta check, no pipeline
  if (n >= 100L) 4L else 5L
}

# --- hand-supplied boxes ---------------------------------------------------
val_box_cr <- function(lam_gen, scale = 1)
  list(lower = c(1e-3, 0), upper = c(5 * lam_gen * scale, 5 * lam_gen * scale))

val_box_dd <- function(lambda0, mu0, K)
  list(lower = c(1e-3, -5 * (lambda0 - mu0) / K, 0, 0),
       upper = c(5 * lambda0, 0, 5 * lambda0, 0))

val_init_cr <- function(kind, lam_gen, mle_pars) {
  if (kind == "far") c(2 * lam_gen, 0.5 * lam_gen)
  else c(max(mle_pars[1L], 1e-3), max(mle_pars[2L], 0))
}

val_init_dd <- function(kind, lambda0, mu0, K, mle_pars) {
  if (kind == "far")
    c(1.5 * lambda0, -0.5 * (lambda0 - mu0) / K, 0.25 * lambda0, 0)
  else as.numeric(mle_pars)
}

# --- per-job timeouts ------------------------------------------------------
# Hard kill by the driver; emphasis max_time is set 30 s below it so a clean
# "time_budget" stop is recorded first.
val_timeout <- function(kind, n, N = 200L) {
  base <- switch(kind,
    cr_fit   = if (n <= 20) 120 else if (n <= 50) 300 else if (n <= 100) 600 else 900,
    dd_fit   = if (n <= 50) 300 else 900,
    pipeline = if (n <= 20) 300 else if (n <= 50) 600 else 1500,
    fhat     = if (n <= 50) 300 else 600,
    ref      = 900,
    600)
  if (N >= 1000L) base <- base * 2
  base
}

# --- estimated cost (seconds, single thread) -------------------------------
# Used for longest-first packing and the budget/drop-order check, never for a
# conclusion.  The smoke tier's measured table takes precedence when it exists:
# the pre-fix formula below under-predicted the smoke run by roughly 3x.
val_load_calibration <- function(tier = "smoke") {
  f <- file.path(VAL_DIR$results, tier, "cost-calibration.csv")
  if (file.exists(f)) utils::read.csv(f, stringsAsFactors = FALSE) else NULL
}

val_est_cost_calib <- function(calib, kind, sampler, n, N = 200L) {
  if (is.null(calib)) return(NA_real_)
  sk <- if (length(sampler) != 1L || is.na(sampler)) "none" else sampler
  d <- calib[calib$kind == kind & calib$sampler == sk, , drop = FALSE]
  if (!nrow(d)) return(NA_real_)
  d <- d[order(d$n), ]
  v <- if (nrow(d) == 1L) d$elapsed[1] else
    stats::approx(d$n, d$elapsed, xout = n, rule = 2)$y
  v * (N / 200)
}

val_est_cost <- function(kind, sampler, n, N = 200L, calib = NULL) {
  v <- val_est_cost_calib(calib, kind, sampler, n, N)
  if (is.finite(v)) return(v)
  .val_est_cost_prefix(kind, sampler, n, N)
}

.val_est_cost_prefix <- function(kind, sampler, n, N = 200L) {
  per_it <- switch(sampler %||% "bdi",
    bdi           = stats::approx(c(20, 50, 100, 200), c(0.13, 0.40, 1.54, 2.51),
                                  xout = n, rule = 2)$y,
    dynamic_fresh = stats::approx(c(20, 50, 100, 200), c(0.04, 0.08, 0.25, 0.59),
                                  xout = n, rule = 2)$y, 0.5)
  its <- 40
  scale <- N / 200
  switch(kind,
    cr_fit   = per_it * its * scale,
    dd_fit   = per_it * its * scale * 1.5,
    pipeline = stats::approx(c(20, 50, 100), c(33, 66, 266), xout = n, rule = 2)$y,
    fhat     = if (identical(sampler, "bdi"))
                 stats::approx(c(50, 100), c(1, 23), xout = n, rule = 2)$y * 3
               else stats::approx(c(50, 100), c(0.8, 51), xout = n, rule = 2)$y * 6,
    ref      = if (identical(sampler, "dd")) 60 else 2,
    10)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# --- decision-rule thresholds (fixed before the main tier) -----------------
VAL_RULES <- list(
  cr = list(delta_ell_median = 0.1, delta_ell_p90 = 0.5,
            bias_se = 0.1, bias_ci_reject = 0.25, mc_sd_se = 0.25,
            e_median = 0.1, e_p90 = 0.5, e_bdi = 1e-6,
            converged_frac = 0.95, ess_min = 20, n_scaling_slope = -0.7),
  dd = list(delta_ell_median = 0.2, delta_ell_p90 = 1.0,
            bias_se = 0.2, bias_ci_reject = 0.5, mc_sd_se = 0.5,
            e_median = 0.2, e_p90 = 1.0, e_bdi = 1e-6,
            converged_frac = 0.95, ess_min = 20, n_scaling_slope = -0.7)
)

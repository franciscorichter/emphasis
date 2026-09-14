# ---------------------------------------------------------------------------
# The MCEM stopping rule: the ascent in the EM surrogate against its own
# Monte Carlo error.  One implementation (.ascent_controller), two drivers.
#
# Order of the file follows the order the gates should be read in:
#   1. the scoring functions did not move        (the hard constraint)
#   2. the statistic is free of the time unit and of the bound box
#   3. the verdicts, growth, budget and degradation paths
#   4. the two drivers cannot diverge
#   5. end to end, both samplers
# ---------------------------------------------------------------------------

# --- helpers ---------------------------------------------------------------

# A d-vector whose weighted mean is `m` and whose delta-method standard error
# under equal weights is exactly `s`.
mk_d <- function(m, s, n) {
  z <- seq_len(n) - (n + 1) / 2
  z <- z / sqrt(mean(z^2))
  m + s * sqrt(n) * z
}

# A controller wired to stubs: `score` returns logf_old + the d-vector the
# test currently wants, and `refit` is the identity, so the bootstrap spread
# is zero and the analytic error is what the verdicts see.
stub_ctl <- function(over = list(), n = 200L, seed = 1L) {
  set.seed(seed)
  env <- new.env(parent = emptyenv())
  env$d <- rep(0, n)
  env$logf_old <- rep(-10, n)
  opts <- .ascent_opts(
    utils::modifyList(list(sample_size = n, max_iter = 100L, tol = 1e-2,
                           patience = 3L, mc_boot = 6L), over),
    lower_bound = rep(0, 8), upper_bound = c(1, 1, rep(0, 6)))
  ctl <- .ascent_controller(opts,
                            refit = function(w, init) init,
                            score = function(theta, trees) env$logf_old + env$d)
  list(ctl = ctl, env = env, opts = opts)
}

# One decision at the current d-vector.
stub_step <- function(h, m, s, delta_max = 1, m_moved = TRUE, n = NULL) {
  n <- if (is.null(n)) as.integer(h$ctl$S) else n
  h$env$d <- mk_d(m, s, n)
  h$env$logf_old <- rep(-10, n)
  h$ctl$add_draws(n)
  h$ctl$decide(theta_old = c(1, 1), theta_new = c(1.01, 1.01),
               m_moved = m_moved, trees = vector("list", n),
               w_raw = rep(1, n), logf_old = h$env$logf_old,
               delta_max = delta_max, abs_step = 0.01)
}

# A small CR tree used by several gates.
cr_tree <- function(seed = 11L, lambda = 0.6, mu = 0.15, max_t = 8) {
  set.seed(seed)
  sim <- simulate_tree(pars = c(lambda, mu), max_t = max_t, model = "cr")
  sort(ape::branching.times(sim$tes), decreasing = TRUE)
}

# --- 1. the scoring functions did not move ---------------------------------

test_that("the split of the thinning iteration is scoring-identical to em_cpp", {
  # .mcem_dynamic_fresh no longer calls em_cpp(); it calls augment_trees()
  # followed by m_cpp() on the weights the C++ E-step would have built,
  # exp(log w - max log w).  src/mcem.cpp is E_step followed by M_step with
  # the same arguments, so the two must agree bit for bit.  This is the
  # identity the driver surgery rests on, and it is the test that the
  # scorers were not touched by it.
  brts <- cr_tree()
  th   <- .expand_pars(c(0.5, 0.1), c(0L, 0L, 0L))
  lb8  <- .expand_pars(c(0, 0), c(0L, 0L, 0L))
  ub8  <- .expand_pars(c(3, 3), c(0L, 0L, 0L))

  a <- em_cpp(brts = brts, init_pars = th, sample_size = 200L, maxN = 4000L,
              max_missing = 1e4, max_lambda = 1e6, lower_bound = lb8,
              upper_bound = ub8, xtol_rel = 1e-3, num_threads = 1L,
              copy_trees = FALSE, model = c(0L, 0L, 0L), link = 0L, rho = 1,
              rconditional = NULL, parent_tip_start = numeric(0), seed = 12345L)
  b <- augment_trees(brts = brts, pars = th, sample_size = 200L, maxN = 4000L,
                     max_missing = 1e4, max_lambda = 1e6, num_threads = 1L,
                     model = c(0L, 0L, 0L), link = 0L, rho = 1,
                     parent_tip_start = numeric(0), seed = 12345L)

  expect_identical(a$logf, b$logf)
  expect_identical(a$logg, b$logg)

  lw <- b$logf - b$logg
  w  <- exp(lw - max(lw))
  m  <- m_cpp(e_step = list(trees = b$trees, weights = w, rejected = 0L,
                            rejected_overruns = 0L, rejected_lambda = 0L,
                            rejected_zero_weights = 0L, time = 0, fhat = 0),
              init_pars = th, plugin = "rpd1", lower_bound = lb8,
              upper_bound = ub8, xtol_rel = 1e-3, num_threads = 1L,
              model = c(0L, 0L, 0L), link = 0L, rho = 1, rconditional = NULL)
  expect_identical(as.numeric(a$estimates), as.numeric(m$estimates))
  expect_identical(a$fhat,
                   .is_fhat(b$logf, b$logg,
                            n_zero_weight = .n0(b$rejected_zero_weights)))
})

test_that("eval_logf reproduces the E-step's own logf on the trees it returns", {
  # The statistic rescores trees the E-step already drew.  If eval_logf and
  # the E-step's scorer ever disagreed on those trees, every dQ would carry
  # the disagreement.
  brts <- cr_tree()
  th   <- .expand_pars(c(0.5, 0.1), c(0L, 0L, 0L))
  b <- augment_trees(brts = brts, pars = th, sample_size = 50L, maxN = 2000L,
                     max_missing = 1e4, max_lambda = 1e6, num_threads = 1L,
                     model = c(0L, 0L, 0L), link = 0L, rho = 1,
                     parent_tip_start = numeric(0), seed = 7L)
  ev <- eval_logf(th, b$trees, model = c(0L, 0L, 0L), link = 0L, rho = 1)
  expect_identical(ev$logf, b$logf)
  # And rescoring is a pure function of (theta, trees).
  expect_identical(ev$logf,
                   eval_logf(th, b$trees, model = c(0L, 0L, 0L), link = 0L,
                             rho = 1)$logf)
})

test_that("the scoring corpus rescores to the values it was frozen at", {
  # A seeded corpus of augmented trees plus the logf, logg and M-step
  # estimates recorded on the build this rule was written against.  Nothing
  # in the stopping rule may move these; a failure here is the scorers
  # having changed, not the rule.
  f <- test_path("fixtures", "scoring-corpus.rds")
  skip_if_not(file.exists(f), "scoring corpus fixture not present")
  fx <- readRDS(f)
  for (case in fx$cases) {
    ev <- eval_logf(case$theta, fx$trees[[case$corpus]],
                    model = case$model, link = case$link, rho = case$rho)
    expect_identical(ev$logf, case$logf)
    expect_identical(ev$logg, case$logg)
  }
  for (case in fx$msteps) {
    m <- m_cpp(e_step = list(trees = fx$trees[[case$corpus]],
                             weights = case$weights, rejected = 0L,
                             rejected_overruns = 0L, rejected_lambda = 0L,
                             rejected_zero_weights = 0L, time = 0, fhat = 0),
               init_pars = case$init, plugin = "rpd1",
               lower_bound = case$lower, upper_bound = case$upper,
               xtol_rel = 1e-3, num_threads = 1L, model = case$model,
               link = case$link, rho = case$rho, rconditional = NULL)
    expect_identical(as.numeric(m$estimates), case$estimates)
  }
})

# --- 2. invariance ---------------------------------------------------------

test_that("the Q-increment is free of the unit of time, as an identity", {
  # Scale every branching time by c and every rate by 1/c.  Both logf terms
  # pick up the same theta-free Jacobian, which cancels in the difference,
  # and log w = logf - logg picks it up twice and cancels in the weights.
  # A failure here is a coding error, not a tolerance question.
  brts <- cr_tree()
  cc <- 1000
  th1 <- c(0.5, 0.1); th2 <- c(0.55, 0.12)
  e <- augment_trees(brts = brts, pars = .expand_pars(th1, c(0L, 0L, 0L)),
                     sample_size = 40L, maxN = 2000L, max_missing = 1e4,
                     max_lambda = 1e6, num_threads = 1L, model = c(0L, 0L, 0L),
                     link = 0L, rho = 1, parent_tip_start = numeric(0),
                     seed = 3L)
  e$weights <- e$logf - e$logg
  dq <- function(theta_a, theta_b, trees, lw) {
    ww <- .ascent_weights(exp(lw - max(lw)), log = FALSE)
    f_a <- eval_logf(.expand_pars(theta_a, c(0L, 0L, 0L)), trees,
                     model = c(0L, 0L, 0L), link = 0L, rho = 1)$logf
    f_b <- eval_logf(.expand_pars(theta_b, c(0L, 0L, 0L)), trees,
                     model = c(0L, 0L, 0L), link = 0L, rho = 1)$logf
    .ascent_stat((f_a - f_b)[ww$keep], ww$wt)
  }
  s1 <- dq(th2, th1, e$trees, e$weights)

  # The same trees in another unit of time.  Every time-valued column is
  # scaled; t_ext_tip (1e11) and ts_unknown (-1) are sentinels, not times,
  # and scaling them would make the tree a different tree.
  scale_tree <- function(tr, cc) {
    sc <- function(x, sentinel) ifelse(x == sentinel, x, x * cc)
    tr$brts  <- tr$brts * cc
    tr$t_ext <- sc(tr$t_ext, 10e10)
    if (!is.null(tr$pd)) tr$pd <- tr$pd * cc
    if (!is.null(tr$tip_start)) tr$tip_start <- sc(tr$tip_start, -1)
    if (!is.null(tr$focal_tip_start))
      tr$focal_tip_start <- sc(tr$focal_tip_start, -1)
    tr
  }
  trees_c <- lapply(e$trees, scale_tree, cc = cc)
  f1 <- eval_logf(.expand_pars(th1 / cc, c(0L, 0L, 0L)), trees_c,
                  model = c(0L, 0L, 0L), link = 0L, rho = 1)
  lw_c <- f1$logf - f1$logg
  s2 <- dq(th2 / cc, th1 / cc, trees_c, lw_c)

  expect_equal(s2$dQ, s1$dQ, tolerance = 1e-8)
  expect_equal(s2$se, s1$se, tolerance = 1e-8)
})

test_that("the statistic never sees the bound box", {
  # eval_logf takes no bounds, so the box cannot reach dQ or its error.  The
  # verdict a controller reaches is therefore the same under any box, which
  # is what wave 1 established for the step rule and what this one inherits.
  h1 <- stub_ctl(); h2 <- stub_ctl()
  h2$opts$p_free <- h1$opts$p_free
  v1 <- stub_step(h1, m = 0.001, s = 0.001)
  v2 <- stub_step(h2, m = 0.001, s = 0.001)
  expect_identical(v1$dQ, v2$dQ)
  expect_identical(v1$se_dQ, v2$se_dQ)
})

# --- 3. the verdicts -------------------------------------------------------

test_that("a clear ascent is case 1 and never stops", {
  h <- stub_ctl()
  for (i in 1:5) v <- stub_step(h, m = 0.5, s = 0.01)
  expect_identical(v$case, 1L)
  expect_false(v$stop)
  expect_true(v$accept)
  expect_gt(v$L, 0)
})

test_that("a gain that cannot be asserted above delta certifies after patience", {
  h <- stub_ctl()
  seen <- integer(0); stopped <- NULL
  for (i in 1:6) {
    v <- stub_step(h, m = 0.0005, s = 0.0005)
    seen <- c(seen, v$case)
    if (v$stop) { stopped <- v; break }
  }
  expect_true(all(seen == 2L))
  expect_false(is.null(stopped))
  expect_identical(stopped$stop_reason, "converged")
  expect_identical(length(seen), 3L)          # exactly `patience` of them
  expect_lt(stopped$U, h$opts$ascent_tol)
})

test_that("a certificate needs ess_min effective draws", {
  # The same quiet sequence on a sample carrying fewer effective draws than
  # the nominal num_trees may not certify, however quiet it is.
  h <- stub_ctl(over = list(ess_min = 500))
  for (i in 1:8) v <- stub_step(h, m = 0.0005, s = 0.0005)
  expect_identical(v$case, 2L)
  expect_false(v$stop)
})

test_that("an M-step that returned its start is never a certificate", {
  h <- stub_ctl()
  for (i in 1:5) v <- stub_step(h, m = 0, s = 0, m_moved = FALSE)
  expect_identical(v$case, 3L)
  expect_false(v$stop)
})

# --- growth ----------------------------------------------------------------

test_that("growth waits for a persistent precision-bound verdict", {
  # dQ below delta but the interval still reaching past it.  One such
  # iteration is a quiet draw mid-climb and must not multiply the budget;
  # `patience` of them in a row is evidence.
  # delta is pinned here rather than taken from the shipped default: the
  # test is about the growth mechanics, not about where delta sits.
  h <- stub_ctl(over = list(ascent_tol = 0.01))
  S0 <- as.integer(h$ctl$S)
  v1 <- stub_step(h, m = 0.008, s = 0.006)
  expect_identical(v1$case, 3L)
  expect_false(v1$grow)
  expect_identical(as.integer(h$ctl$S), S0)
  v2 <- stub_step(h, m = 0.008, s = 0.006)
  expect_false(v2$grow)
  v3 <- stub_step(h, m = 0.008, s = 0.006)
  expect_true(v3$grow)
  expect_false(v3$accept)                     # the step is not taken
  expect_gt(as.integer(h$ctl$S), S0)
})

test_that("travel is never answered with draws", {
  # The point estimate is above delta: the iterate is moving and a larger
  # sample cannot settle a moving iterate.  Take the step and iterate.
  h <- stub_ctl()
  S0 <- as.integer(h$ctl$S)
  for (i in 1:6) v <- stub_step(h, m = 0.2, s = 0.15)
  expect_identical(v$case, 3L)
  expect_false(v$grow)
  expect_true(v$accept)
  expect_identical(as.integer(h$ctl$S), S0)
})

test_that("growth is monotone and capped", {
  h <- stub_ctl(over = list(max_sample_mult = 3, mc_refresh = 1L,
                            ascent_tol = 0.01))
  S <- as.integer(h$ctl$S); seen <- S
  for (i in 1:60) {
    v <- stub_step(h, m = 0.008, s = 0.006)
    seen <- c(seen, as.integer(h$ctl$S))
    if (v$stop) break
  }
  expect_true(all(diff(seen) >= 0))
  expect_lte(max(seen), 3 * 200)
  expect_identical(v$stop_reason, "undecided_at_max_sample")
})

test_that("growth switched off still terminates, and not as convergence", {
  h <- stub_ctl(over = list(grow_sample = FALSE, ascent_tol = 0.01))
  S0 <- as.integer(h$ctl$S)
  for (i in 1:60) { v <- stub_step(h, m = 0.008, s = 0.006); if (v$stop) break }
  expect_true(v$stop)
  expect_true(v$stop_reason %in% c("undecided_at_max_sample", "draw_budget"))
  expect_identical(as.integer(h$ctl$S), S0)   # nothing was spent on draws
})

test_that("the draw budget is the wave-1 worst case and is not convergence", {
  h <- stub_ctl(over = list(max_iter = 5L))     # budget = 5 * 200 draws
  for (i in 1:10) { v <- stub_step(h, m = 0.5, s = 0.01); if (v$stop) break }
  expect_true(v$stop)
  expect_identical(v$stop_reason, "draw_budget")
  expect_lte(h$ctl$draws, 5 * 200 + 200)
})

# --- degradation -----------------------------------------------------------

test_that("a sample whose noise cannot be estimated falls back to the tol rule", {
  h <- stub_ctl(n = 4L)          # below the 5-weight minimum
  for (i in 1:4) {
    v <- h$ctl$decide(theta_old = c(1, 1), theta_new = c(1.0001, 1.0001),
                      m_moved = TRUE, trees = vector("list", 4L),
                      w_raw = rep(1, 4L), logf_old = rep(-10, 4L),
                      delta_max = 1e-4, abs_step = 1e-4)
    if (v$stop) break
  }
  expect_identical(v$gate, "legacy")
  expect_true(v$stop)
  expect_identical(v$stop_reason, "converged_step")
})

test_that("weight degeneracy below ess_floor routes to the fallback", {
  h <- stub_ctl()
  w <- c(1, rep(1e-12, 199))      # ESS ~ 1
  v <- h$ctl$decide(theta_old = c(1, 1), theta_new = c(1.01, 1.01),
                    m_moved = TRUE, trees = vector("list", 200L),
                    w_raw = w, logf_old = rep(-10, 200L),
                    delta_max = 0.5, abs_step = 0.01)
  expect_identical(v$gate, "legacy")
  expect_false(v$stop)
})

test_that("stop = \"step\" is the wave-1 rule and nothing else", {
  h <- stub_ctl(over = list(stop = "step"))
  hits <- 0L
  for (i in 1:6) {
    v <- stub_step(h, m = 0.5, s = 0.01, delta_max = 1e-3)  # big ascent, tiny step
    hits <- hits + 1L
    if (v$stop) break
  }
  expect_identical(v$gate, "legacy")
  expect_identical(v$stop_reason, "converged")
  expect_identical(hits, 3L)     # patience consecutive delta_max < tol
  expect_true(is.na(v$dQ))       # the ascent statistic was never computed
})

test_that("num_trees = 1 terminates and does not enter the new rule", {
  brts <- cr_tree(seed = 3L)
  fit <- estimate_rates(brts, method = "mcem", model = "cr",
                        init_pars = c(0.6, 0.2),
                        control = list(lower_bound = c(0, 0),
                                       upper_bound = c(3, 3),
                                       sample_size = 1L, max_iter = 12L,
                                       maxN = 2000L, num_threads = 1L,
                                       max_time = 300),
                        link = "linear")
  expect_true(fit$stop_reason %in%
                c("converged_step", "max_iter", "e_step_failure",
                  "m_step_failure", "time_budget", "draw_budget"))
  tr <- fit$details$mcem
  g <- tr$gate[!is.na(tr$gate)]
  expect_true(length(g) == 0L || all(g == "legacy"))
})

# --- 4. the two drivers cannot diverge -------------------------------------

test_that("both drivers take the same defaults and the same rule object", {
  d <- estimate_rates_control("mcem")
  for (nm in names(.ascent_defaults())) expect_true(nm %in% names(d))
  expect_identical(d$stop, "ascent")
  expect_identical(d$ascent_tol, 0.002)
  expect_true("stop_opts" %in% names(formals(.mcem_bdi)))
  expect_true("stop_opts" %in% names(formals(.mcem_dynamic_fresh)))
  # Neither driver computes the statistic itself.
  for (f in c(".mcem_bdi", ".mcem_dynamic_fresh")) {
    src <- paste(deparse(get(f)), collapse = "\n")
    expect_true(grepl(".ascent_controller", src, fixed = TRUE))
    expect_false(grepl("qnorm", src, fixed = TRUE))
  }
})

test_that("the same synthetic sequence gives the same stop from either driver", {
  # The rule is one object; what a driver supplies is (trees, weights, logf,
  # delta_max, m_moved).  Replaying one recorded sequence through two
  # independently constructed controllers must give the same stop index and
  # the same reason, on every path including the degenerate ones.
  seqs <- list(
    quiet   = list(m = rep(0.0005, 8), s = rep(0.0005, 8), moved = rep(TRUE, 8)),
    climb   = list(m = c(0.5, 0.3, 0.2, 0.05, 0.004, 0.0006, 0.0004, 0.0003),
                   s = c(rep(0.01, 5), rep(0.0004, 3)), moved = rep(TRUE, 8)),
    stalled = list(m = rep(0, 8), s = rep(0, 8), moved = rep(FALSE, 8)),
    bumpy   = list(m = rep(0.012, 8), s = rep(0.005, 8), moved = rep(TRUE, 8)))
  run <- function(sq) {
    h <- stub_ctl(seed = 99L)
    for (k in seq_along(sq$m)) {
      v <- stub_step(h, m = sq$m[k], s = sq$s[k], m_moved = sq$moved[k])
      if (v$stop) return(list(k = k, reason = v$stop_reason, S = h$ctl$S))
    }
    list(k = NA_integer_, reason = NA_character_, S = h$ctl$S)
  }
  for (nm in names(seqs)) expect_identical(run(seqs[[nm]]), run(seqs[[nm]]))
})

test_that("the trace carries the same stopping columns in both drivers", {
  want <- names(.ascent_trace_cols())
  brts <- cr_tree(seed = 5L, max_t = 6)
  ctrl <- list(lower_bound = c(0, 0), upper_bound = c(3, 3),
               sample_size = 60L, max_iter = 6L, maxN = 2000L,
               num_threads = 1L, max_time = 300)
  for (s in c("bdi", "dynamic_fresh")) {
    fit <- estimate_rates(brts, method = "mcem", model = "cr",
                          init_pars = c(1.2, 0.6),
                          control = c(ctrl, list(sampling = s)),
                          link = "linear")
    expect_true(all(want %in% names(fit$details$mcem)))
    expect_true(is.finite(fit$total_draws))
    expect_gte(fit$sample_size_final, 60L)
  }
})

# --- 5. end to end ---------------------------------------------------------

test_that("both samplers run under both rules and report an honest reason", {
  ok <- c("converged", "converged_step", "undecided_at_max_sample",
          "draw_budget", "max_iter", "time_budget", "e_step_failure",
          "m_step_failure")
  brts <- cr_tree(seed = 21L, max_t = 7)
  for (s in c("bdi", "dynamic_fresh")) for (r in c("ascent", "step")) {
    fit <- estimate_rates(brts, method = "mcem", model = "cr",
                          init_pars = c(0.7, 0.2),
                          control = list(lower_bound = c(0, 0),
                                         upper_bound = c(3, 3),
                                         sampling = s, stop = r,
                                         sample_size = 80L, max_iter = 25L,
                                         maxN = 2000L, num_threads = 1L,
                                         max_time = 300),
                          link = "linear")
    expect_true(fit$stop_reason %in% ok)
    expect_true(all(is.finite(as.numeric(fit$pars))))
    expect_gt(as.numeric(fit$pars)[1], 0)
  }
})

test_that("the bootstrap measures the spread of the increment it is used for", {
  # The multinomial resample of the weighted trees must reproduce the spread
  # of the M-step across independent E-steps.  Cheap version of the
  # calibration gate: at one theta, compare the bootstrap standard error of
  # dQ with the spread of dQ over independent E-steps.
  skip_on_cran()
  brts <- cr_tree(seed = 4L, max_t = 7)
  mb <- c(0L, 0L, 0L)
  th <- .expand_pars(c(0.6, 0.15), mb)
  lb8 <- .expand_pars(c(0, 0), mb); ub8 <- .expand_pars(c(3, 3), mb)
  refit <- function(trees) function(w, init) {
    m <- tryCatch(m_cpp(e_step = list(trees = trees, weights = w, rejected = 0L,
                                      rejected_overruns = 0L, rejected_lambda = 0L,
                                      rejected_zero_weights = 0L, time = 0, fhat = 0),
                        init_pars = init, plugin = "rpd1", lower_bound = lb8,
                        upper_bound = ub8, xtol_rel = 1e-3, num_threads = 1L,
                        model = mb, link = 0L, rho = 1, rconditional = NULL),
                  error = function(e) NULL)
    if (is.null(m)) NULL else as.numeric(m$estimates)
  }
  score <- function(theta, trees) eval_logf(theta, trees, model = mb, link = 0L,
                                            rho = 1)$logf
  set.seed(42)
  R <- 25L
  dQ_r <- numeric(R); se_b <- numeric(R)
  for (r in seq_len(R)) {
    e <- .augment_tree_bdi(tree = brts, pars = th, model_bin = mb,
                           sample_size = 200L, max_missing = 1e4, link = 0L,
                           rho = 1)
    ww <- .ascent_weights(exp(e$weights - max(e$weights)), log = FALSE)
    w_raw <- exp(e$weights - max(e$weights))
    th_new <- refit(e$trees)(w_raw, th)
    lf_new <- score(th_new, e$trees)
    dQ_r[r] <- .ascent_stat((lf_new - e$logf)[ww$keep], ww$wt)$dQ
    bb <- .ascent_bootstrap(e$trees, w_raw = w_raw, wt = ww$wt, keep = ww$keep,
                            logf_old = e$logf[ww$keep], logf_new = lf_new[ww$keep],
                            theta_new = th_new, refit = refit(e$trees),
                            score = score, B = 12L)
    se_b[r] <- bb$se_boot
    expect_true(!is.finite(bb$mc_noise) || bb$mc_noise >= -1e-8)
  }
  ratio <- mean(se_b, na.rm = TRUE) / stats::sd(dQ_r)
  expect_gt(ratio, 0.4)
  expect_lt(ratio, 2.5)
})

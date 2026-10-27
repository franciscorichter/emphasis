# Latent marks (dev/latent-marks/design.md): an unseen state every lineage
# carries, inherited at splits, that multiplies its rates.  The likelihood sums
# the marks out by pruning on the lineage forest (Model::loglik_marks).

nk  <- emphasis:::.resolve_model("nk")
lgm <- function(pars, trees, mb, S = 1L, log_clam = 0, log_cmu = 0, C = 1, P = numeric(0), pi = 1,
                rho = 1) {
  emphasis:::eval_logf_marks(pars, trees, model = mb, link = 1L, rho = rho, S = S,
                             log_clam = log_clam, log_cmu = log_cmu, C = C, P = P, pi = pi)
}
aug <- function(phy, pars, mb, n, seed, rho = 1) {
  brts <- emphasis:::.extract_brts(phy)
  emphasis:::augment_trees(brts, emphasis:::.expand_pars(pars, mb), sample_size = n, maxN = 50L * n,
                           max_missing = 500L, max_lambda = 1e6, num_threads = 1L, model = mb,
                           link = 1L, rho = rho, parent_tip_start = emphasis:::.pts(brts), seed = seed,
                           parent_id = emphasis:::.pid(brts))$trees
}

test_that("with one mark the likelihood is eval_logf's, on observed and augmented trees", {
  set.seed(4)
  phy <- ape::rphylo(25L, 0.7, 0)
  fr  <- emphasis:::.observed_frame(phy)
  for (m in c("nk", "ned", "nedc", "nd", "dd", "k")) {
    mb <- emphasis:::.resolve_model(m)
    np <- length(emphasis:::.par_names(mb))
    p  <- emphasis:::.expand_pars(c(log(0.5), rep(0.02, np / 2 - 1), log(0.1), rep(-0.03, np / 2 - 1)), mb)
    expect_equal(lgm(p, list(fr), mb), emphasis:::eval_logf(p, list(fr), model = mb, link = 1L, rho = 1)$logf,
                 tolerance = 1e-10, info = m)
  }
  p  <- c(log(0.6), -0.01, 0.08, log(0.3), 0.01, -0.05)
  for (rho in c(1, 0.8)) {
    tr <- aug(phy, p, nk, 10L, seed = 5L, rho = rho)
    pe <- emphasis:::.expand_pars(p, nk)
    expect_equal(lgm(pe, tr, nk, rho = rho),
                 emphasis:::eval_logf(pe, tr, model = nk, link = 1L, rho = rho)$logf, tolerance = 1e-10)
  }
})

test_that("marks that do not change the rates leave the likelihood unchanged", {
  set.seed(6)
  phy <- ape::rphylo(15L, 0.7, 0)
  pe  <- emphasis:::.expand_pars(c(log(0.6), -0.01, 0.08, log(0.3), 0.01, -0.05), nk)
  tr  <- aug(phy, c(log(0.6), -0.01, 0.08, log(0.3), 0.01, -0.05), nk, 5L, seed = 2L)
  S  <- 3L
  Cm <- array(stats::runif(S^3), c(S, S, S))                  # Cm[m, a, b]
  for (m in seq_len(S)) Cm[m, , ] <- Cm[m, , ] / sum(Cm[m, , ])
  expect_equal(lgm(pe, tr, nk, S = S, log_clam = rep(0, S), log_cmu = rep(0, S),
                   C = as.vector(aperm(Cm, c(3, 2, 1))), pi = rep(1 / S, S)),
               lgm(pe, tr, nk), tolerance = 1e-10)
})

# The complete-data log density of one mark history, from the lineage table:
# the one-mark likelihood plus, per stretch, the multipliers' event terms and
# their change to the integrals, plus the kernel at every split and the root.
brute_marks <- function(pe, tree, S, lc, lm, Cbar, pi) {
  tab <- emphasis:::lineage_table_cpp(tree)
  tab <- tab[order(tab$lineage, tab$t0), ]
  l0  <- lgm(pe, list(tree), nk)
  eta_l <- pe[1] + pe[2] * tab$N + pe[11] * tab$K
  eta_m <- pe[5] + pe[6] * tab$N + pe[12] * tab$K
  dt <- tab$t1 - tab$t0
  # stretch of each row: the splits of its lineage that end an earlier row
  tab$j <- stats::ave(tab$event == 1L, tab$lineage, FUN = function(e) c(0L, cumsum(e)[-length(e)]))
  key <- paste(tab$lineage, tab$j)
  st  <- unique(key)
  first_stretch <- function(l) paste(l, 0L)
  born <- tapply(tab$t0, tab$lineage, min)
  sp  <- tab[tab$event == 1L, ]                                   # splits: lineage, stretch, time
  dau <- vapply(seq_len(nrow(sp)), function(r) {
    cand <- as.integer(names(born))[abs(born - sp$t1[r]) < 1e-12 & as.integer(names(born)) != sp$lineage[r]]
    cand[1]
  }, 1L)
  grid <- as.matrix(expand.grid(rep(list(seq_len(S)), length(st) + 1L)))
  vals <- apply(grid, 1, function(g) {
    m <- stats::setNames(g[-1], st); r <- g[1]
    mr <- m[key]
    v <- l0 + log(pi[r]) + log(Cbar[r, m[first_stretch(0L)], m[first_stretch(1L)]])
    v <- v - sum((exp(lc[mr]) - 1) * exp(eta_l) * dt) - sum((exp(lm[mr]) - 1) * exp(eta_m) * dt)
    v <- v + sum(lc[mr[tab$event == 1L]]) + sum(lm[mr[tab$event == 2L]])
    for (s in seq_len(nrow(sp))) {
      a <- m[paste(sp$lineage[s], sp$j[s])]
      v <- v + log(Cbar[a, m[paste(sp$lineage[s], sp$j[s] + 1L)], m[first_stretch(dau[s])]])
    }
    v
  })
  mx <- max(vals); mx + log(sum(exp(vals - mx)))
}

test_that("pruning equals the sum over every mark history", {
  set.seed(9)
  phy <- ape::rphylo(4L, 0.7, 0)
  p   <- c(log(0.6), -0.02, 0.1, log(0.4), 0.01, -0.05)
  pe  <- emphasis:::.expand_pars(p, nk)
  trees <- c(list(emphasis:::.observed_frame(phy)), aug(phy, p, nk, 20L, seed = 7L))
  nst <- vapply(trees, function(tr) {
    tab <- emphasis:::lineage_table_cpp(tr); length(unique(tab$lineage)) + sum(tab$event == 1L)
  }, 1)
  trees <- trees[nst <= 11]
  expect_true(any(vapply(trees, function(tr) any(emphasis:::lineage_table_cpp(tr)$event == 2L), TRUE)))
  S  <- 2L
  lc <- c(0, 0.7); lm <- c(0, -0.4); pi <- c(0.3, 0.7)
  C  <- array(c(0.5, 0.1, 0.2, 0.2,  0.1, 0.3, 0.1, 0.5), c(2, 2, 2))   # C[a, b, m]
  Cm <- aperm(C, c(3, 1, 2))                                           # Cm[m, a, b]
  Cbar <- (Cm + aperm(Cm, c(1, 3, 2))) / 2
  for (tr in unique(trees)[seq_len(min(5, length(unique(trees))))]) {
    got <- lgm(pe, list(tr), nk, S = S, log_clam = lc, log_cmu = lm,
               C = as.vector(aperm(Cm, c(3, 2, 1))), pi = pi)
    expect_equal(got, brute_marks(pe, tr, S, lc, lm, Cbar, pi), tolerance = 1e-9)
  }
  # independent daughters: P[m, a]
  P  <- matrix(c(0.8, 0.3, 0.2, 0.7), 2)
  Ci <- array(0, c(2, 2, 2)); for (m in 1:2) Ci[m, , ] <- outer(P[m, ], P[m, ])
  tr <- trees[[1]]
  expect_equal(lgm(pe, list(tr), nk, S = S, log_clam = lc, log_cmu = lm, P = as.vector(t(P)), pi = pi),
               brute_marks(pe, tr, S, lc, lm, Ci, pi), tolerance = 1e-9)
})

test_that("the linear link and a wrong kernel size are refused", {
  phy <- ape::rphylo(6L, 0.7, 0)
  fr  <- emphasis:::.observed_frame(phy)
  pe  <- emphasis:::.expand_pars(c(0.5, 0, 0.01, 0.1, 0, 0), nk)
  expect_error(emphasis:::eval_logf_marks(pe, list(fr), model = nk, link = 0L, rho = 1, S = 1L,
                                          log_clam = 0, log_cmu = 0, C = 1, P = numeric(0), pi = 1),
               "exponential link")
  expect_error(lgm(pe, list(fr), nk, S = 2L, log_clam = c(0, 0), log_cmu = c(0, 0), C = 1, pi = c(0.5, 0.5)),
               "do not match S")
})

test_that("the latent-rate preset reduces to one mark at sigma = 0 and converges in the grid", {
  set.seed(12)
  phy <- ape::rphylo(30L, 0.7, 0)
  p   <- c(log(0.6), -0.01, 0.05, log(0.2), 0, 0)
  pe  <- emphasis:::.expand_pars(p, nk)
  tr  <- aug(phy, p, nk, 3L, seed = 4L)
  l1  <- lgm(pe, tr, nk)
  expect_equal(emphasis:::.eval_logf_marks(pe, tr, nk, emphasis:::.latent_rate_marks(0)), l1, tolerance = 1e-12)
  # a walk too narrow to move a rate is the one-mark likelihood to rounding
  expect_equal(emphasis:::.eval_logf_marks(pe, tr, nk, emphasis:::.latent_rate_marks(1e-4, S = 21L)), l1,
               tolerance = 1e-5)
  # the grid: 21 and 61 points at sigma = 0.3 agree to a hundredth of a nat
  a <- emphasis:::.eval_logf_marks(pe, tr, nk, emphasis:::.latent_rate_marks(0.3, S = 21L))
  b <- emphasis:::.eval_logf_marks(pe, tr, nk, emphasis:::.latent_rate_marks(0.3, S = 61L))
  expect_true(all(abs(a - b) < 0.01))
  expect_true(all(a != l1))
})

# ---- the simulator (phase 2) ----

sim1 <- function(pe, mb, T, seed, ...) {
  emphasis:::simulate_div_tree_cpp(pe, mb, T, 5000L, 0L, 1L, seed = seed, ...)
}

test_that("without marks the simulator draws the same tree as before, seed for seed", {
  cr <- emphasis:::.resolve_model("cr")
  pe <- emphasis:::.expand_pars(c(log(0.4), log(0.1)), cr)
  for (s in 1:5) {
    a <- sim1(pe, cr, 6, s)
    b <- sim1(pe, cr, 6, s, marks_S = 1L, log_clam = 0, log_cmu = 0, marks_C = 1, marks_P = numeric(0), marks_pi = 1)
    expect_identical(a$Ltable, b$Ltable); expect_true(all(a$marks == 0L))
  }
})

tips_of <- function(n, pe, mb, T, seed0, ...) {
  vapply(seq_len(n), function(i) {
    r <- sim1(pe, mb, T, seed0 + i, ...)
    if (r$status == "done") sum(r$Ltable[, 4] == -1) else NA_real_
  }, 1)
}

test_that("marks that leave the rates alone leave the trees' distribution alone", {
  skip_on_cran()
  cr <- emphasis:::.resolve_model("cr")
  pe <- emphasis:::.expand_pars(c(log(0.4), log(0.1)), cr)
  a <- tips_of(400, pe, cr, 6, 1000L)
  b <- tips_of(400, pe, cr, 6, 5000L, marks_S = 3L, log_clam = rep(0, 3), log_cmu = rep(0, 3),
               marks_P = rep(1 / 3, 9), marks_pi = rep(1 / 3, 3))
  expect_gt(suppressWarnings(stats::ks.test(a[!is.na(a)], b[!is.na(b)])$p.value), 0.01)
})

test_that("a kernel that sends every lineage to a doubled rate simulates the doubled model", {
  skip_on_cran()
  cr <- emphasis:::.resolve_model("cr")
  pe <- emphasis:::.expand_pars(c(log(0.25), log(0.05)), cr)
  p2 <- emphasis:::.expand_pars(c(log(0.5), log(0.1)), cr)
  # mark 1 doubles both rates; every daughter, and both crown lineages, get mark 1
  a <- tips_of(400, p2, cr, 6, 1000L)
  b <- tips_of(400, pe, cr, 6, 5000L, marks_S = 2L, log_clam = c(0, log(2)), log_cmu = c(0, log(2)),
               marks_P = c(0, 1, 0, 1), marks_pi = c(1, 0))
  expect_gt(suppressWarnings(stats::ks.test(a[!is.na(a)], b[!is.na(b)])$p.value), 0.01)
  r <- sim1(pe, cr, 6, 7L, marks_S = 2L, log_clam = c(0, log(2)), log_cmu = c(0, log(2)),
            marks_P = c(0, 1, 0, 1), marks_pi = c(1, 0))
  expect_true(all(r$marks == 1L))
})

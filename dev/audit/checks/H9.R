## H9 — Are the M / D covariates evaluated by the estimator the covariates
##      that the forward simulator uses to generate data?
##
## Self-contained. Runs in < 30 s. num_threads = 1 throughout.
##
## Parts
##  A  topology-free 5-tip check of pd at the first observed split
##  B  forward-simulated nd tree with mu = 0 (complete tree == observed tree):
##     true M(t), E_focal, D_focal from the L-table vs pd/n, e_s - M in the
##     augmented data frame (max_missing = 0, so nothing is augmented)
##  C  complete-data log-likelihood on that tree: exact per-lineage value
##     (D constant on each inter-event segment) vs eval_logf, linear and
##     exponential links, as a function of beta_D
##  D  cr / dd (rho = 1): logf invariant to the pd column, all three links
##  E  ep_exp running-sum convention inside loglik vs the pd convention

.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressMessages(library(emphasis))
set.seed(9)

aug <- function(brts, pars8, model_bin, link, max_missing = 0L) {
  a <- emphasis:::augment_trees(brts, pars = as.numeric(pars8),
                                sample_size = 1L, maxN = 200L,
                                max_missing = as.integer(max_missing),
                                max_lambda = 1e6, num_threads = 1L,
                                model = as.integer(model_bin),
                                link = as.integer(link), rho = 1)
  stopifnot(length(a$trees) == 1L)
  a$trees[[1]]
}
logf <- function(pars8, df, model_bin, link)
  emphasis:::eval_logf(as.numeric(pars8), list(df),
                       model = as.integer(model_bin), link = as.integer(link),
                       rho = 1)$logf

cat("==== Part A: 5-tip tree, first split, topology-free ====\n")
brtsA <- c(1, 0.7, 0.4, 0.2)               # crown age 1, splits at fwd 0.3 0.6 0.8
dfA <- aug(brtsA, c(0.5, 0, 0, 0.1, 0, 0, 0, 0), c(1, 0, 1), 0)
print(dfA)
cat(sprintf("node 0 (fwd t = 0.3, n = 2): pkg pd = %.3f, pkg M = %.4f\n",
            dfA$pd[1], dfA$pd[1] / dfA$n[1]))
cat(sprintf("  true P(0.3^-) = 2 crown lineages x 0.3 = 0.600 (no topology needed), M = 0.300\n"))
cat(sprintf("  ratio pkg/true = %.3f\n", dfA$pd[1] / 0.6))
cat("  pkg pd at node i = (i+1)*s_i:", all.equal(dfA$pd, seq_along(dfA$pd) * dfA$brts), "\n")

cat("\n==== Part B: forward nd tree (mu = 0) -> true vs package covariates ====\n")
## nd, linear link, compact pars c(b0, bN, bD, g0, gN, gD); mu = 0 so the
## complete tree is the observed tree and no lineage is ever missing.
pars_c <- c(0.9, -0.02, 0.25, 0, 0, 0)
sim <- NULL
for (k in 1:50) {
  s <- simulate_tree(pars = pars_c, max_t = 3, model = "nd", link = "linear",
                     max_lin = 200, useDDD = FALSE)
  if (s$status == "done" && nrow(s$L) >= 8 && nrow(s$L) <= 40) { sim <- s; break }
}
stopifnot(!is.null(sim))
L <- sim$L
Tc <- 3
n_tips <- nrow(L)
cat(sprintf("tree: %d tips, crown age %g, all extant: %s\n", n_tips, Tc, all(L[, 4] == -1)))

## forward times
L_birth <- Tc - L[, 1]
ord <- order(L_birth)
L <- L[ord, ]; L_birth <- L_birth[ord]
ev <- L[-(1:2), , drop = FALSE]; ev_t <- L_birth[-(1:2)]      # observed events
brtsB <- sort(Tc - ev_t, decreasing = TRUE); brtsB <- c(Tc, brtsB)

## true tip_start of lineage `lab` just before forward time t
ts_at <- function(lab, t) {
  own <- L_birth[L[, 3] == lab]
  kids <- L_birth[L[, 2] == lab & L_birth < t]
  max(c(own, kids))
}
true_cov <- t(sapply(seq_along(ev_t), function(i) {
  t <- ev_t[i]
  alive <- L[L_birth < t, 3]
  ts <- sapply(alive, ts_at, t = t)
  P <- sum(t - ts); N <- length(alive); M <- P / N
  Ef <- t - ts_at(ev[i, 2], t)
  c(t = t, N = N, P = P, M = M, E_focal = Ef, D_focal = Ef - M)
}))

dfB <- aug(brtsB, emphasis:::.expand_pars(pars_c, c(1, 0, 1)), c(1, 0, 1), 0)
stopifnot(nrow(dfB) == n_tips - 1, all(dfB$parent_id == -1))
obs <- dfB[-nrow(dfB), ]                    # drop closing sentinel
stopifnot(isTRUE(all.equal(obs$brts, true_cov[, "t"], tolerance = 1e-5)))
M_pkg <- obs$pd / obs$n
E_pkg <- M_pkg                              # e_s() returns M when parent_id == -1
cmp <- data.frame(t = round(true_cov[, "t"], 3), N = true_cov[, "N"],
                  P_true = round(true_cov[, "P"], 3), P_pkg = round(obs$pd, 3),
                  M_true = round(true_cov[, "M"], 3), M_pkg = round(M_pkg, 3),
                  D_true = round(true_cov[, "D_focal"], 3), D_pkg = E_pkg - M_pkg)
print(cmp)
cat(sprintf("max |M_pkg - M_true| = %.3f (mean M_true = %.3f); mean |D_true| = %.3f, D_pkg == 0 at all %d events\n",
            max(abs(M_pkg - true_cov[, "M"])), mean(true_cov[, "M"]),
            mean(abs(true_cov[, "D_focal"])), nrow(obs)))
cat("N matches:", all(obs$n == true_cov[, "N"]), "\n")

cat("\n==== Part C: exact complete-data loglik vs eval_logf (mu = 0) ====\n")
## Segments: (0,t1), (t1,t2), ..., (tk, T). On each, N and every D_s are
## constant, so log f = sum_events log lambda_focal - sum_seg dt * sum_s lambda_s.
seg_ends <- c(ev_t, Tc)
true_loglik <- function(p8, link) {
  rate <- function(N, D) {
    eta <- p8[1] + p8[2] * N + p8[4] * D
    if (link == 0) max(0, eta) else exp(eta)
  }
  ll <- 0
  for (j in seq_along(seg_ends)) {
    t0 <- if (j == 1) 0 else seg_ends[j - 1]; t1 <- seg_ends[j]
    tm <- (t0 + t1) / 2
    alive <- L[L_birth < t1 - 1e-9, 3]
    ts <- sapply(alive, ts_at, t = tm)
    N <- length(alive); M <- sum(tm - ts) / N
    D <- (tm - ts) - M
    ll <- ll - (t1 - t0) * sum(sapply(D, function(d) rate(N, d)))
    if (j <= length(ev_t)) {
      foc <- which(alive == ev[j, 2])
      ll <- ll + log(rate(N, D[foc]))
    }
  }
  ll
}
res <- NULL
for (link in 0:1) {
  ## exp link: gamma_0 = -30 gives mu = 9e-14, i.e. no augmentation and a
  ## negligible (< 1e-12) mu contribution that the exact formula ignores
  g0 <- if (link == 1) -30 else 0
  dfC <- aug(brtsB, emphasis:::.expand_pars(replace(pars_c, 4, g0), c(1, 0, 1)), c(1, 0, 1), link)
  for (bD in c(-0.3, 0, 0.25, 0.6)) {
    p8 <- emphasis:::.expand_pars(replace(replace(pars_c, 3, bD), 4, g0), c(1, 0, 1))
    res <- rbind(res, data.frame(link = c("linear", "exp")[link + 1], beta_D = bD,
                                 logf_true = true_loglik(p8, link),
                                 logf_pkg = logf(p8, dfC, c(1, 0, 1), link)))
  }
}
res$diff <- res$logf_pkg - res$logf_true
print(res, row.names = FALSE)
lin <- res[res$link == "linear", ]
cat(sprintf("linear: range of logf_pkg over beta_D = %.3g (flat); range of logf_true = %.3f\n",
            diff(range(lin$logf_pkg)), diff(range(lin$logf_true))))
cat(sprintf("at beta_D = 0 (model reduces to dd): diff linear = %.3g, exp = %.3g\n",
            res$diff[res$link == "linear" & res$beta_D == 0],
            res$diff[res$link == "exp" & res$beta_D == 0]))

cat("\n==== Part D: cr / dd, rho = 1: logf invariant to pd ====\n")
for (mb in list(cr = c(0, 0, 0), dd = c(1, 0, 0))) for (link in 0:2) {
  pc <- if (sum(mb) == 0) c(0.6, 0.1) else c(0.6, -0.01, 0.1, 0)
  if (link == 1 && sum(mb) == 0) pc <- c(log(0.6), log(0.1))
  if (link == 1 && sum(mb) == 1) pc <- c(log(0.6), -0.01, log(0.1), 0)
  if (link == 2 && sum(mb) == 1) pc <- c(0.6, 0.05, 0.1, 0.02)
  p8 <- emphasis:::.expand_pars(pc, mb)
  d0 <- aug(brtsB, p8, mb, link, max_missing = 20L)
  d1 <- d0; d1$pd <- 0; d2 <- d0; d2$pd <- 1e3 * runif(nrow(d0))
  d3 <- d0; d3$tip_start <- runif(nrow(d0)); d3$focal_tip_start <- runif(nrow(d0))
  v <- c(logf(p8, d0, mb, link), logf(p8, d1, mb, link), logf(p8, d2, mb, link), logf(p8, d3, mb, link))
  cat(sprintf("%s link=%d n_aug=%d: logf = %.6f, pd->0: %.6f, pd->rand: %.6f, ts->rand: %.6f  identical=%s\n",
              names(which(sapply(list(cr = c(0,0,0), dd = c(1,0,0)), identical, mb))), link,
              sum(d0$parent_id >= 0 & d0$t_ext != 0), v[1], v[2], v[3], v[4],
              isTRUE(all.equal(v, rep(v[1], 4), tolerance = 1e-12))))
}

cat("\n==== Part E: two ts conventions inside loglik (D + exponential) ====\n")
## In loglik(), ep_exp path: sum_exp_bE starts at tree[0].n (= 2 crown lineages at
## ts = 0) and adds exp(-bD * brts) per observed node -> observed daughters born at
## brts, parents never reset.  The M used in the same integral is pd/n with pd
## counting the daughters since t = 0 and the crown lineages not at all.
## Show that the exp path is NOT flat in beta_D even though D_pkg = 0 at events,
## and differs from the truth (Part C already tabulates it).
cat("see Part C table: exp-link diff varies with beta_D ->",
    sprintf("%s", paste(round(res$diff[res$link == "exp"], 3), collapse = " ")), "\n")
cat("DONE\n")

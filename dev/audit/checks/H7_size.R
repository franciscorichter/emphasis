## H7 (b) — does the envelope defect persist on a larger tree and under dd?
## (i) 20-tip CR tree: first-missing-event survival vs closed form; mean #missing C++ vs
##     dominating-envelope replica; fhat gap C++ vs replica (both scored by eval_logf).
## (ii) 7-tip tree, dd (linear, beta_N < 0): P(no missing) C++ vs closed form.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({ library(emphasis); library(ape) })
set.seed(11)
T_EXT_TIP <- 1e11; T_EXT_EXT <- 0
is_missing <- function(df) !(df$t_ext == T_EXT_EXT | df$t_ext == T_EXT_TIP | df$t_ext == 5e10)

run_case <- function(brts_age, lam, mu, bN = 0, Nsamp = 4000L, Nrep = 2000L, label = "") {
  T <- brts_age[1]; s <- T - brts_age[-1]
  lamN <- function(n) pmax(0, lam + bN * n)
  seg_int <- function(n, a, b) n * lamN(n) * ((b - a) - (1/mu) * (exp(-mu*(T-b)) - exp(-mu*(T-a))))
  H <- function(t) { knots <- c(0, s, T); nvec <- 2 + seq_along(knots) - 1; h <- 0
    for (i in seq_len(length(knots) - 1)) { a <- knots[i]; b <- min(knots[i+1], t)
      if (b > a) h <- h + seg_int(nvec[i], a, b); if (knots[i+1] >= t) break }; h }
  p8 <- c(lam, bN, 0, 0, mu, 0, 0, 0)
  raw <- emphasis:::augment_trees(brts = brts_age, pars = p8, sample_size = Nsamp, maxN = 100L*Nsamp,
                                  max_missing = 500L, max_lambda = 1e6, num_threads = 1L,
                                  model = c(0L, as.integer(bN != 0), 0L), link = 0L, rho = 1)
  nmiss <- sapply(raw$trees, function(df) sum(is_missing(df)))
  t1 <- sapply(raw$trees, function(df) { m <- is_missing(df); if (any(m)) min(df$brts[m]) else Inf })
  cat(sprintf("\n== %s: %d tips, T = %.2f, lambda = %.2f, mu = %.2f, beta_N = %.3f; rejected: overrun %d lambda %d zero %d\n",
              label, length(brts_age) + 1, T, lam, mu, bN, raw$rejected_overruns, raw$rejected_lambda, raw$rejected_zero_weights))
  cat(sprintf("   P(no missing): closed %.4f | C++ %.4f (se %.4f)\n", exp(-H(T)), mean(nmiss == 0), sqrt(mean(nmiss==0)*(1-mean(nmiss==0))/Nsamp)))
  qs <- quantile(s, c(0.25, 0.5, 0.75, 1))
  for (tt in qs) { emp <- mean(t1 > tt); cl <- exp(-H(tt))
    cat(sprintf("   P(t1 > %5.2f): closed %.4f | C++ %.4f (se %.4f)  z = %5.1f\n", tt, cl, emp, sqrt(emp*(1-emp)/Nsamp), (emp-cl)/sqrt(emp*(1-emp)/Nsamp))) }
  ## dominating-envelope replica (time-only sampler; CR/dd, rho = 1)
  recount <- function(tr) { tr <- tr[order(tr$brts), ]; n <- numeric(nrow(tr)); n[1] <- 2
    if (nrow(tr) > 1) for (i in 2:nrow(tr)) n[i] <- n[i-1] + ifelse(tr$t_ext[i-1] == T_EXT_EXT, -1, 1); tr$n <- n; tr }
  lower_n <- function(tr, t) { i <- which(tr$brts >= t); if (length(i) == 0) nrow(tr) else i[1] }
  nh <- function(tr, t) { n <- tr$n[lower_n(tr, t)]; n * lamN(n) * (1 - exp(-mu * (T - t))) }
  next_bt <- function(tr, cbt) { i <- which(tr$brts > cbt); if (length(i)) tr$brts[i[1]] else T }
  rtrunc_exp <- function(upper) { x <- rexp(1, mu); while (x > upper) x <- rexp(1, mu); x }
  replica <- function() {
    tr <- data.frame(brts = c(s, T), n = 2 + 0:length(s), t_ext = T_EXT_TIP); cbt <- 0
    while (cbt < T) { nb <- next_bt(tr, cbt); n <- tr$n[lower_n(tr, nb)]
      lmax <- n * lamN(n) * (1 - exp(-mu * (T - cbt)))
      tstar <- if (lmax > 0) cbt - log(runif(1)) / lmax else nb
      if (tstar < nb) { pt <- nh(tr, tstar) / lmax; stopifnot(pt <= 1 + 1e-12)
        if (runif(1) < pt) { text <- tstar + rtrunc_exp(T - tstar)
          tr <- recount(rbind(tr, data.frame(brts = c(tstar, text), n = 0, t_ext = c(text, T_EXT_EXT)))) } }
      cbt <- min(tstar, nb) }
    tr$pd <- 0; tr$tip_start <- 0; tr$focal_tip_start <- 0; tr$id <- -1L; tr$parent_id <- -1L; rownames(tr) <- NULL; tr }
  dfs <- replicate(Nrep, replica(), simplify = FALSE)
  nm_fix <- sapply(dfs, function(d) sum(d$t_ext == T_EXT_EXT))
  cat(sprintf("   mean #missing: C++ %.3f (se %.3f) | dominating replica %.3f (se %.3f)\n",
              mean(nmiss), sd(nmiss)/sqrt(Nsamp), mean(nm_fix), sd(nm_fix)/sqrt(Nrep)))
  fh <- function(logf, logg, nz = 0) { lw <- logf - logg; m <- max(lw); w <- exp(lw - m)
    c(fhat = log(mean(w)) + m - log(1 + nz/length(w)), se = sd(w)/(mean(w)*sqrt(length(w))), ess = sum(w)^2/sum(w^2)) }
  ev <- emphasis:::eval_logf(p8, dfs, model = c(0L, as.integer(bN != 0), 0L), link = 0L, rho = 1)
  fc <- fh(raw$logf, raw$logg, raw$rejected_zero_weights); ff <- fh(ev$logf, ev$logg)
  cat(sprintf("   fhat: C++ %.4f (se %.4f, ESS %.0f) | dominating replica %.4f (se %.4f, ESS %.0f) | diff %.4f\n",
              fc["fhat"], fc["se"], fc["ess"], ff["fhat"], ff["se"], ff["ess"], fc["fhat"] - ff["fhat"]))
  invisible(NULL)
}

## (i) 20-tip CR tree
phy <- ape::rphylo(20, 0.4, 0.2)
ba  <- sort(as.numeric(ape::branching.times(phy)), decreasing = TRUE)
run_case(ba, 0.4, 0.2, label = "CR 20 tips")
run_case(ba, 0.4, 0.05, label = "CR 20 tips, small mu")
## (ii) 7-tip tree, dd linear
run_case(c(6, 4.5, 3.0, 2.0, 1.2, 0.5), 0.6, 0.2, bN = -0.04, label = "DD 7 tips")

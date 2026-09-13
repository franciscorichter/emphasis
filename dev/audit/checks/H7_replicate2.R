## H7 replication, part 2: does the bias survive on larger / random trees?
##   (D) 15-tip tree, N=20000: first-event-time KS + survival z (exact closed form, no replica)
##   (E) 30-tip tree: same test; mean #missing C++ vs dominating-envelope replica; fhat vs DDD
##   (F) 12-tip tree, mu-grid at fixed lambda, 3 x 4000 trees: C++ fhat - DDD (theta dependence)
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({ library(emphasis); library(DDD); library(ape) })
T_EXT_TIP <- 1e11; T_EXT_EXT <- 0
is_missing <- function(df) !(df$t_ext == T_EXT_EXT | df$t_ext == T_EXT_TIP | df$t_ext == 5e10)
closed <- function(brts_age, lam, mu) {
  T <- brts_age[1]; s <- T - brts_age[-1]
  seg_int <- function(n, a, b) n * lam * ((b - a) - (1/mu) * (exp(-mu*(T-b)) - exp(-mu*(T-a))))
  H <- function(t) { knots <- c(0, s, T); nvec <- 2 + seq_along(knots) - 1; h <- 0
    for (i in seq_len(length(knots) - 1)) { a <- knots[i]; b <- min(knots[i+1], t)
      if (b > a) h <- h + seg_int(nvec[i], a, b); if (knots[i+1] >= t) break }
    h }
  list(H = H, s = s, T = T)
}
run_cpp <- function(brts_age, pars8, N) emphasis:::augment_trees(brts = brts_age, pars = pars8, sample_size = N,
  maxN = 50L*N, max_missing = 400L, max_lambda = 1e6, num_threads = 1L, model = c(0L,0L,0L), link = 0L, rho = 1)
first_event_test <- function(tag, raw, cf, N, at) {
  nmiss <- sapply(raw$trees, function(df) sum(is_missing(df)))
  t1 <- sapply(raw$trees, function(df) { m <- is_missing(df); if (any(m)) min(df$brts[m]) else Inf })
  obs <- t1[is.finite(t1)]
  Ftr <- function(t) (1 - exp(-sapply(t, cf$H))) / (1 - exp(-cf$H(cf$T)))
  ks <- suppressWarnings(ks.test(obs, Ftr))
  cat(sprintf("  %s: N=%d, mean #missing %.3f, first-event KS D=%.4f p=%.1e; rej ovr/lam/zero=%d/%d/%d\n", tag, N, mean(nmiss),
              ks$statistic, ks$p.value, raw$rejected_overruns, raw$rejected_lambda, raw$rejected_zero_weights))
  for (tt in at) { emp <- mean(t1 > tt); cl <- exp(-cf$H(tt))
    cat(sprintf("     P(t1 > %5.2f): closed %.4f | C++ %.4f (se %.4f) z=%+.1f\n", tt, cl, emp, sqrt(emp*(1-emp)/N), (emp-cl)/sqrt(emp*(1-emp)/N))) }
  nmiss
}
## dominating-envelope replica (n of the segment (cbt,nb] x survival factor at cbt)
replica_dom <- function(s, T, lam, mu) {
  b <- c(s, T); e <- rep(T_EXT_TIP, length(b))
  recount <- function() { n <- numeric(length(b)); n[1] <- 2
    if (length(b) > 1) for (i in 2:length(b)) n[i] <- n[i-1] + ifelse(e[i-1] == T_EXT_EXT, -1, 1); n }
  n <- recount(); surv <- function(t) 1 - exp(-mu*(T-t))
  n_seg <- function(t) { i <- which(b > t); i <- if (length(i)) i[1] else length(b); n[i] }
  cbt <- 0
  while (cbt < T) {
    nb <- { i <- which(b > cbt); if (length(i)) b[i[1]] else T }
    lmax <- n_seg(cbt) * lam * surv(cbt)
    tstar <- if (lmax > 0) cbt - log(runif(1))/lmax else nb
    if (tstar < nb) {
      pt <- n_seg(cbt) * lam * surv(tstar) / lmax
      if (runif(1) < pt) { x <- rexp(1, mu); while (x > T - tstar) x <- rexp(1, mu); te <- tstar + x
        b <- c(b, tstar, te); e <- c(e, te, T_EXT_EXT); o <- order(b); b <- b[o]; e <- e[o]; n <- recount() }
    }
    cbt <- min(tstar, nb)
  }
  data.frame(brts = b, n = n, t_ext = e, pd = 0, tip_start = 0, focal_tip_start = 0, id = -1L, parent_id = -1L)
}
fhat_se <- function(logf, logg, nzero = 0) { lw <- logf - logg; m <- max(lw); w <- exp(lw - m)
  c(fhat = log(mean(w)) + m - log(1 + nzero/length(w)), se = sd(w)/(mean(w)*sqrt(length(w))), ess = sum(w)^2/sum(w^2)) }

## ---- (D) -------------------------------------------------------------------
set.seed(11); brtsA <- sort(as.numeric(branching.times(rphylo(15, 0.5, 0.3))), decreasing = TRUE)
lam <- 0.5; mu <- 0.3; cf <- closed(brtsA, lam, mu)
cat(sprintf("\n[D] 15-tip tree (crown %.2f), lambda=%.2f mu=%.2f\n", cf$T, lam, mu))
first_event_test("C++", run_cpp(brtsA, c(lam,0,0,0,mu,0,0,0), 20000L), cf, 20000L, quantile(cf$s, c(.25,.5,.75,1)))

## ---- (E) -------------------------------------------------------------------
set.seed(5); brtsE <- sort(as.numeric(branching.times(rphylo(30, 0.4, 0.2))), decreasing = TRUE)
lam <- 0.4; mu <- 0.2; cfE <- closed(brtsE, lam, mu)
cat(sprintf("\n[E] 30-tip tree (crown %.2f), lambda=%.2f mu=%.2f\n", cfE$T, lam, mu))
rawE <- run_cpp(brtsE, c(lam,0,0,0,mu,0,0,0), 10000L)
nmE <- first_event_test("C++", rawE, cfE, 10000L, quantile(cfE$s, c(.25,.5,.75,1)))
set.seed(31); Ndom <- 800L
dom <- replicate(Ndom, replica_dom(cfE$s, cfE$T, lam, mu), simplify = FALSE)
nm_dom <- sapply(dom, function(d) sum(d$t_ext == T_EXT_EXT))
cat(sprintf("     mean #missing: C++ %.3f (se %.3f) | dominating replica %.3f (se %.3f) | z=%+.1f\n",
            mean(nmE), sd(nmE)/sqrt(length(nmE)), mean(nm_dom), sd(nm_dom)/sqrt(Ndom),
            (mean(nmE)-mean(nm_dom))/sqrt(var(nmE)/length(nmE)+var(nm_dom)/Ndom)))
ddd <- DDD::bd_loglik(pars1 = c(lam, mu, 0, 0), pars2 = c(0, 0, 1, 0, 2), brts = brtsE, missnumspec = 0)
fc <- fhat_se(rawE$logf, rawE$logg, rawE$rejected_zero_weights)
ev <- emphasis:::eval_logf(c(lam,0,0,0,mu,0,0,0), dom, model = c(0L,0L,0L), link = 0L, rho = 1)
fd <- fhat_se(ev$logf, ev$logg)
cat(sprintf("     fhat - DDD: C++ %+.4f (se %.4f, ESS %.0f) | dominating replica %+.4f (se %.4f, ESS %.0f)\n",
            fc["fhat"]-ddd, fc["se"], fc["ess"], fd["fhat"]-ddd, fd["se"], fd["ess"]))

## ---- (F) -------------------------------------------------------------------
set.seed(3); brtsF <- sort(as.numeric(branching.times(rphylo(12, 0.4, 0.2))), decreasing = TRUE)
cat(sprintf("\n[F] 12-tip tree (crown %.2f): C++ fhat - DDD, 3 runs x 4000 trees, lambda fixed 0.4\n", brtsF[1]))
for (mu in c(0.05, 0.2, 0.35)) {
  lam <- 0.4; ddd <- DDD::bd_loglik(pars1 = c(lam, mu, 0, 0), pars2 = c(0, 0, 1, 0, 2), brts = brtsF, missnumspec = 0)
  g <- sapply(1:3, function(i) { r <- run_cpp(brtsF, c(lam,0,0,0,mu,0,0,0), 4000L); f <- fhat_se(r$logf, r$logg, r$rejected_zero_weights); c(f["fhat"]-ddd, f["se"], f["ess"]) })
  cat(sprintf("    mu=%.2f: gaps %s | pooled %+.4f (se %.4f) | ESS %s\n", mu, paste(sprintf("%+.4f", g[1,]), collapse=" "),
              mean(g[1,]), sqrt(sum(g[2,]^2))/3, paste(sprintf("%.0f", g[3,]), collapse=" ")))
}

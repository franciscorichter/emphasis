## H7 replication (independent). Varies what the verifier did not:
##   (A) a different tree (15-tip ape::rphylo), lambda=.5 mu=.3, linear AND exponential link:
##       P(no missing) closed form vs C++ thinning sampler; survival of first event.
##   (B) mu-grid at fixed lambda on a 20-tip tree: C++ fhat - DDD::bd_loglik (exact reference),
##       is the gap mu-dependent beyond MC error?
##   (C) lean R replica of do_augment_tree_cont on the 7-tip tree with the two defects
##       switched on/off separately (pre-event n at interval start; envelope collapse
##       after rejection) to see which one carries the bias.
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

run_cpp <- function(brts_age, pars8, link, N, model = c(0L,0L,0L)) {
  emphasis:::augment_trees(brts = brts_age, pars = pars8, sample_size = N, maxN = 50L*N,
                           max_missing = 300L, max_lambda = 1e6, num_threads = 1L,
                           model = model, link = link, rho = 1)
}
report <- function(tag, raw, cf, N) {
  nmiss <- sapply(raw$trees, function(df) sum(is_missing(df)))
  t1 <- sapply(raw$trees, function(df) { m <- is_missing(df); if (any(m)) min(df$brts[m]) else Inf })
  p0 <- mean(nmiss == 0); se <- sqrt(p0*(1-p0)/N); p0c <- exp(-cf$H(cf$T))
  cat(sprintf("  %-10s P(0): closed %.4f | C++ %.4f (se %.4f) z=%+.1f | mean #missing %.3f | rej overrun/lambda/zero = %d/%d/%d\n",
              tag, p0c, p0, se, (p0-p0c)/se, mean(nmiss),
              raw$rejected_overruns, raw$rejected_lambda, raw$rejected_zero_weights))
  zs <- sapply(cf$s, function(tt) { emp <- mean(t1 > tt); cl <- exp(-cf$H(tt)); (emp-cl)/sqrt(emp*(1-emp)/N) })
  cat(sprintf("  %-10s survival z at each branching time: %s\n", "", paste(sprintf("%+.1f", zs), collapse = " ")))
  invisible(nmiss)
}

## ---- (A) 15-tip tree, both links ------------------------------------------
set.seed(11)
phy <- ape::rphylo(15, 0.5, 0.3)
brtsA <- sort(as.numeric(ape::branching.times(phy)), decreasing = TRUE)
lam <- 0.5; mu <- 0.3; N <- 6000L
cat(sprintf("\n[A] 15-tip rphylo tree, crown age %.2f, lambda=%.2f mu=%.2f, N=%d\n", brtsA[1], lam, mu, N))
cf <- closed(brtsA, lam, mu)
rawL <- run_cpp(brtsA, c(lam,0,0,0,mu,0,0,0), 0L, N)
report("linear", rawL, cf, N)
rawE <- run_cpp(brtsA, c(log(lam),0,0,0,log(mu),0,0,0), 1L, N)
report("exp-link", rawE, cf, N)
## sanity: logg of empty augmentations equals -H(T) for both links
for (r in list(rawL, rawE)) {
  nm <- sapply(r$trees, function(df) sum(is_missing(df)))
  if (any(nm == 0)) cat(sprintf("     logg(empty) = %.5f vs -H(T) = %.5f\n", r$logg[nm == 0][1], -cf$H(cf$T)))
}

## ---- (B) mu-grid at fixed lambda, 20-tip tree, C++ fhat vs DDD -------------
set.seed(23)
phy20 <- ape::rphylo(20, 0.4, 0.2)
brtsB <- sort(as.numeric(ape::branching.times(phy20)), decreasing = TRUE)
fhat_se <- function(logf, logg, nzero = 0) {
  lw <- logf - logg; m <- max(lw); w <- exp(lw - m)
  c(fhat = log(mean(w)) + m - log(1 + nzero/length(w)), se = sd(w)/(mean(w)*sqrt(length(w))), ess = sum(w)^2/sum(w^2))
}
cat(sprintf("\n[B] 20-tip rphylo tree, crown age %.2f; C++ fhat - DDD::bd_loglik(cond=0,btorph=1,soc=2), 2 independent runs of 3000 trees\n", brtsB[1]))
cat(sprintf("    %6s %6s | %9s %6s %5s | %9s %6s %5s\n", "lambda", "mu", "gap1", "se", "ESS", "gap2", "se", "ESS"))
gridB <- rbind(c(0.4, 0.02), c(0.4, 0.1), c(0.4, 0.2), c(0.4, 0.35), c(0.25, 0.2), c(0.6, 0.2))
resB <- NULL
for (k in seq_len(nrow(gridB))) {
  lam <- gridB[k,1]; mu <- gridB[k,2]
  ddd <- DDD::bd_loglik(pars1 = c(lam, mu, 0, 0), pars2 = c(0, 0, 1, 0, 2), brts = brtsB, missnumspec = 0)
  g <- sapply(1:2, function(i) { r <- run_cpp(brtsB, c(lam,0,0,0,mu,0,0,0), 0L, 3000L)
    f <- fhat_se(r$logf, r$logg, r$rejected_zero_weights); c(f["fhat"] - ddd, f["se"], f["ess"]) })
  cat(sprintf("    %6.2f %6.2f | %9.4f %6.4f %5.0f | %9.4f %6.4f %5.0f\n", lam, mu, g[1,1], g[2,1], g[3,1], g[1,2], g[2,2], g[3,2]))
  resB <- rbind(resB, c(lam, mu, mean(g[1,]), sqrt(sum(g[2,]^2))/2))
}
cat(sprintf("    pooled gap at mu=0.02: %+.4f (se %.4f); at mu=0.35: %+.4f (se %.4f); difference z = %.1f\n",
            resB[1,3], resB[1,4], resB[4,3], resB[4,4], (resB[4,3]-resB[1,3])/sqrt(resB[1,4]^2+resB[4,4]^2)))

## ---- (C) lean replica, defects separately, 7-tip tree -----------------------
brts7 <- c(6, 4.5, 3.0, 2.0, 1.2, 0.5); lam <- 0.4; mu <- 0.2
cf7 <- closed(brts7, lam, mu); T <- cf7$T; s <- cf7$s
replica <- function(pre_event_n = TRUE, collapse = TRUE) {
  b <- c(s, T); e <- rep(T_EXT_TIP, length(b))              # nodes sorted by brts
  recount <- function() { n <- numeric(length(b)); n[1] <- 2
    if (length(b) > 1) for (i in 2:length(b)) n[i] <- n[i-1] + ifelse(e[i-1] == T_EXT_EXT, -1, 1); n }
  n <- recount()
  surv <- function(t) 1 - exp(-mu*(T-t))
  nh_at <- function(t) { i <- which(b >= t); i <- if (length(i)) i[1] else length(b); n[i] * lam * surv(t) }   # C++ lower_bound
  n_seg <- function(t) { i <- which(b > t); i <- if (length(i)) i[1] else length(b); n[i] }                   # upper_bound (post-event)
  start_rate <- function(t) if (pre_event_n) nh_at(t) else n_seg(t) * lam * surv(t)
  cbt <- 0; lambda2 <- 0; dirty <- TRUE; fresh <- TRUE; lmax <- 0; gt1 <- 0L; nc <- 0L
  while (cbt < T) {
    nb <- { i <- which(b > cbt); if (length(i)) b[i[1]] else T }
    if (collapse) {                       # C++ semantics: after a rejection lambda1 <- lambda2
      lambda1 <- if (!dirty && !fresh) lambda2 else start_rate(cbt)
      lambda2 <- nh_at(nb); lmax <- max(lambda1, lambda2)
    } else if (fresh) {                   # keep the interval envelope across rejections
      lambda1 <- start_rate(cbt); lambda2 <- nh_at(nb); lmax <- max(lambda1, lambda2)
    }
    tstar <- if (lmax > 0) cbt - log(runif(1))/lmax else nb
    dirty <- FALSE; fresh <- FALSE
    if (tstar < nb) {
      nc <- nc + 1L; pt <- nh_at(tstar)/lmax; if (pt > 1) gt1 <- gt1 + 1L
      if (runif(1) < pt) {
        x <- rexp(1, mu); while (x > T - tstar) x <- rexp(1, mu); te <- tstar + x
        b <- c(b, tstar, te); e <- c(e, te, T_EXT_EXT); o <- order(b); b <- b[o]; e <- e[o]; n <- recount()
        dirty <- TRUE; fresh <- TRUE
      }
    } else fresh <- TRUE
    cbt <- min(tstar, nb)
  }
  c(nmiss = sum(e == T_EXT_EXT), gt1 = gt1, nc = nc)
}
set.seed(99); Nrep <- 3000L
cat(sprintf("\n[C] 7-tip tree lambda=.4 mu=.2, lean replica, closed P(0) = %.4f, Nrep=%d\n", exp(-cf7$H(T)), Nrep))
for (cfg in list(c(TRUE,TRUE), c(FALSE,TRUE), c(TRUE,FALSE), c(FALSE,FALSE))) {
  r <- replicate(Nrep, replica(cfg[1], cfg[2]))
  p0 <- mean(r["nmiss",] == 0)
  cat(sprintf("  pre_event_n=%-5s collapse=%-5s | pt>1: %5.1f%% | P(0) %.4f (se %.4f) | mean #missing %.3f\n",
              cfg[1], cfg[2], 100*sum(r["gt1",])/sum(r["nc",]), p0, sqrt(p0*(1-p0)/Nrep), mean(r["nmiss",])))
}

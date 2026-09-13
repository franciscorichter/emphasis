## H7 — thinning envelope in do_augment_tree_cont is not a dominating rate.
## Self-contained. Tests, on a fixed 7-tip CR tree (rho = 1):
##   (1) P(no missing lineage) from the installed C++ sampler vs the closed form
##       exp(-int nh) that sampling_prob charges (logg of an empty augmentation);
##   (2) the first-missing-event time distribution vs 1 - exp(-H(t)) (KS);
##   (3) an R replica of the C++ loop (same lower_bound semantics) that counts
##       pt > 1, reproduces the C++ P(0), and a replica with a dominating
##       envelope that reproduces the closed form;
##   (4) fhat(theta) from the C++ sampler vs fhat from the dominating-envelope
##       replica (both scored by emphasis:::eval_logf), against DDD::bd_loglik
##       across a theta grid — is the H7 gap theta-dependent?
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({ library(emphasis); library(DDD) })
set.seed(7)

## ---- tree, parameters -------------------------------------------------------
brts_age <- c(6, 4.5, 3.0, 2.0, 1.2, 0.5)          # 7 tips, crown age 6
T  <- brts_age[1]
s  <- T - brts_age[-1]                             # forward branching times
lam <- 0.4; mu <- 0.2
pars8 <- c(lam, 0, 0, 0, mu, 0, 0, 0)              # linear link: lambda = beta0, mu = gamma0
Nsamp <- 5000L
T_EXT_TIP <- 1e11; T_EXT_EXT <- 0

## ---- closed forms ----------------------------------------------------------
## n(u) on the OBSERVED tree = 2 + #{s_i < u}; nh(u) = n(u) lam (1 - exp(-mu (T-u)))
seg_int <- function(n, a, b) n * lam * ((b - a) - (1/mu) * (exp(-mu*(T-b)) - exp(-mu*(T-a))))
H <- function(t) {                                   # cumulative hazard on [0,t]
  knots <- c(0, s, T); nvec <- 2 + seq_along(knots) - 1
  h <- 0
  for (i in seq_len(length(knots) - 1)) {
    a <- knots[i]; b <- min(knots[i+1], t)
    if (b > a) h <- h + seg_int(nvec[i], a, b)
    if (knots[i+1] >= t) break
  }
  h
}
p0_closed <- exp(-H(T))

## ---- (1)+(2) installed C++ sampler ------------------------------------------
raw <- emphasis:::augment_trees(brts = brts_age, pars = pars8, sample_size = Nsamp,
                                maxN = 50L * Nsamp, max_missing = 200L, max_lambda = 1e6,
                                num_threads = 1L, model = c(0L,0L,0L), link = 0L, rho = 1)
is_missing <- function(df) !(df$t_ext == T_EXT_EXT | df$t_ext == T_EXT_TIP | df$t_ext == 5e10)
nmiss <- sapply(raw$trees, function(df) sum(is_missing(df)))
t1    <- sapply(raw$trees, function(df) { m <- is_missing(df); if (any(m)) min(df$brts[m]) else Inf })
p0_cpp <- mean(nmiss == 0)
se0 <- sqrt(p0_cpp * (1 - p0_cpp) / Nsamp)
cat(sprintf("\n[1] P(no missing): closed form exp(-int nh) = %.4f | C++ empirical = %.4f (se %.4f) | z = %.1f\n",
            p0_closed, p0_cpp, se0, (p0_cpp - p0_closed)/se0))
## logg charged for the empty augmentation must equal -H(T): confirms closed form == sampling_prob
lg0 <- raw$logg[nmiss == 0]
cat(sprintf("    logg of empty augmentations: %.6f (all equal: %s) | -H(T) = %.6f\n",
            lg0[1], all(abs(lg0 - lg0[1]) < 1e-9), -H(T)))

## first-event time: KS against truncated CDF F(t)/F(T), F = 1 - exp(-H)
obs <- t1[is.finite(t1)]
Ftr <- function(t) (1 - exp(-sapply(t, H))) / (1 - exp(-H(T)))
ks  <- suppressWarnings(ks.test(obs, Ftr))
cat(sprintf("[2] first missing event time: n = %d, KS D = %.4f, p = %.2e\n", length(obs), ks$statistic, ks$p.value))
## where does the mismatch sit? survival P(t1 > t) at each observed branching time
for (tt in s) {
  emp <- mean(t1 > tt); cl <- exp(-H(tt))
  cat(sprintf("    P(t1 > %.2f): closed %.4f | C++ %.4f (se %.4f)\n", tt, cl, emp, sqrt(emp*(1-emp)/Nsamp)))
}

## ---- (3) R replica of do_augment_tree_cont for CR, rho = 1 ------------------
## node list: brts, n (count on the segment ENDING at the node), t_ext
init_tree <- function() data.frame(brts = c(s, T), n = 2 + 0:length(s), t_ext = T_EXT_TIP)
recount <- function(tr) {                              # n_after chain, as insert_species maintains it
  tr <- tr[order(tr$brts), ]
  n <- numeric(nrow(tr)); n[1] <- 2
  if (nrow(tr) > 1) for (i in 2:nrow(tr)) n[i] <- n[i-1] + ifelse(tr$t_ext[i-1] == T_EXT_EXT, -1, 1)
  tr$n <- n; tr
}
lower_n <- function(tr, t) {                           # C++ lower_bound_node: first node with brts >= t, capped
  i <- which(tr$brts >= t); if (length(i) == 0) nrow(tr) else i[1]
}
nh <- function(tr, t) tr$n[lower_n(tr, t)] * lam * (1 - exp(-mu * (T - t)))
next_bt <- function(tr, cbt) { i <- which(tr$brts > cbt); if (length(i)) tr$brts[i[1]] else T }
rtrunc_exp <- function(upper) { x <- rexp(1, mu); while (x > upper) x <- rexp(1, mu); x }

replica <- function(fixed = FALSE) {
  tr <- init_tree(); cbt <- 0; lambda2 <- 0; dirty <- TRUE; n_gt1 <- 0L; n_cand <- 0L
  while (cbt < T) {
    nb <- next_bt(tr, cbt)
    if (!fixed) {
      lambda1 <- if (!dirty) lambda2 else max(0, nh(tr, cbt))
      lambda2 <- max(0, nh(tr, nb))
      lmax <- max(lambda1, lambda2)
    } else {
      ## dominating envelope: n of the segment (cbt, nb] times the survival factor at cbt
      lmax <- tr$n[lower_n(tr, nb)] * lam * (1 - exp(-mu * (T - cbt)))
    }
    tstar <- if (lmax > 0) cbt - log(runif(1)) / lmax else nb
    dirty <- FALSE
    if (tstar < nb) {
      n_cand <- n_cand + 1L
      pt <- max(0, nh(tr, tstar)) / lmax
      if (pt > 1) n_gt1 <- n_gt1 + 1L
      if (runif(1) < pt) {
        text <- tstar + rtrunc_exp(T - tstar)
        tr <- rbind(tr, data.frame(brts = c(tstar, text), n = 0, t_ext = c(text, T_EXT_EXT)))
        tr <- recount(tr); dirty <- TRUE
      }
    }
    cbt <- min(tstar, nb)
  }
  list(tree = tr, n_gt1 = n_gt1, n_cand = n_cand)
}
Nrep <- 3000L
rep_bug <- replicate(Nrep, replica(FALSE), simplify = FALSE)
rep_fix <- replicate(Nrep, replica(TRUE),  simplify = FALSE)
nm_bug <- sapply(rep_bug, function(r) sum(r$tree$t_ext == T_EXT_EXT))
nm_fix <- sapply(rep_fix, function(r) sum(r$tree$t_ext == T_EXT_EXT))
gt1 <- sum(sapply(rep_bug, `[[`, "n_gt1")); cand <- sum(sapply(rep_bug, `[[`, "n_cand"))
gt1f <- sum(sapply(rep_fix, `[[`, "n_gt1"))
cat(sprintf("\n[3] replica (C++ semantics): candidates with pt > 1: %d of %d (%.1f%%); dominating replica: %d\n",
            gt1, cand, 100*gt1/cand, gt1f))
cat(sprintf("    P(no missing): C++ %.4f | buggy replica %.4f (se %.4f) | dominating replica %.4f (se %.4f) | closed %.4f\n",
            p0_cpp, mean(nm_bug==0), sqrt(mean(nm_bug==0)*(1-mean(nm_bug==0))/Nrep),
            mean(nm_fix==0), sqrt(mean(nm_fix==0)*(1-mean(nm_fix==0))/Nrep), p0_closed))
cat(sprintf("    mean #missing: C++ %.3f | buggy replica %.3f | dominating replica %.3f\n",
            mean(nmiss), mean(nm_bug), mean(nm_fix)))

## ---- (4) fhat vs DDD across a theta grid -------------------------------------
to_df <- function(tr) { tr <- recount(tr); tr$pd <- 0; tr$tip_start <- 0; tr$focal_tip_start <- 0
  tr$id <- -1L; tr$parent_id <- -1L; rownames(tr) <- NULL; tr }
fhat_se <- function(logf, logg, nzero = 0) {
  lw <- logf - logg; m <- max(lw); w <- exp(lw - m)
  c(fhat = log(mean(w)) + m - log(1 + nzero/length(w)),
    se = sd(w) / (mean(w) * sqrt(length(w))),
    ess = sum(w)^2 / sum(w^2))
}
grid <- rbind(c(0.3, 0.1), c(0.4, 0.2), c(0.5, 0.3), c(0.6, 0.45), c(0.35, 0.3))
cat("\n[4] fhat(theta) - DDD::bd_loglik(cond=0, btorph=1, soc=2):  constant offset expected\n")
cat(sprintf("    %6s %6s | %10s %6s %6s | %10s %6s %6s | %8s\n", "lambda", "mu",
            "C++ gap", "se", "ESS", "dom gap", "se", "ESS", "C++-dom"))
gaps <- NULL
for (k in seq_len(nrow(grid))) {
  lam <<- grid[k,1]; mu <<- grid[k,2]; p8 <- c(lam,0,0,0,mu,0,0,0)
  rc <- emphasis:::augment_trees(brts = brts_age, pars = p8, sample_size = 4000L, maxN = 200000L,
                                 max_missing = 200L, max_lambda = 1e6, num_threads = 1L,
                                 model = c(0L,0L,0L), link = 0L, rho = 1)
  fc <- fhat_se(rc$logf, rc$logg, rc$rejected_zero_weights)
  rf <- replicate(2500L, replica(TRUE), simplify = FALSE)
  dfs <- lapply(rf, function(r) to_df(r$tree))
  ev <- emphasis:::eval_logf(p8, dfs, model = c(0L,0L,0L), link = 0L, rho = 1)
  ff <- fhat_se(ev$logf, ev$logg)
  ddd <- DDD::bd_loglik(pars1 = c(lam, mu, 0, 0), pars2 = c(0, 0, 1, 0, 2), brts = brts_age, missnumspec = 0)  # installed DDD: c(tdmodel, cond, btorph, verbose, soc)
  cat(sprintf("    %6.2f %6.2f | %10.4f %6.4f %6.0f | %10.4f %6.4f %6.0f | %8.4f\n",
              lam, mu, fc["fhat"] - ddd, fc["se"], fc["ess"], ff["fhat"] - ddd, ff["se"], ff["ess"],
              fc["fhat"] - ff["fhat"]))
  gaps <- rbind(gaps, c(lam, mu, fc["fhat"] - ddd, fc["se"], ff["fhat"] - ddd, ff["se"]))
}
cat(sprintf("    spread of C++ gap across grid: %.4f | spread of dominating-replica gap: %.4f\n",
            diff(range(gaps[,3])), diff(range(gaps[,5]))))

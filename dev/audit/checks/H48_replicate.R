## H48 independent replication: vary tree, sampler model/link, rho, max_missing.
## Self-contained; ~60 s.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({ library(emphasis); library(ape) })

T_EXT_TIP <- 1e11; T_EXT_UNS <- 5e10; T_EXT_EXT <- 0.0
is_ext  <- function(df) df$t_ext == T_EXT_EXT
is_tip  <- function(df) df$t_ext == T_EXT_TIP
is_uns  <- function(df) df$t_ext == T_EXT_UNS
is_miss <- function(df) !(is_ext(df) | is_tip(df) | is_uns(df))
is_aug  <- function(df) is_miss(df) | is_uns(df)
hand_P  <- function(df, tm, ts) { alive <- !is_ext(df) & df$brts <= tm & df$t_ext > tm; sum(tm - ts[alive]) }

## --- a different tree: simulated BD tree, 18-30 extant tips -------------------
set.seed(4821)
repeat {
  tr <- ape::rlineage(0.35, 0.15, Tmax = 8)
  tr <- ape::drop.fossil(tr)
  if (!is.null(tr) && Ntip(tr) >= 18 && Ntip(tr) <= 30) break
}
brts <- sort(as.numeric(ape::branching.times(tr)), decreasing = TRUE)
crown <- brts[1]; first_split_fwd <- crown - brts[2]
cat(sprintf("tree: %d tips, crown age %.3f, first non-crown split at forward time %.3f\n",
            Ntip(tr), crown, first_split_fwd))

invariants <- function(trees, label) {
  aug <- do.call(rbind, lapply(seq_along(trees), function(k) {
    d <- trees[[k]]; a <- which(is_aug(d)); if (!length(a)) return(NULL)
    data.frame(draw = k, brts = d$brts[a], parent_id = d$parent_id[a], tip_start = d$tip_start[a],
               unsampled = is_uns(d)[a])
  }))
  m1 <- aug$parent_id == -1
  cat(sprintf("[%s] %d draws, %d augmented lineages, %d (%.1f%%) parent_id == -1; %d unsampled-extant among -1\n",
              label, length(trees), nrow(aug), sum(m1), 100 * mean(m1), sum(aug$unsampled[m1])))
  cat(sprintf("   parent -1 => tip_start == 0: %s ; parent >=0 => tip_start == brts: %s ; all -1 born before first split (%.3f): %s\n",
              all(aug$tip_start[m1] == 0), all(abs(aug$tip_start - aug$brts)[!m1] < 1e-12),
              first_split_fwd, all(aug$brts[m1] < first_split_fwd + 1e-12)))
  # pd reproduced by stored tip_start; extinction node keeps t_spec; over-count == t_spec
  chk <- t(sapply(which(sapply(trees, function(d) any(is_aug(d) & d$parent_id == -1))), function(k) {
    d <- trees[[k]]
    P_code <- sapply(d$brts, function(tm) hand_P(d, tm, d$tip_start))
    ts_true <- d$tip_start; ts_true[is_aug(d)] <- d$brts[is_aug(d)]
    P_true <- sapply(d$brts, function(tm) hand_P(d, tm, ts_true))
    # extinction nodes of parent -1 lineages: tip_start == t_spec ?
    a <- which(is_miss(d) & d$parent_id == -1)
    e_ok <- all(sapply(a, function(i) { e <- which(is_ext(d) & d$id == d$id[i]); length(e) == 1 && abs(d$tip_start[e] - d$brts[i]) < 1e-12 }))
    c(pd_vs_code = max(abs(d$pd - P_code)), pd_over_max = max(d$pd - P_true), pd_over_min = min(d$pd - P_true),
      sum_tspec_m1 = sum(d$brts[is_aug(d) & d$parent_id == -1]), ext_ok = e_ok)
  }))
  cat(sprintf("   pd == formula(stored tip_start): max|diff| %.2g ; pd - P_true in [%.2g, %.3f]; max over-count <= sum(t_spec of -1 lineages alive): %s ; extinction nodes keep t_spec: %s\n",
              max(chk[, "pd_vs_code"]), min(chk[, "pd_over_min"]), max(chk[, "pd_over_max"]),
              all(chk[, "pd_over_max"] <= chk[, "sum_tspec_m1"] + 1e-9), all(chk[, "ext_ok"] == 1)))
  invisible(aug)
}

## --- 1. sampler run under the user-facing nd model, exp link (pd enters nh_rate) ----
p_nd_exp <- c(log(0.3), -0.01, 0, 0.2, log(0.15), 0, 0, 0.1)
raw1 <- emphasis:::augment_trees(brts, p_nd_exp, sample_size = 300L, maxN = 100000L, max_missing = 50L,
                                 max_lambda = 500, num_threads = 1L, model = c(1L,0L,1L), link = 1L, rho = 1.0)
invisible(invariants(raw1$trees, "nd/exp, rho=1"))
ev1 <- emphasis:::eval_logf(p_nd_exp, raw1$trees, model = c(1L,0L,1L), link = 1L, rho = 1.0)
cat(sprintf("   eval_logf reproduces E-step for nd/exp: max|dlogf| %.2g, max|dlogg| %.2g\n",
            max(abs(ev1$logf - raw1$logf)), max(abs(ev1$logg - raw1$logg))))

## --- 2. gaussian link, rho = 0.8 (unsampled extant path also gets parent -1) -----------
p_nd_g <- c(0.3, 0.05, 0, 0.2, 0.15, 0.05, 0, 0.1)
raw2 <- emphasis:::augment_trees(brts, p_nd_g, sample_size = 300L, maxN = 100000L, max_missing = 50L,
                                 max_lambda = 500, num_threads = 1L, model = c(1L,0L,1L), link = 2L, rho = 0.8)
invisible(invariants(raw2$trees, "nd/gauss, rho=0.8"))

## --- 3. linear link, nd, max_missing large: is dlogf still 0 as with max_missing = 1? ---
p_nd_lin <- c(0.3, -0.005, 0, 0.1, 0.15, 0, 0, 0.05)
raw3 <- emphasis:::augment_trees(brts, p_nd_lin, sample_size = 300L, maxN = 100000L, max_missing = 50L,
                                 max_lambda = 500, num_threads = 1L, model = c(1L,0L,1L), link = 0L, rho = 1.0)
invisible(invariants(raw3$trees, "nd/linear, rho=1"))
fix_both <- function(df) { a <- which(is_aug(df)); df$tip_start[a] <- df$brts[a]
                           df$pd <- sapply(df$brts, function(tm) hand_P(df, tm, df$tip_start)); df }
score <- function(df, pars, model, link, rho = 1.0) {
  r <- emphasis:::eval_logf(pars, list(df), model = as.integer(model), link = as.integer(link), rho = rho)
  c(logf = r$logf, logg = r$logg)
}
delta <- function(trees, pars, model, link, rho = 1.0) {
  sel <- which(sapply(trees, function(d) any(is_aug(d) & d$parent_id == -1)))
  d <- t(sapply(trees[sel], function(df) score(df, pars, model, link, rho) - score(fix_both(df), pars, model, link, rho)))
  sprintf("n=%d draws with a -1 lineage: max|dlogf| = %.4f (mean %.4f), max|dlogg| = %.4f ; -1 lineages per draw max %d",
          length(sel), max(abs(d[, "logf"])), mean(abs(d[, "logf"])), max(abs(d[, "logg"])),
          max(sapply(trees[sel], function(dd) sum(is_aug(dd) & dd$parent_id == -1))))
}
cat("nd/linear   :", delta(raw3$trees, p_nd_lin, c(1,0,1), 0), "\n")
cat("nd/exp      :", delta(raw1$trees, p_nd_exp, c(1,0,1), 1), "\n")
cat("nd/gauss r.8:", delta(raw2$trees, p_nd_g,   c(1,0,1), 2, 0.8), "\n")
# control: same trees scored under cr / dd (beta_D = gamma_D = 0) must give 0
cat("control dd/linear on nd/linear trees:", delta(raw3$trees, c(0.3, -0.005, 0, 0, 0.15, 0, 0, 0), c(1,0,0), 0), "\n")

## --- 4. prevalence depends on the tree: a tree whose first non-crown split is very early ---
brts_early <- c(6, 5.97, 4.5, 3.8, 3.1, 2.0, 1.4, 0.9, 0.4)
raw4 <- emphasis:::augment_trees(brts_early, c(0.2,0,0,0,0.5,0,0,0), sample_size = 300L, maxN = 100000L,
                                 max_missing = 50L, max_lambda = 500, num_threads = 1L,
                                 model = c(0L,0L,0L), link = 0L, rho = 1.0)
first_split_fwd <- 0.03
invisible(invariants(raw4$trees, "cr, first split at 0.03"))

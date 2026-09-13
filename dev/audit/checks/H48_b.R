## H48 (b): isolate which field carries the effect, reproduce the C++ D+exp loglik
## in R to pin the pendant-start used by the running sums vs by P(t), and check
## that eval_logf on the returned data frame reproduces the E-step's own logf/logg.
## Self-contained; ~20 s.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages(library(emphasis))

T_EXT_TIP <- 1e11; T_EXT_UNS <- 5e10; T_EXT_EXT <- 0.0
is_ext  <- function(df) df$t_ext == T_EXT_EXT
is_tip  <- function(df) df$t_ext == T_EXT_TIP
is_uns  <- function(df) df$t_ext == T_EXT_UNS
is_miss <- function(df) !(is_ext(df) | is_tip(df) | is_uns(df))
hand_P  <- function(df, tm, ts) { alive <- !is_ext(df) & df$brts <= tm & df$t_ext > tm; sum(tm - ts[alive]) }

brts  <- c(4, 3.2, 2.5, 1.8, 1.1, 0.5)
pars8 <- c(0.2, 0, 0, 0, 0.5, 0, 0, 0)
raw <- emphasis:::augment_trees(brts, pars8, sample_size = 400L, maxN = 200000L,
                                max_missing = 1L, max_lambda = 500, num_threads = 1L,
                                model = c(0L, 0L, 0L), link = 0L, rho = 1.0)
trees <- raw$trees
n_miss <- sapply(trees, function(d) sum(is_miss(d)))
pid    <- sapply(trees, function(d) { m <- which(is_miss(d)); if (length(m)) d$parent_id[m[1]] else NA })
sel_m1 <- which(n_miss == 1 & pid == -1)
cat("parent_id == -1 single-missing draws:", length(sel_m1), "\n")

## 0. eval_logf on the returned df reproduces the E-step's own logf/logg (cr linear)
ev <- emphasis:::eval_logf(pars8, trees, model = c(0L,0L,0L), link = 0L, rho = 1.0)
cat(sprintf("eval_logf reproduces E-step logf: max|diff| = %.3g ; logg: %.3g\n",
            max(abs(ev$logf - raw$logf)), max(abs(ev$logg - raw$logg))))

score <- function(df, pars, model, link) {
  r <- emphasis:::eval_logf(pars, list(df), model = as.integer(model), link = as.integer(link), rho = 1.0)
  c(logf = r$logf, logg = r$logg)
}
fix_ts_only <- function(df) { m <- which(is_miss(df)); df$tip_start[m] <- df$brts[m]; df }
fix_pd_only <- function(df) { m <- which(is_miss(df)); ts <- df$tip_start; ts[m] <- df$brts[m]
                              df$pd <- sapply(df$brts, function(tm) hand_P(df, tm, ts)); df }
fix_both    <- function(df) fix_pd_only(fix_ts_only(df))

cases <- list(
  m_linear = list(pars = c(0.3, 0, 0.1, 0, 0.3, 0, 0.05, 0),                     model = c(0,1,0), link = 0),
  nd_exp   = list(pars = c(log(0.3), -0.01, 0, 0.2, log(0.3), 0, 0, 0.1),        model = c(1,0,1), link = 1),
  nd_gauss = list(pars = c(0.3, 0.1, 0, 0.3, 0.3, 0.1, 0, 0.2),                  model = c(1,0,1), link = 2),
  d_gauss_pdfixed = list(pars = c(0.3, 0, 0, 0.3, 0.3, 0, 0, 0.2),               model = c(0,0,1), link = 2)
)
cat("\n== which field carries the effect (max |as returned - corrected| over parent_id == -1 draws) ==\n")
for (nm in names(cases)) {
  cs <- cases[[nm]]
  for (fx in c("fix_ts_only", "fix_pd_only", "fix_both")) {
    d <- t(sapply(trees[sel_m1], function(df) score(df, cs$pars, cs$model, cs$link) -
                                             score(get(fx)(df), cs$pars, cs$model, cs$link)))
    cat(sprintf("%-16s %-12s max|dlogf| = %.4g   max|dlogg| = %.4g\n", nm, fx,
                max(abs(d[, "logf"])), max(abs(d[, "logg"]))))
  }
}

## 1. R re-implementation of Model::loglik for the D + exponential path (model.hpp:359-471)
##    with the pendant start used by the running sums made explicit.
loglik_ep_exp_R <- function(df, p, ts_runsum) {
  # ts_runsum: per-node pendant start used when a lineage is added to / removed from the running sums
  N <- nrow(df); prev <- 0; inte <- 0; loglam <- 0; logmu <- 0
  sb <- df$n[1]; sg <- df$n[1]
  for (i in seq_len(N)) {
    nd <- df[i, ]; dt <- nd$brts - prev
    M  <- if (nd$n > 0) nd$pd / nd$n else 0
    if (dt > 0) {
      A_lam <- p[1] + p[2]*nd$n + (p[3] - p[4])*M
      A_mu  <- p[5] + p[6]*nd$n + (p[7] - p[8])*M
      ei <- function(A, b, t1, t2) if (abs(b) < 1e-12) exp(A)*(t2 - t1) else exp(A)*(exp(b*t2) - exp(b*t1))/b
      inte <- inte + sb*ei(A_lam, p[4], prev, nd$brts) + sg*ei(A_mu, p[8], prev, nd$brts)
    }
    E <- if (is_ext(nd)) nd$brts - nd$tip_start else if (nd$parent_id >= 0) nd$brts - nd$focal_tip_start else M
    D <- E - M
    lam <- exp(p[1] + p[2]*nd$n + p[3]*M + p[4]*D)
    mu  <- exp(p[5] + p[6]*nd$n + p[7]*M + p[8]*D)
    if (is_ext(nd)) logmu <- logmu + log(max(mu, 1e-300)) else if (i != N) loglam <- loglam + log(lam)
    if (is_miss(nd) || is_uns(nd)) { sb <- sb + exp(-p[4]*ts_runsum[i]); sg <- sg + exp(-p[8]*ts_runsum[i]) }
    else if (is_ext(nd))           { sb <- sb - exp(-p[4]*ts_runsum[i]); sg <- sg - exp(-p[8]*ts_runsum[i]) }
    else if (i != N)               { sb <- sb + exp(-p[4]*ts_runsum[i]); sg <- sg + exp(-p[8]*ts_runsum[i]) }
    prev <- nd$brts
  }
  logmu + loglam - inte
}
p <- cases$nd_exp$pars
cat("\n== D+exp running sums: which pendant start does C++ use? (R reimplementation vs eval_logf) ==\n")
chk <- t(sapply(trees[sel_m1], function(df) {
  cpp <- score(df, p, c(1,0,1), 1)[["logf"]]
  # (A) speciation nodes add exp(-b*brts), extinction nodes subtract exp(-b*tip_start) [= t_spec for augmented]
  ts_A <- ifelse(is_ext(df), df$tip_start, df$brts)
  # (B) every node uses its stored tip_start (0 for the parent_id == -1 speciation node)
  ts_B <- df$tip_start
  c(cpp_minus_A = cpp - loglik_ep_exp_R(df, p, ts_A), cpp_minus_B = cpp - loglik_ep_exp_R(df, p, ts_B))
}))
cat(sprintf("C++ - R(A: add at brts=t_spec, subtract at t_spec): max|diff| = %.3g\n", max(abs(chk[, "cpp_minus_A"]))))
cat(sprintf("C++ - R(B: add/subtract at stored tip_start=0):      max|diff| = %.3g\n", max(abs(chk[, "cpp_minus_B"]))))
cat("=> the running sums are self-consistent for the augmented lineage (start = t_spec on both sides);\n",
    "   the inconsistency is with P(t)/M, which credits the same lineage with start 0.\n")

## 2. Size of the M distortion at the nodes following a parent_id == -1 lineage
dist <- do.call(rbind, lapply(trees[sel_m1], function(df) {
  m <- which(is_miss(df)); tsp <- df$brts[m]
  ts_true <- df$tip_start; ts_true[m] <- tsp
  P_true <- sapply(df$brts, function(tm) hand_P(df, tm, ts_true))
  k <- which(abs(df$pd - P_true) > 1e-12)
  data.frame(t_spec = tsp, M_code = df$pd[k] / df$n[k], M_true = P_true[k] / df$n[k])
}))
cat(sprintf("\nM at over-counted nodes (%d nodes): mean M_code = %.4f, mean M_true = %.4f; absolute over-count of P = t_spec (median %.3f);\n  relative over-count of M: median = %.1f%%, 90th pct = %.1f%% (M_true = 0 at the lineage's own birth node, where pd_code = t_spec)\n",
            nrow(dist), mean(dist$M_code), mean(dist$M_true), median(dist$t_spec),
            100*median(dist$M_code/dist$M_true - 1), 100*quantile(dist$M_code/dist$M_true - 1, 0.9)))

## 3. Discontinuity: the same lineage born just before vs just after the first observed split
##    (parent_id == -1 vs parent_id >= 0) gets a different pendant start.
base <- trees[[sel_m1[1]]]
m <- which(is_miss(base)); e <- which(is_ext(base))
mk <- function(df, tsp, pid) {
  df <- df[-c(m, e), ]
  # rebuild: insert speciation at tsp, extinction at t_ext (kept), recompute n via the same rule as C++
  t_ext <- base$t_ext[m]
  sp <- base[m, ]; sp$brts <- tsp; sp$parent_id <- pid; sp$tip_start <- if (pid == -1) 0 else tsp
  ex <- base[e, ]; ex$tip_start <- tsp; ex$parent_id <- pid
  df <- rbind(df, sp, ex); df <- df[order(df$brts), ]
  alive <- function(t) 2 + sum(df$brts < t & !is_ext(df) & df$t_ext > t)   # n on the segment ending at t
  df$n <- sapply(seq_len(nrow(df)), function(i) alive(df$brts[i]) + 0)  # lineages alive just before node i
  df$n <- 2 + sapply(df$brts, function(t) sum(!is_ext(df) & df$brts < t & df$t_ext > t))
  df$pd <- sapply(df$brts, function(tm) hand_P(df, tm, df$tip_start))
  df
}
base$t_ext[m] <- 3.0; base$brts[e] <- 3.0            # make the lineage live until forward time 3.0
lo <- mk(base, 0.79, -1L); hi <- mk(base, 0.81, 0L)
k_lo <- which(lo$brts == 1.5); k_hi <- which(hi$brts == 1.5)   # observed node at forward time 1.5 (age 2.5)
cat(sprintf("\nsame lineage born at 0.79 (parent_id=-1) vs 0.81 (parent_id=0), alive at the observed node t=1.5:\n  pd = %.4f vs %.4f, M = %.4f vs %.4f (n = %d); true pendant sum differs by only 0.02\n",
            lo$pd[k_lo], hi$pd[k_hi], lo$pd[k_lo]/lo$n[k_lo], hi$pd[k_hi]/hi$n[k_hi], as.integer(lo$n[k_lo])))

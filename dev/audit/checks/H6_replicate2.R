## H6 replication part 2: (1) PIT under mu_draw in a light-augmentation regime
## (few missing, no zero-weight rejections) to show the mild non-uniformity in
## part 1 is acceptance-filter selection, not the lifetime draw; (2) decompose
## the compensator gap for beta_M-only into H48 (pd_final vs pd_sampler at the
## segment end) and H6 proper (M frozen at the segment end vs M(t) integrated).
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages(library(emphasis))
options(width = 120)
src <- readLines("/Users/pancho/Code/emphasis/dev/audit/checks/H6.R")
cut <- grep("^## ---- driver", src); eval(parse(text = src[1:(cut - 1)]))
set.seed(7)
crown <- 6; brts <- c(crown, sort(runif(10, 0.2, crown - 0.3), decreasing = TRUE)); TT <- crown  # 12 tips
draw <- function(p8, model_bin, link, N)
  emphasis:::augment_trees(brts = brts, pars = p8, sample_size = as.integer(N), maxN = 200L*N,
                           max_missing = 10000L, max_lambda = 2000, num_threads = 1L,
                           model = as.integer(model_bin), link = as.integer(link), rho = 1)
ptrunc <- function(l, mu, r) (1 - exp(-mu*l)) / (1 - exp(-mu*r))
ks <- function(u) suppressWarnings(ks.test(u, "punif")$p.value)

## (1) light regime: lambda small, gamma_M large -> mu_credit far from mu_draw, few missing
p8 <- c(log(0.25), 0, 0, 0, log(0.15), 0, 1.0, 0); mb <- c(0L,1L,0L); link <- 1L
E <- draw(p8, mb, link, 300)
ev <- do.call(rbind, lapply(seq_along(E$trees), function(i) { qt <- q_true_R(E$trees[[i]], p8, mb, link); qt$ev }))
cat(sprintf("(1) light regime: trees 300, augmented nodes %d, mean #missing %.2f, lw range [%.1f, %.1f]\n",
            nrow(ev), nrow(ev)/300, min(E$logf - E$logg), max(E$logf - E$logg)))
cat(sprintf("    PIT KS p: mu_draw %.3g | mu_credit %.3g;  mean life obs %.4f (se %.4f), E[mu_draw] %.4f, E[mu_credit] %.4f; mu_credit range [%.3f, %.3f]\n",
            ks(ptrunc(ev$life, ev$mu_draw, TT-ev$t_spec)), ks(ptrunc(ev$life, ev$mu_credit, TT-ev$t_spec)),
            mean(ev$life), sd(ev$life)/sqrt(nrow(ev)),
            mean(1/ev$mu_draw - (TT-ev$t_spec)*exp(-ev$mu_draw*(TT-ev$t_spec))/(1-exp(-ev$mu_draw*(TT-ev$t_spec)))),
            mean(1/ev$mu_credit - (TT-ev$t_spec)*exp(-ev$mu_credit*(TT-ev$t_spec))/(1-exp(-ev$mu_credit*(TT-ev$t_spec)))),
            min(ev$mu_credit), max(ev$mu_credit)))

## (2) decomposition, beta_M only (mu constant => survival integral exact; only lambda(t) matters)
p8 <- c(log(0.5), 0, 0.25, 0, log(0.2), 0, 0, 0); mb <- c(0L,1L,0L); link <- 1L
E <- draw(p8, mb, link, 150)
res <- t(sapply(seq_along(E$trees), function(i) {
  df <- E$trees[[i]]
  qc <- q_cpp_R(df, p8, mb, link); qt <- q_true_R(df, p8, mb, link)
  ## compensator with the C++ formula (M frozen at segment-end node) but pd = pd_sampler(brts)
  ext <- is_ext(df$t_ext); mis <- is_mis(df$t_ext); sp <- !ext
  ts_s <- ifelse(ext, NA, ifelse(mis, df$brts, 0))
  pd_s <- function(t) { k <- sp & df$brts <= t & df$t_ext > t; sum(t - ts_s[k]) }
  mu <- exp(p8[5]); prev <- 0; inte_s <- 0
  for (j in seq_len(nrow(df))) {
    nd <- df[j, ]; dt <- nd$brts - prev
    lam <- exp(p8[1] + p8[3] * pd_s(nd$brts)/nd$n)
    inte_s <- inte_s + nd$n*lam*(dt - (1/mu)*(exp(-mu*(TT-nd$brts)) - exp(-mu*(TT-prev))))
    prev <- nd$brts
  }
  c(cpp_minus_true = qc$inte - qt$inte, h48_part = qc$inte - inte_s, freeze_part = inte_s - qt$inte)
}))
cat(sprintf("(2) beta_M-only compensator gap (150 trees): total %+.4f | H48 part (pd_final vs pd_sampler) %+.4f | H6 proper (M frozen at segment end vs int M(t)) %+.4f\n",
            mean(res[,1]), mean(res[,2]), mean(res[,3])))
cat(sprintf("    sd: total %.4f, H48 %.4f, freeze %.4f;  freeze part always >= 0? %s\n",
            sd(res[,1]), sd(res[,2]), sd(res[,3]), all(res[,3] >= -1e-12)))

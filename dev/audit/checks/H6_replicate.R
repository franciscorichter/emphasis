## H6 replication — independent variation of dev/audit/checks/H6.R.
## Reuses the verifier's q_cpp_R / q_true_R (sourced from H6.R up to its driver)
## but changes: the tree (random 16-tip, crown age 8), the parameter sets, adds
## an N-only exponential model with gamma_N != 0 (boundary of the claim: pd/D
## absent -> gap must be 0 even though mu varies with n), and adds a
## reconstruction-free test of mechanism (a): PIT of each lifetime under
## mu_draw (pd = 0) vs under mu_credit, KS against U(0,1).
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages(library(emphasis))
options(width = 120)
src <- readLines("/Users/pancho/Code/emphasis/dev/audit/checks/H6.R")
cut <- grep("^## ---- driver", src)
eval(parse(text = src[1:(cut - 1)]))

set.seed(20260913)
crown <- 8
brts <- c(crown, sort(runif(14, 0.2, crown - 0.3), decreasing = TRUE))   # 16 tips
cat("brts:", round(brts, 3), "\n")
TT <- crown

draw <- function(p8, model_bin, link, N) {
  emphasis:::augment_trees(brts = brts, pars = p8, sample_size = as.integer(N),
                           maxN = 200L*N, max_missing = 10000L, max_lambda = 2000,
                           num_threads = 1L, model = as.integer(model_bin),
                           link = as.integer(link), rho = 1)
}
etrunc <- function(mu, r) 1/mu - r*exp(-mu*r)/(1 - exp(-mu*r))
ptrunc <- function(l, mu, r) (1 - exp(-mu*l)) / (1 - exp(-mu*r))
ess <- function(lw) { w <- exp(lw - max(lw)); sum(w)^2/sum(w^2) }

analyse <- function(label, p8, model_bin, link, N = 120) {
  E <- draw(p8, model_bin, link, N)
  d_cpp <- d6 <- d_inte <- numeric(N); ev_all <- data.frame()
  for (i in seq_len(N)) {
    d <- E$trees[[i]]
    qc <- q_cpp_R(d, p8, model_bin, link); qt <- q_true_R(d, p8, model_bin, link)
    d_cpp[i] <- qc$logg - E$logg[i]
    parent_true <- if (nrow(qt$ev)) sum(ifelse(qt$ev$K > 0, -log(qt$ev$K), 0)) else 0
    d6[i] <- E$logg[i] - (qt$logq - parent_true + sum(qc$parent_terms))
    d_inte[i] <- qc$inte - qt$inte
    if (nrow(qt$ev)) ev_all <- rbind(ev_all, cbind(tree = i, qt$ev))
  }
  lw_cpp <- E$logf - E$logg; lw_fix <- E$logf - (E$logg - d6)
  cat(sprintf("\n=== %s  pars8=(%s) model_bin=(%s) link=%d N=%d mean#missing=%.2f\n", label,
              paste(signif(p8,3), collapse=","), paste(model_bin, collapse=","), link, N,
              mean(sapply(E$trees, function(d) sum(is_mis(d$t_ext))))))
  cat(sprintf("  transcription max|q_cpp_R - logg_cpp| = %.2e\n", max(abs(d_cpp))))
  cat(sprintf("  H6-only gap logg_cpp - logq_true: mean %+.4f sd %.4f max|.| %.3g\n", mean(d6), sd(d6), max(abs(d6))))
  cat(sprintf("  compensator inte_cpp - int nh:    mean %+.4f sd %.4f max|.| %.3g\n", mean(d_inte), sd(d_inte), max(abs(d_inte))))
  if (nrow(ev_all)) {
    ev <- ev_all
    cat(sprintf("  nodes=%d  max|lam_nh-lam_credit|=%.3g  max|mu_draw-mu_nh|=%.3g  max|mu_draw-mu_credit|=%.3g\n",
                nrow(ev), max(abs(ev$lam_nh-ev$lam_credit)), max(abs(ev$mu_draw-ev$mu_nh)), max(abs(ev$mu_draw-ev$mu_credit))))
    ## reconstruction-free mechanism test: PIT of lifetimes
    u_draw <- ptrunc(ev$life, ev$mu_draw, TT - ev$t_spec)
    u_cred <- ptrunc(ev$life, ev$mu_credit, TT - ev$t_spec)
    u_nh   <- ptrunc(ev$life, ev$mu_nh, TT - ev$t_spec)
    ks <- function(u) suppressWarnings(ks.test(u, "punif")$p.value)
    cat(sprintf("  lifetime PIT KS p-values: under mu_draw(pd=0) %.3g | under mu_nh %.3g | under mu_credit %.3g\n",
                ks(u_draw), ks(u_nh), ks(u_cred)))
    cat(sprintf("  mean lifetime obs %.4f (se %.4f) | E[mu_draw] %.4f | E[mu_nh] %.4f | E[mu_credit] %.4f\n",
                mean(ev$life), sd(ev$life)/sqrt(nrow(ev)), mean(etrunc(ev$mu_draw, TT-ev$t_spec)),
                mean(etrunc(ev$mu_nh, TT-ev$t_spec)), mean(etrunc(ev$mu_credit, TT-ev$t_spec))))
  }
  gap <- lmexp(lw_cpp) - lmexp(lw_fix)
  bs <- replicate(200, { j <- sample(N, replace = TRUE); lmexp(lw_cpp[j]) - lmexp(lw_fix[j]) })
  cat(sprintf("  fhat C++ %.3f (ESS %.0f) | H6-corrected %.3f (ESS %.0f) | gap %+.3f (bs se %.3f)\n",
              lmexp(lw_cpp), ess(lw_cpp), lmexp(lw_fix), ess(lw_fix), gap, sd(bs)))
  invisible(NULL)
}

## boundary: N-only, exponential link, mu depends on n (no pd, no D) -> claim says gap == 0
analyse("N-only exp, gamma_N != 0", c(log(0.7), -0.02, 0, 0, log(0.2), 0.03, 0, 0), c(1L,0L,0L), 1L)
## M-only linear link (verifier used only exp for M) -- moderate coefficients
analyse("M-only linear",          c(0.6, 0, 0.08, 0, 0.2, 0, 0.12, 0),            c(0L,1L,0L), 0L)
## M-only exp, mu-side only (mechanism a vs c)
analyse("M-only exp gamma_M only", c(log(0.5), 0, 0, 0, log(0.2), 0, 0.6, 0),      c(0L,1L,0L), 1L)
## M-only exp, lambda-side only (mechanism: compensator + H48 only)
analyse("M-only exp beta_M only",  c(log(0.5), 0, 0.2, 0, log(0.2), 0, 0, 0),      c(0L,1L,0L), 1L)
## D-only exp (ep_exp compensator branch)
analyse("D-only exp",              c(log(0.5), 0, 0, 0.15, log(0.2), 0, 0, 0.3),   c(0L,0L,1L), 1L)
## full nd-like: N + D linear
analyse("N+D linear",              c(0.7, -0.01, 0, 0.08, 0.2, 0, 0, 0.1),         c(1L,0L,1L), 0L)

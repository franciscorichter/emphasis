# (b) Seam: weight conventions across BDI (mean-1), thinning E-step
# (max-scaled) and the M-step objective (sum_w scaling of the conditional).
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib-wave1", .libPaths()))
library(emphasis)

set.seed(11)
brts <- sort(ape::branching.times(ape::rcoal(20)), decreasing = TRUE)
# dd model: BDI weights are non-degenerate (CR weights are exactly constant)
pars8 <- function(l, bN, m) c(l, bN, 0, 0, m, 0, 0, 0)
mb <- c(1L, 0L, 0L)
lb8 <- c(0, -1, 0, 0, 0, -1, 0, 0); ub8 <- c(3, 0, 0, 0, 3, 0, 0, 0)
init <- pars8(0.8, -0.01, 0.2)

# one BDI E-step -> trees + log-weights
e <- emphasis:::.augment_tree_bdi(tree = brts, pars = init,
                                  model_bin = mb,
                                  sample_size = 60L, max_missing = 1e4L,
                                  link = 0L, rho = 1)
lw <- e$weights
ok <- is.finite(lw) & is.finite(e$logf)
trees <- e$trees[ok]; lw <- lw[ok]
cat("n trees:", length(trees), " sd(lw):", stats::sd(lw), "\n")

mk <- function(w) list(trees = trees, weights = w, rejected = 0L,
                       rejected_overruns = 0L, rejected_lambda = 0L,
                       rejected_zero_weights = 0L, time = 0, fhat = 0)

w_mean1 <- { v <- exp(lw - max(lw)); v / sum(v) * length(v) }   # BDI convention
w_max   <- exp(lw - max(lw))                                    # thinning convention
cat("sum_w mean1 =", sum(w_mean1), "  sum_w max-scaled =", sum(w_max), "\n")

# exact CR crown-survival conditional: log P(both crown clades survive to T)
Tc <- max(brts)
logP <- function(p8) {
  lam <- p8[1]; mu <- p8[5]
  if (lam <= 0) return(-700)
  d <- lam - mu
  ps <- if (abs(d) < 1e-12) (lam * Tc) / (1 + lam * Tc) else
    (d) / (lam - mu * exp(-d * Tc)) * (1 - exp(-d * Tc)) * lam / d
  ps <- 1 - (mu * (1 - exp(-d * Tc))) / (lam - mu * exp(-d * Tc))
  log(max(min(ps, 1) ^ 2, 1e-300))
}

run <- function(w, cond) {
  r <- emphasis:::m_cpp(e_step = mk(w), init_pars = init, plugin = "rpd1",
             lower_bound = lb8, upper_bound = ub8, xtol_rel = 1e-6,
             num_threads = 1L, model = mb, link = 0L, rho = 1,
             rconditional = cond)
  as.numeric(r$estimates)[c(1, 2, 5)]
}

cat("\n## unconditioned: scale invariance of the weights\n")
u1 <- run(w_mean1, NULL); u2 <- run(w_max, NULL); u3 <- run(7 * w_mean1, NULL)
cat("  mean-1      :", u1, "\n  max-scaled  :", u2, "\n  7x mean-1   :", u3, "\n")
cat("  max|mean1 - maxscaled| =", max(abs(u1 - u2)), "\n")
cat("  max|mean1 - 7x|        =", max(abs(u1 - u3)), "\n")

cat("\n## conditioned: same three weightings\n")
c1 <- run(w_mean1, logP); c2 <- run(w_max, logP); c3 <- run(7 * w_mean1, logP)
cat("  mean-1      :", c1, "\n  max-scaled  :", c2, "\n  7x mean-1   :", c3, "\n")
cat("  max|mean1 - maxscaled| =", max(abs(c1 - c2)), "\n")
cat("  max|mean1 - 7x|        =", max(abs(c1 - c3)), "\n")
cat("  conditioned moves the estimate:", max(abs(c1 - u1)), "\n")

cat("\n## same on the PRE-FIX build (for contrast)\n")
saveRDS(list(trees = trees, lw = lw, init = init, lb8 = lb8, ub8 = ub8, Tc = Tc),
        "/tmp/emphasis_review_estep.rds")

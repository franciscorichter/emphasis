# H99 — "No test in R or C++ exercises link_int = 2 (gaussian), the link every
# README example recommends."
#
# Part A: coverage fact — count references to the gaussian link under tests/.
# Part B: does the gaussian link, as compiled, implement README:92
#         f(beta_0, eta) = beta_0 * exp(-(eta - 1)^2 / 2), with beta_0 the peak
#         (reached at eta_cov = 1)?  Deterministic checks through eval_logf on
#         hand-built trees, an R reference likelihood on augmented trees, and the
#         CR mapping lambda = beta_0 * e^{-1/2}.
# Part C: forward simulator gaussian branch vs the linear branch at the mapped
#         parameters (same law -> same tip-count distribution).
# Part D: smoke — estimate_rates(link = "gaussian") runs end to end (thinning MCEM).

.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({ library(emphasis); library(ape) })
pkg <- "/Users/pancho/Code/emphasis"

cat("emphasis", as.character(packageVersion("emphasis")), "\n\n")

# ---------------------------------------------------------------- Part A ----
cat("== A. test coverage of link 2 ==\n")
test_files <- list.files(file.path(pkg, "tests"), recursive = TRUE, full.names = TRUE)
hits <- unlist(lapply(test_files, function(f) {
  l <- readLines(f, warn = FALSE)
  i <- grep("gaussian|link *= *2|link_int *= *2|LinkType::gaussian", l)
  if (length(i)) paste0(basename(f), ":", i, ": ", l[i]) else character(0)
}))
cat("files scanned:", length(test_files), " hits for gaussian/link=2:", length(hits), "\n")
print(hits)
src_tests <- list.files(file.path(pkg, "src"), pattern = "test", recursive = TRUE)
cat("C++ test files under src/:", length(src_tests), "\n")
readme <- readLines(file.path(pkg, "README.md"))
cat("README lines using link = \"gaussian\":", paste(grep('link *= *"gaussian"', readme), collapse = ","), "\n")
skips <- sum(grepl("^\\s*skip\\(", readLines(file.path(pkg, "tests/testthat/test-de.R"))))
cat("skip() calls in test-de.R (C++ integration tests):", skips, "\n\n")

# ---------------------------------------------------------------- Part B ----
cat("== B. gaussian link implements README:92 and beta_0 is the peak ==\n")
TIP <- 10e10
ok <- function(lbl, a, b, tol = 1e-9) {
  d <- max(abs(a - b))
  cat(sprintf("  %-58s |diff| = %.2e  %s\n", lbl, d, if (d < tol) "PASS" else "FAIL"))
  invisible(d < tol)
}
g <- function(b0, eta) b0 * exp(-(eta - 1)^2 / 2)      # README:92, written independently

# B1: 2-tip tree, n == 2 throughout, model dd (use_N). With beta_N = gamma_N = 1/2,
#     eta_cov = 1 for every node -> rates are exactly beta_0, gamma_0 (the "peak").
T2 <- 3.7; b0 <- 0.8; g0 <- 0.3
df2 <- data.frame(brts = T2, n = 2, t_ext = TIP, pd = 2 * T2)
lf <- emphasis:::eval_logf(c(b0, 0.5, 0, 0, g0, 0.5, 0, 0), list(df2), c(1L, 0L, 0L), 2L)$logf
ok("B1 eta_cov=1: logf == -2T(beta_0+gamma_0)", lf, -2 * T2 * (b0 + g0))
lf0 <- emphasis:::eval_logf(c(b0, 0, 0, 0, g0, 0, 0, 0), list(df2), c(1L, 0L, 0L), 2L)$logf
ok("B1' eta_cov=0: logf == -2T e^{-1/2}(beta_0+gamma_0)", lf0, -2 * T2 * exp(-0.5) * (b0 + g0))
# Peak property: scanning beta_N, logf (= -2T(lambda+mu)) is minimised where lambda is max, i.e. beta_N = 1/2
bN <- seq(0, 1, by = 0.05)
lf_scan <- sapply(bN, function(b) emphasis:::eval_logf(c(b0, b, 0, 0, 0, 0, 0, 0), list(df2), c(1L, 0L, 0L), 2L)$logf)
cat(sprintf("  B1'' argmax of implied lambda(2) over beta_N grid: %.2f (expected 0.50); lambda there = %.4f (beta_0 = %.2f)\n",
            bN[which.min(lf_scan)], -min(lf_scan) / (2 * T2), b0))

# B2: 3-tip tree brts(backward) = c(T, t1): one speciation at n = 2, then n = 3 to the present.
T3 <- 4; t1 <- 1.5
df3 <- data.frame(brts = c(T3 - t1, T3), n = c(2, 3), t_ext = c(TIP, TIP), pd = c(2 * (T3 - t1), 2 * (T3 - t1) + 3 * t1))
pars3 <- c(0.9, 0.3, 0, 0, 0.2, -0.1, 0, 0)
lam <- function(n) g(pars3[1], pars3[2] * n); mu <- function(n) g(pars3[5], pars3[6] * n)
ref3 <- log(lam(2)) - (T3 - t1) * 2 * (lam(2) + mu(2)) - t1 * 3 * (lam(3) + mu(3))
lf3 <- emphasis:::eval_logf(pars3, list(df3), c(1L, 0L, 0L), 2L)$logf
ok("B2 3-tip tree, gaussian dd vs README-formula reference", lf3, ref3)

# B3: CR mapping on augmented trees: gaussian(beta_0, gamma_0) == linear(beta_0 e^-1/2, gamma_0 e^-1/2)
brts <- c(4, 2.5, 1.2, 0.6)
aug <- emphasis:::augment_trees(brts, c(0.5, 0, 0, 0, 0.1, 0, 0, 0), 30L, 300L, 1000L, 100, 1L, c(0L, 0L, 0L), 0L)
pg <- c(0.5, 0, 0, 0, 0.1, 0, 0, 0)
pl <- pg * exp(-0.5)
lg <- emphasis:::eval_logf(pg, aug$trees, c(0L, 0L, 0L), 2L)$logf
ll <- emphasis:::eval_logf(pl, aug$trees, c(0L, 0L, 0L), 0L)$logf
ok(sprintf("B3 CR: gaussian vs linear at mapped pars (%d aug. trees)", length(lg)), lg, ll)

# B4: R reference likelihood (piecewise-constant path of model.hpp:loglik) for
#     model c(1,1,0) (N and M covariates) under the README gaussian formula,
#     on augmented trees with extinct lineages.
ref_logf <- function(df, p) {
  M <- ifelse(df$n > 0, df$pd / df$n, 0)
  lam <- g(p[1], p[2] * df$n + p[3] * M)
  mu  <- g(p[5], p[6] * df$n + p[7] * M)
  dt  <- diff(c(0, df$brts))
  inte <- sum(dt * df$n * (lam + mu))
  ext <- df$t_ext == 0
  n_r <- nrow(df)
  ev  <- sum(log(mu[ext])) + sum(log(lam[!ext & seq_len(n_r) != n_r]))
  ev - inte
}
p4 <- c(0.7, 0.15, 0.05, 0, 0.25, 0.1, -0.2, 0)
aug4 <- emphasis:::augment_trees(brts, c(0.6, 0, 0, 0, 0.2, 0, 0, 0), 30L, 300L, 1000L, 100, 1L, c(0L, 0L, 0L), 0L)
n_ext <- sum(sapply(aug4$trees, function(d) sum(d$t_ext == 0)))
lf4 <- emphasis:::eval_logf(p4, aug4$trees, c(1L, 1L, 0L), 2L)$logf
rf4 <- sapply(aug4$trees, ref_logf, p = p4)
ok(sprintf("B4 dd+M gaussian on %d aug. trees (%d extinct lineages) vs R reference", length(lf4), n_ext), lf4, rf4)
cat("\n")

# ---------------------------------------------------------------- Part C ----
cat("== C. forward simulator: gaussian vs linear at mapped parameters ==\n")
sim_tips <- function(pars, link, R = 300, max_t = 3) {
  out <- replicate(R, {
    s <- emphasis:::simulate_div_tree_cpp(pars, c(0L, 0L, 0L), max_t, 2000L, 1L, link)
    if (s$status == "done") sum(s$Ltable[, 4] == -1) else NA_real_
  })
  out
}
tg <- sim_tips(c(1.0, 0, 0, 0, 0.2, 0, 0, 0), 2L)
tl <- sim_tips(c(1.0, 0, 0, 0, 0.2, 0, 0, 0) * exp(-0.5), 0L)
cat(sprintf("  done fraction: gaussian %.3f  linear %.3f\n", mean(!is.na(tg)), mean(!is.na(tl))))
cat(sprintf("  mean tips (done): gaussian %.2f (se %.2f)  linear %.2f (se %.2f)\n",
            mean(tg, na.rm = TRUE), sd(tg, na.rm = TRUE) / sqrt(sum(!is.na(tg))),
            mean(tl, na.rm = TRUE), sd(tl, na.rm = TRUE) / sqrt(sum(!is.na(tl)))))
wt <- wilcox.test(tg, tl); cat(sprintf("  Wilcoxon p = %.3f (same law expected -> not small)\n", wt$p.value))
# unconditional theory for comparison: E[N(T)] = 2 exp((lam-mu) T) with lam = e^-1/2, mu = 0.2 e^-1/2
cat(sprintf("  unconditional E[N(T)] = %.2f (conditioned-on-done means above must exceed it)\n", 2 * exp((1 - 0.2) * exp(-0.5) * 3)))
cat("\n")

# ---------------------------------------------------------------- Part D ----
cat("== D. smoke: estimate_rates(link = 'gaussian') end to end ==\n")
set.seed(7)
tr <- NULL
while (is.null(tr)) {
  x <- ape::rlineage(0.5, 0.1, Tmax = 6)
  x <- tryCatch(ape::drop.fossil(x), error = function(e) NULL)
  if (!is.null(x) && Ntip(x) >= 12 && Ntip(x) <= 40) tr <- x
}
cat("  tree tips:", Ntip(tr), "\n")
t0 <- Sys.time()
fit <- tryCatch(
  estimate_rates(tr, method = "mcem", model = "cr", link = "gaussian",
                 control = list(lower_bound = c(0.05, 0.001), upper_bound = c(3, 1.5),
                                max_iter = 4L, max_time = 90, sample_size = 30L,
                                maxN = 300L, num_threads = 1L)),
  error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL })
cat(sprintf("  elapsed %.1fs\n", as.numeric(Sys.time() - t0, units = "secs")))
if (!is.null(fit)) {
  print(str(fit, max.level = 1))
  p <- tryCatch(fit$pars, error = function(e) NULL); if (is.null(p)) p <- tryCatch(fit$estimates, error = function(e) NULL)
  cat("  pars:", paste(format(unlist(p), digits = 4), collapse = " "), "\n")
  cat("  details names:", paste(names(fit$details), collapse = ", "), "\n")
  it <- fit$details$iterations; if (is.null(it)) it <- fit$details$n_iter
  cat("  iterations recorded:", if (is.null(it)) "?" else it, "\n")
  cat("  sampler:", if (is.null(fit$details$sampler)) "?" else fit$details$sampler, "\n")
  cat(sprintf("  truth on gaussian scale: beta_0 = 0.5 e^{1/2} = %.3f, gamma_0 = 0.1 e^{1/2} = %.3f (4 iterations only)\n",
              0.5 * exp(0.5), 0.1 * exp(0.5)))
}

# What a gaussian test in the R layer would see (cf. H73): auto_bounds centre sign.
cat("\n== E. auto_bounds under the gaussian link (R layer, cf. H73) ==\n")
ab <- tryCatch(auto_bounds(tr, model = "cr", link = "gaussian", verbose = FALSE),
               error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL })
if (!is.null(ab)) {
  cat("  lower_bound:", paste(format(ab$lower_bound, digits = 4), collapse = " "), "\n")
  cat("  upper_bound:", paste(format(ab$upper_bound, digits = 4), collapse = " "), "\n")
  cat("  beta_0 lower bound > 0 ?", ab$lower_bound[1] > 0, "\n")
}
cat("\ndone\n")

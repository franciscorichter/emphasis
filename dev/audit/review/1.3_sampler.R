# Item 1.3 review: sampler-level attack (POST-FIX build)
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib-wave1", .libPaths()))
suppressWarnings(library(emphasis))
aug <- emphasis:::.augment_tree_bdi
print(args(aug))
brts11 <- c(5, 4.745201, 4.530461, 4.067871, 3.622029, 2.929002,
            1.468698, 1.386875, 1.302139, 0.044729)
bd <- function(brts, lam, mu) DDD::bd_loglik(pars1 = c(lam, mu, 0, 0), pars2 = c(0, 0, 1, 0, 2),
                                            brts = brts, missnumspec = 0)
chk <- function(brts, pars, link = 0L, ss = 30L, seed = 1, mm = 1e4L, label = "") {
  set.seed(seed)
  t0 <- proc.time()[3]
  a <- withCallingHandlers(
    tryCatch(aug(brts, pars = if (link == 1L) log(pars) else pars, model_bin = c(0L, 0L, 0L),
                 sample_size = ss, link = link, rho = 1, max_missing = mm),
             error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }),
    warning = function(w) { cat("  WARN:", conditionMessage(w), "\n"); invokeRestart("muffleWarning") })
  el <- proc.time()[3] - t0
  if (is.null(a)) return(invisible(NULL))
  ref <- bd(brts, pars[1], pars[2])
  nmiss <- vapply(a$trees, function(tr) sum(!is.na(tr$t_ext) & is.finite(tr$t_ext)) , numeric(1))
  cat(sprintf("%-32s lam=%g mu=%g link=%d ss=%d: n_valid=%d n_att=%d rej=%d rej_mm=%d  sd(lw)=%.1e  fhat-bd=%.2e  missing/draw: min %d med %d max %d  %.1fs\n",
              label, pars[1], pars[2], link, ss, a$n_valid, a$n_attempts, a$n_rejected, a$n_rejected_max_missing,
              if (length(a$weights) > 1) sd(a$weights) else NA, a$fhat - ref,
              min(nmiss), median(nmiss), max(nmiss), el))
  invisible(a)
}

cat("== A. one tree structure (to know the columns) ==\n")
set.seed(1); a <- aug(brts11, pars = c(0.3, 0.5), model_bin = c(0L,0L,0L), sample_size = 2L, link = 0L, rho = 1)
str(a$trees[[1]]); print(names(a))

cat("\n== B. mu >> lam, mu = lam, mu slightly > lam, on brts11 ==\n")
chk(brts11, c(0.5, 0.3), label = "baseline mu<lam")
chk(brts11, c(0.3, 0.5), label = "mu>lam")
chk(brts11, c(0.4, 0.4), label = "critical")
chk(brts11, c(0.4, 0.4 + 1e-13), label = "critical inside switch")
chk(brts11, c(0.4, 0.4 + 1e-10), label = "just outside switch")
chk(brts11, c(0.1, 1.0), label = "mu = 10 lam")
chk(brts11, c(0.1, 3.0), label = "mu = 30 lam")
chk(brts11, c(0.05, 8), label = "mu=8 tp=5 (dtau=40)")
chk(brts11, c(1e-6, 1e-6), label = "tiny critical")
chk(brts11, c(1e-6, 2e-6), label = "tiny mu>lam")
chk(brts11, c(0, 0.5), label = "lam=0")
chk(brts11, c(0.5, 0), label = "mu=0")
chk(brts11, c(0, 0), label = "both 0")
chk(brts11, c(0.3, 0.5), link = 1L, label = "exp link mu>lam")
chk(brts11, c(0.4, 0.4), link = 1L, label = "exp link critical")
chk(brts11, c(0.3, 0.5), ss = 1L, label = "sample_size 1")

cat("\n== C. tiny trees ==\n")
chk(c(5), c(0.3, 0.5), label = "2 tips mu>lam")
chk(c(5), c(0.4, 0.4), label = "2 tips critical")
chk(c(5, 4.99999), c(0.3, 0.5), label = "3 tips near-simultaneous")
chk(c(5, 1e-9), c(0.3, 0.5), label = "3 tips split at ~tp")

cat("\n== D. larger / longer trees ==\n")
set.seed(7)
b60 <- sort(c(30, runif(58, 0, 30)), decreasing = TRUE)
chk(b60, c(0.3, 0.2), label = "60 tips tp=30 mu<lam")
chk(b60, c(0.2, 0.3), label = "60 tips tp=30 mu>lam")
chk(b60, c(0.25, 0.25), label = "60 tips tp=30 critical")
chk(b60, c(0.1, 1.0), label = "60 tips tp=30 mu=10lam")
b100 <- sort(c(100, runif(98, 0, 100)), decreasing = TRUE)
chk(b100, c(0.05, 0.08), label = "100 tips tp=100 mu>lam")
chk(b100, c(0.05, 0.05), label = "100 tips tp=100 critical")
chk(b100, c(0.05, 8), label = "100 tips tp=100 mu=8 (overflow)")
chk(b100, c(0.05, 5), label = "100 tips tp=100 mu=5 (dtau=495)")
chk(b100, c(0.05, 7.3), label = "100 tips tp=100 mu=7.3 (dtau=725)")

cat("\n== E. zero-length / duplicated event times (spurious events) at mu>lam ==\n")
set.seed(3)
a <- aug(brts11, pars = c(0.3, 0.5), model_bin = c(0L,0L,0L), sample_size = 200L, link = 0L, rho = 1)
dup <- 0L; zero_len <- 0L; at_brts <- 0L
for (tr in a$trees) {
  tb <- tr$brts; te <- tr$t_ext
  if (any(duplicated(tb[!is.na(tb)]))) dup <- dup + 1L
  if (any(te - tb <= 1e-12, na.rm = TRUE)) zero_len <- zero_len + 1L
  miss_b <- tb[!is.na(te) & is.finite(te)]
  if (any(abs(outer(miss_b, brts11, "-")) < 1e-12)) at_brts <- at_brts + 1L
}
cat(sprintf("200 draws: trees with duplicated brts=%d, zero-length missing lineages=%d, missing births at observed brts=%d; sd(lw)=%.1e\n",
            dup, zero_len, at_brts, sd(a$weights)))

cat("\n== F. thinning sampler cross-check at mu>lam (estimate_likelihood / other sampler) ==\n")
print(grep("augment|likelihood|sim_tree", getNamespaceExports("emphasis"), value = TRUE))
print(grep("augment|likelihood", ls(asNamespace("emphasis"), all.names = TRUE), value = TRUE))

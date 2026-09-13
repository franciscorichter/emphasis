## H1 replication: vary tree/generator, model, link, sampler, method; test bias-vs-noise at the cut.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({library(emphasis); library(TreeSim); library(DDD); library(ape)})
lam <- 0.2; mu <- 0.05
p8 <- function(l, m, bN = 0, gN = 0) c(l, bN, 0, 0, m, gN, 0, 0)
aug <- function(brts, pars, N, maxN, model = c(0L,0L,0L), link = 0L) tryCatch(
  emphasis:::augment_trees(brts, pars, as.integer(N), as.integer(maxN), 10000L, 1e6, 1L, as.integer(model), as.integer(link), 1.0),
  error = function(e) e)
fh <- function(r) { lw <- r$logf - r$logg; m <- max(lw); log(sum(exp(lw - m))) + m - log(length(lw) + r$rejected_zero_weights) }
rep_line <- function(r, lab) if (inherits(r, "error")) cat(sprintf("%-34s ERROR: %s\n", lab, sub("with rejection reasons: ", "", conditionMessage(r)))) else {
  lw <- r$logf - r$logg
  cat(sprintf("%-34s kept=%d zero_w=%d overrun=%d lam=%d rej=%d  lw min=%.1f med=%.1f max=%.1f  fhat=%.2f\n", lab, length(lw), r$rejected_zero_weights,
              r$rejected_overruns, r$rejected_lambda, r$rejected, min(lw), median(lw), max(lw), fh(r))) }

cat("=== 1. different generator: ape::rlineage 40-tip tree, linear link, CR ===\n")
set.seed(7); repeat { t <- ape::rlineage(lam, mu, Tmax = 25); t2 <- ape::drop.fossil(t); if (Ntip(t2) >= 35 && Ntip(t2) <= 60) break }
b <- sort(as.numeric(branching.times(t2)), decreasing = TRUE); n <- length(b) + 1
cat("n_tips =", n, " crown =", round(b[1], 2), "\n")
for (s in c(1, 1e3, 1e5, 1e7, 1e8)) rep_line(aug(b * s, p8(lam / s, mu / s), 20, 200), sprintf("rlineage s=%g pred shift %.0f", s, -(n - 2) * log(s)))

cat("\n=== 2. exponential link (log-pars), same tree ===\n")
for (s in c(1, 1e5, 1e8)) rep_line(aug(b * s, p8(log(lam / s), log(mu / s)), 20, 200, link = 1L), sprintf("explink s=%g", s))

cat("\n=== 3. dd model (N-dependent), TreeSim 40-tip tree, seed 41 ===\n")
set.seed(41); tt <- TreeSim::sim.bd.taxa(40, 1, 0.3, 0.05, complete = FALSE)[[1]]
b40 <- sort(as.numeric(branching.times(tt)), decreasing = TRUE); n40 <- 40
for (s in c(1, 1e5, 1e8)) rep_line(aug(b40 * s, p8(0.3 / s, 0.05 / s, bN = -0.002 / s), 20, 200, model = c(1L,0L,0L)), sprintf("dd s=%g", s))

cat("\n=== 4. bias vs noise at the cut: 8 replicates s=1 vs 8 at straddle s* (rlineage tree) ===\n")
r1 <- aug(b, p8(lam, mu), 40, 400); med1 <- median(r1$logf - r1$logg)
f1 <- replicate(8, fh(aug(b, p8(lam, mu), 40, 400)))
for (target in c(-743, -745)) {
  s_str <- exp((target - med1) / (-(n - 2)))
  rs <- lapply(1:8, function(i) aug(b * s_str, p8(lam / s_str, mu / s_str), 40, 4000))
  ok <- !sapply(rs, inherits, "error")
  fs <- sapply(rs[ok], fh); zw <- sapply(rs[ok], `[[`, "rejected_zero_weights"); mn <- sapply(rs[ok], function(r) min(r$logf - r$logg))
  ## upper bound on the bias from dropping weights <= exp(-745.13) but counting them in S:
  bound <- sapply(rs[ok], function(r) { lw <- r$logf - r$logg; m <- max(lw); log(1 - min(1, r$rejected_zero_weights * exp(-745.13 - m) / sum(exp(lw - m)))) })
  cat(sprintf("target med lw %d: s*=%.3g  runs ok=%d/8  zero_w=%s  min lw=%.2f..%.2f\n", target, s_str, sum(ok), paste(zw, collapse=","), min(mn), max(mn)))
  cat(sprintf("   fhat(s=1): mean %.2f sd %.2f ; predicted at s*: %.2f ; observed fhat(s*): mean %.2f sd %.2f ; gap mean %.2f (se %.2f) ; theoretical truncation bias >= %.3f\n",
              mean(f1), sd(f1), mean(f1) - (n - 2) * log(s_str), mean(fs), sd(fs), mean(fs) - (mean(f1) - (n - 2) * log(s_str)),
              sqrt(sd(fs)^2/length(fs) + sd(f1)^2/length(f1)), min(bound)))
}

cat("\n=== 5. cem and gam methods on the rlineage tree rescaled by 1e7 (natural-unit fit shown first) ===\n")
lb <- c(0.001, 0); ub <- c(1, 0.5)
run <- function(meth, bb, l, u, ctrl) { t0 <- Sys.time()
  ws <- character(0)
  f <- withCallingHandlers(tryCatch(estimate_rates(bb, method = meth, model = "cr", init_pars = c(l[1] * 20, l[1] * 5), control = c(list(lower_bound = l, upper_bound = u, max_time = 90), ctrl)),
       error = function(e) { cat("   ", meth, "ERROR:", conditionMessage(e), "\n"); NULL }), warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") })
  cat(sprintf("   %s (%.0fs): pars=%s loglik=%s warnings=%d %s\n", meth, as.numeric(Sys.time() - t0, units = "secs"),
              if (is.null(f)) "NULL" else paste(signif(f$pars, 3), collapse = ","), if (is.null(f)) "NULL" else format(f$loglik, digits = 6), length(ws),
              if (length(ws)) paste0("[", substr(ws[1], 1, 90), "...]") else "")) }
s <- 1e7
run("cem", b, lb, ub, list(num_particles = 10L, max_iter = 2L, maxN = 10L, num_trees = 1L))
run("cem", b * s, lb / s, ub / s, list(num_particles = 10L, max_iter = 2L, maxN = 10L, num_trees = 1L))
run("gam", b, lb, ub, list(sample_size = 5L, maxN = 50L, max_iter = 2L))
run("gam", b * s, lb / s, ub / s, list(sample_size = 5L, maxN = 50L, max_iter = 2L))
run("mcem", b * s, lb / s, ub / s, list(sample_size = 5L, maxN = 50L, max_iter = 2L))
run("mcem", b * s, lb / s, ub / s, list(sampling = "dynamic_fresh", sample_size = 5L, maxN = 50L, max_iter = 2L))

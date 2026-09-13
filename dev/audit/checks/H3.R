## H3 — CEM shared-tree mode (Mode 2) evaluates the IS denominator at the
## target particle theta_j instead of at the density that drew the tree.
##
## Part A: per-particle fhat, Mode 1 vs Mode 2 vs mixture-corrected Mode 2,
##         against DDD::bd_loglik (theta-independent offset removed).
## Part B: small end-to-end emphasis_cem Mode 1 vs Mode 2 vs DDD::bd_ML.
##
## Self-contained; ~1-2 min.  Rscript H3.R [A|B|AB]  (default AB)
PART <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "AB"

.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({ library(emphasis); library(ape); library(DDD) })

set.seed(3)                                  # R-side only (tree); C++ RNG is clock-seeded
phy  <- ape::rphylo(n = 20, birth = 0.5, death = 0.2)
brts <- sort(as.numeric(ape::branching.times(phy)), decreasing = TRUE)
cat(sprintf("tree: %d tips, crown age %.3f\n", length(brts) + 1, brts[1]))

model <- c(0L, 0L, 0L); link <- 0L; rho <- 1
pars8 <- function(lam, mu) c(lam, 0, 0, 0, mu, 0, 0, 0)
ref   <- function(lam, mu)                  # Nee crown-age lik, no conditioning, brts
  DDD::bd_loglik(pars1 = c(lam, mu, 0, 0), pars2 = c(0, 0, 1, 2, 0),
                 brts = brts, missnumspec = 0)

lme <- function(lw) { m <- max(lw); log(mean(exp(lw - m))) + m }

## particle set: same lambda, mu spread so the proposals differ
grid <- data.frame(lam = 0.5, mu = c(0.05, 0.15, 0.25, 0.35, 0.45))
K    <- nrow(grid)
N    <- 1500L                                # trees per particle
R    <- 4L                                   # replicate pools

draw <- function(k) {
  raw <- emphasis:::.simulate_particle(brts, pars8(grid$lam[k], grid$mu[k]), model, link,
                                       sample_size = N, maxN = 200L * N,
                                       max_missing = 1e4, max_lambda = 1e6,
                                       num_threads = 1L, rho = rho)
  stopifnot(!is.null(raw), length(raw$logf) == N, raw$rejected_zero_weights == 0)
  raw
}

if (grepl("A", PART)) {
res <- array(NA_real_, c(R, K, 3), dimnames = list(NULL, NULL, c("mode1", "mode2", "mode2_mix")))
consist <- numeric(0)
for (r in seq_len(R)) {
  raws  <- lapply(seq_len(K), draw)
  pool  <- do.call(c, lapply(raws, `[[`, "trees"))
  who   <- rep(seq_len(K), each = N)
  ## log q(z | theta_k) for every pooled tree at every particle (K x |pool|)
  LG <- LF <- matrix(NA_real_, K, length(pool))
  for (j in seq_len(K)) {
    ev <- emphasis:::eval_logf(pars8(grid$lam[j], grid$mu[j]), pool, model, link, rho)
    LF[j, ] <- ev$logf; LG[j, ] <- ev$logg
  }
  ## sanity: eval_logf's logg at the drawing particle == augment_trees' logg
  for (k in seq_len(K))
    consist <- c(consist, max(abs(LG[k, who == k] - raws[[k]]$logg)),
                          max(abs(LF[k, who == k] - raws[[k]]$logf)))
  for (j in seq_len(K)) {
    res[r, j, "mode1"]     <- lme(LF[j, who == j] - LG[j, who == j])         # .eval_independent
    res[r, j, "mode2"]     <- lme(LF[j, ] - LG[j, ])                          # .eval_shared (de.R:340-357)
    mix                    <- apply(LG, 2, function(v) { m <- max(v); log(mean(exp(v - m))) + m })
    res[r, j, "mode2_mix"] <- lme(LF[j, ] - mix)                              # proposed fix
  }
}
cat(sprintf("max |eval_logf - augment_trees| on logf/logg at the drawing particle: %.2e\n", max(consist)))

refv <- mapply(ref, grid$lam, grid$mu)
mn   <- apply(res, c(2, 3), mean); se <- apply(res, c(2, 3), sd) / sqrt(R)
## remove the theta-independent labelling offset using Mode 1 at particle 1
off  <- mn[1, "mode1"] - refv[1]
out  <- data.frame(mu = grid$mu, ref = refv - refv[1],
                   mode1 = mn[, "mode1"] - off - refv[1],
                   mode2 = mn[, "mode2"] - off - refv[1],
                   mode2_mix = mn[, "mode2_mix"] - off - refv[1])
out$bias_mode2 <- mn[, "mode2"] - mn[, "mode1"]
out$se_mode2   <- se[, "mode2"]; out$se_mode1 <- se[, "mode1"]
cat("\nPart A: relative loglik (offset removed, particle 1 = 0). ref = DDD::bd_loglik\n")
print(round(out, 3))
cat(sprintf("\noffset mode1 - bd_loglik at each particle (should be constant): %s\n",
            paste(round(mn[, "mode1"] - refv, 3), collapse = " ")))
cat(sprintf("offset mode2 - bd_loglik at each particle: %s\n",
            paste(round(mn[, "mode2"] - refv, 3), collapse = " ")))
cat(sprintf("offset mode2_mix - bd_loglik at each particle: %s\n",
            paste(round(mn[, "mode2_mix"] - refv, 3), collapse = " ")))
cat(sprintf("argmax over particles: ref %d, mode1 %d, mode2 %d, mode2_mix %d\n",
            which.max(refv), which.max(mn[, "mode1"]), which.max(mn[, "mode2"]),
            which.max(mn[, "mode2_mix"])))
}

## ---------------------------------------------------------------------------
## Part B: end-to-end CEM, Mode 1 vs Mode 2, vs DDD::bd_ML
if (grepl("B", PART)) {
cat("\nPart B: emphasis_cem end to end (4 replicates per mode)\n")
ml <- DDD::bd_ML(brts = brts, initparsopt = c(0.5, 0.2), idparsopt = 1:2,
                 parsfix = c(0, 0), idparsfix = 3:4, cond = 0, btorph = 1, soc = 2,
                 verbose = FALSE)
cat(sprintf("DDD::bd_ML: lambda %.3f mu %.3f loglik %.3f\n", ml$lambda0, ml$mu0, ml$loglik))
fit_cem <- function(shared) {
  fit <- estimate_rates(brts, method = "cem", model = "cr",
                        control = list(lower_bound = c(0.05, 0.0), upper_bound = c(1.0, 0.5),
                                       num_particles = 20L, num_trees = 10L, max_iter = 12L,
                                       maxN = 2000L,
                                       shared_trees = shared, max_time = 90, num_threads = 1L))
  c(fit$pars, loglik = fit$loglik)
}
for (sh in c(FALSE, TRUE)) {
  m <- t(replicate(4, fit_cem(sh)))
  m <- rbind(m, mean = colMeans(m), sd = apply(m, 2, sd))
  cat(sprintf("shared_trees = %s\n", sh)); print(round(m, 3))
}
}

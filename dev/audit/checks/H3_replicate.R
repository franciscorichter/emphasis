## H3 replication — vary what the verifier did not:
##   * a different tree (30 tips, seed 11, birth 0.4, death 0.1)
##   * particles that vary lambda (mu fixed) instead of mu
##   * a dd-model (model = c(1,0,0)) particle set as well
##   * the cheap fix from the sketch (cached drawing density) alongside the mixture
## Scores identical pooled trees four ways; mode1 is the reference (.eval_independent).
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({ library(emphasis); library(ape); library(DDD) })

set.seed(11)
phy  <- ape::rphylo(n = 30, birth = 0.4, death = 0.1)
brts <- sort(as.numeric(ape::branching.times(phy)), decreasing = TRUE)
cat(sprintf("tree: %d tips, crown age %.3f\n", length(brts) + 1, brts[1]))
lme <- function(lw) { m <- max(lw); log(mean(exp(lw - m))) + m }

run_set <- function(label, model, P, N = 1000L, R = 4L, ref = NULL) {
  K <- nrow(P)
  draw <- function(k) {
    raw <- emphasis:::.simulate_particle(brts, as.numeric(P[k, ]), model, 0L,
                                         sample_size = N, maxN = 200L * N,
                                         max_missing = 1e4, max_lambda = 1e6,
                                         num_threads = 1L, rho = 1)
    stopifnot(!is.null(raw), length(raw$logf) == N)
    raw
  }
  modes <- c("mode1", "mode2", "fix_cached", "fix_mix")
  res <- array(NA_real_, c(R, K, 4), dimnames = list(NULL, NULL, modes))
  nz <- 0L
  for (r in seq_len(R)) {
    raws <- lapply(seq_len(K), draw)
    nz   <- nz + sum(sapply(raws, `[[`, "rejected_zero_weights"))
    pool <- do.call(c, lapply(raws, `[[`, "trees"))
    who  <- rep(seq_len(K), each = N)
    logq_draw <- do.call(c, lapply(raws, `[[`, "logg"))     # pop$log_q, cached but unused
    LG <- LF <- matrix(NA_real_, K, length(pool))
    for (j in seq_len(K)) {
      ev <- emphasis:::eval_logf(as.numeric(P[j, ]), pool, model, 0L, 1)
      LF[j, ] <- ev$logf; LG[j, ] <- ev$logg
    }
    mix <- apply(LG, 2, function(v) { m <- max(v); log(mean(exp(v - m))) + m })
    for (j in seq_len(K)) {
      res[r, j, "mode1"]      <- lme(LF[j, who == j] - LG[j, who == j])
      res[r, j, "mode2"]      <- lme(LF[j, ] - LG[j, ])
      res[r, j, "fix_cached"] <- lme(LF[j, ] - logq_draw)
      res[r, j, "fix_mix"]    <- lme(LF[j, ] - mix)
    }
  }
  mn <- apply(res, c(2, 3), mean); se <- apply(res, c(2, 3), sd) / sqrt(R)
  cat(sprintf("\n== %s  (N=%d per particle, R=%d pools, zero-weight rejections=%d) ==\n",
              label, N, R, nz))
  out <- data.frame(P[, c(1, 5)], mode1 = mn[, "mode1"], se1 = se[, "mode1"],
                    bias_mode2 = mn[, "mode2"] - mn[, "mode1"], se2 = se[, "mode2"],
                    bias_cached = mn[, "fix_cached"] - mn[, "mode1"], se_c = se[, "fix_cached"],
                    bias_mix = mn[, "fix_mix"] - mn[, "mode1"], se_m = se[, "fix_mix"])
  names(out)[1:2] <- c("lam", "mu")
  if (!is.null(ref)) out$ref_rel <- ref - ref[1]
  print(round(out, 3))
  cat(sprintf("argmax: mode1 %d, mode2 %d, fix_cached %d, fix_mix %d%s\n",
              which.max(mn[, "mode1"]), which.max(mn[, "mode2"]),
              which.max(mn[, "fix_cached"]), which.max(mn[, "fix_mix"]),
              if (is.null(ref)) "" else sprintf(", DDD ref %d", which.max(ref))))
  invisible(res)
}

## Set 1: CR, vary lambda at fixed mu
P1 <- t(sapply(c(0.25, 0.35, 0.45, 0.55, 0.7), function(l) c(l, 0, 0, 0, 0.1, 0, 0, 0)))
ref1 <- sapply(P1[, 1], function(l)
  DDD::bd_loglik(pars1 = c(l, 0.1, 0, 0), pars2 = c(0, 0, 1, 2, 0), brts = brts, missnumspec = 0))
run_set("CR, vary lambda (mu = 0.1)", c(0L, 0L, 0L), P1, ref = ref1)

## Set 2: DD (model = c(1,0,0), linear link): lambda = b0 + bN*N, vary bN
P2 <- t(sapply(c(0, -0.004, -0.008, -0.012), function(bn) c(0.6, bn, 0, 0, 0.1, 0, 0, 0)))
run_set("DD linear, vary beta_N (b0 = 0.6, mu = 0.1)", c(1L, 0L, 0L), P2, N = 800L)

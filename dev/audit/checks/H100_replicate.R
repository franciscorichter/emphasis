## H100 replication: vary tree (20 tips, seed 11), theta grid, rho values {1, .6, .35},
## and add the exponential link (link = 1) as a second link.  Self-contained.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({ library(emphasis); library(ape); library(DDD) })
options(width = 130); t0 <- Sys.time()

p1_rho <- function(t, lam, mu, rho) {
  r <- lam - mu; E <- exp(-r * t); den <- rho * lam + (lam * (1 - rho) - mu) * E
  rho * r^2 * E / den^2
}
nee_rho <- function(lam, mu, rho, brts)
  2 * log(p1_rho(brts[1], lam, mu, rho)) + sum(log(lam) + log(p1_rho(brts[-1], lam, mu, rho)))

set.seed(11)
phy  <- ape::rphylo(20, birth = 0.4, death = 0.2)
brts <- sort(as.numeric(ape::branching.times(phy)), decreasing = TRUE)
n    <- length(brts) + 1L
cat(sprintf("Tree: n_tips = %d, crown age = %.4f\n", n, brts[1]))
grid <- rbind(c(0.6, 0.2), c(0.3, 0.25), c(0.5, 0.4))
off <- apply(grid, 1, function(p) nee_rho(p[1], p[2], 1, brts) -
               DDD::bd_loglik(c(p[1], p[2], 0, 0), c(0, 0, 1, 0, 2), brts, 0))
cat(sprintf("1a. nee_rho(rho=1) - DDD::bd_loglik: spread over 3 theta = %.2e\n", diff(range(off))))

pars8 <- function(l, m, link = 0L) if (link == 0L) c(l, 0, 0, 0, m, 0, 0, 0) else c(log(l), 0, 0, 0, log(m), 0, 0, 0)
aug <- function(brts, l, m, N, rho, link = 0L, maxN = 20L * N)
  emphasis:::augment_trees(brts, pars8(l, m, link), as.integer(N), as.integer(maxN), 100000L, 1e6, 1L,
                           c(0L, 0L, 0L), link, rho)
fhat_of <- function(r) { lw <- r$logf - r$logg; m <- max(lw)
  log(sum(exp(lw - m))) + m - log(length(lw) + r$rejected_zero_weights) }
ess_of <- function(r) { w <- exp(r$logf - r$logg - max(r$logf - r$logg)); sum(w)^2 / sum(w^2) }

Cv <- apply(grid, 1, function(p)
  emphasis:::.augment_tree_bdi(brts, pars = p, model_bin = c(0L, 0L, 0L), sample_size = 30L,
                               link = 0L, rho = 1.0)$fhat - nee_rho(p[1], p[2], 1, brts))
cat(sprintf("2. C (BDI fhat - nee at rho=1) per theta: %s\n", paste(sprintf("%.2e", Cv), collapse = " ")))
C <- mean(Cv)

reps <- 6L
res <- list()
for (rho in c(1, 0.6, 0.35)) for (i in seq_len(nrow(grid))) {
  p <- grid[i, ]; N <- if (rho < 0.5) 2000L else 1000L
  g <- t(sapply(seq_len(reps), function(k) {
    r <- aug(brts, p[1], p[2], N, rho)
    c(fhat = fhat_of(r), ess = ess_of(r), zero = r$rejected_zero_weights,
      n_unsamp = mean(sapply(r$trees, function(tr) sum(tr$t_ext == 5e10))))
  }))
  res[[length(res) + 1]] <- data.frame(rho = rho, lam = p[1], mu = p[2], N = N,
    gap = mean(g[, "fhat"]) - (nee_rho(p[1], p[2], rho, brts) + C),
    se = sd(g[, "fhat"]) / sqrt(reps), ess = mean(g[, "ess"]), zero_w = mean(g[, "zero"]),
    unsamp = mean(g[, "n_unsamp"]))
}
res <- do.call(rbind, res)
print(res, digits = 4, row.names = FALSE)
g1 <- res[res$rho == 1, ]; g0 <- res[res$rho < 1, ]
d <- merge(g0, g1, by = c("lam", "mu"), suffixes = c("", "_r1"))
d$z <- (d$gap - d$gap_r1) / sqrt(d$se^2 + d$se_r1^2)
print(d[, c("rho", "lam", "mu", "gap", "gap_r1", "z")], digits = 3, row.names = FALSE)
cat(sprintf("   max |z| = %.2f; for scale n*log(0.6) = %.2f, n*log(0.35) = %.2f\n",
            max(abs(d$z)), n * log(0.6), n * log(0.35)))

## other link: exponential, one cell, rho = 0.6 -> must match linear-link fhat within IS error
cat("\n2c. link = 1 (exponential) vs link = 0 at rho = 0.6, lam .6 mu .2, N = 1000, 6 reps\n")
f0 <- replicate(6, fhat_of(aug(brts, 0.6, 0.2, 1000L, 0.6, 0L)))
f1 <- replicate(6, fhat_of(aug(brts, 0.6, 0.2, 1000L, 0.6, 1L)))
cat(sprintf("   link0 %.4f (se %.4f) | link1 %.4f (se %.4f) | z = %.2f | ref+C = %.4f\n",
            mean(f0), sd(f0)/sqrt(6), mean(f1), sd(f1)/sqrt(6),
            (mean(f1)-mean(f0))/sqrt(var(f0)/6 + var(f1)/6), nee_rho(0.6, 0.2, 0.6, brts) + C))

## 3. pin the identities on this tree at rho = 0.35, both links
for (link in c(0L, 1L)) {
  rho <- 0.35; lam <- 0.5; mu <- 0.4
  r <- aug(brts, lam, mu, 200L, rho, link)
  e_r <- emphasis:::eval_logf(pars8(lam, mu, link), r$trees, model = c(0L, 0L, 0L), link = link, rho = rho)
  e_1 <- emphasis:::eval_logf(pars8(lam, mu, link), r$trees, model = c(0L, 0L, 0L), link = link, rho = 1.0)
  pred <- sapply(r$trees, function(tr) n * log(rho) + sum(tr$t_ext == 5e10) * log(1 - rho))
  cat(sprintf("3a. link %d: logf(rho)-logf(1) vs n log rho + n_unsamp log(1-rho): max|diff| = %.2e; unsamp/tree %s; any t_ext==5e10 at rho=1? ",
              link, max(abs((e_r$logf - e_1$logf) - pred)), paste(range(sapply(r$trees, function(tr) sum(tr$t_ext == 5e10))), collapse = "..")))
  r1 <- aug(brts, lam, mu, 200L, 1.0, link)
  cat(sprintf("%s\n", any(sapply(r1$trees, function(tr) any(tr$t_ext == 5e10)))))
  hand_logg <- function(tr) {
    o <- order(tr$brts); tr <- tr[o, ]
    Tt <- tr$brts[nrow(tr)]; prev <- 0; inte <- 0; logg <- 0; tips <- tr$n[1]; Ne <- 0
    for (i in seq_len(nrow(tr))) {
      s <- tr$brts[i]; dt <- s - prev
      inte <- inte + tr$n[i] * lam * (dt - (rho / mu) * (exp(-mu * (Tt - s)) - exp(-mu * (Tt - prev))))
      is_ext <- tr$t_ext[i] == 0; is_uns <- tr$t_ext[i] == 5e10; is_tip <- tr$t_ext[i] == 1e11
      is_mis <- !(is_ext || is_uns || is_tip)
      tips <- tips + is_tip; Ne <- Ne - is_ext
      if (is_mis) { logg <- logg + log(tr$n[i] * mu * lam) - mu * (tr$t_ext[i] - s) - log(2 * tips + Ne); Ne <- Ne + 1 }
      else if (is_uns) { logg <- logg + log(tr$n[i] * lam * (1 - rho)) - mu * (Tt - s) - log(2 * tips + Ne); Ne <- Ne + 1 }
      prev <- s
    }
    logg - inte
  }
  cat(sprintf("3b. link %d: hand logg vs returned: max|diff| = %.2e\n", link, max(abs(sapply(r$trees, hand_logg) - r$logg))))
}
cat(sprintf("\nelapsed %.0f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

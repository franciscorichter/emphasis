.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({ library(emphasis); library(ape) })
p1_rho <- function(t, lam, mu, rho) { r <- lam - mu; E <- exp(-r * t); den <- rho * lam + (lam * (1 - rho) - mu) * E; rho * r^2 * E / den^2 }
nee_rho <- function(lam, mu, rho, brts) 2 * log(p1_rho(brts[1], lam, mu, rho)) + sum(log(lam) + log(p1_rho(brts[-1], lam, mu, rho)))
set.seed(11); phy <- ape::rphylo(20, birth = 0.4, death = 0.2)
brts <- sort(as.numeric(ape::branching.times(phy)), decreasing = TRUE)
pars8 <- function(l, m, link) if (link == 0L) c(l,0,0,0,m,0,0,0) else c(log(l),0,0,0,log(m),0,0,0)
aug <- function(l, m, N, rho, link) emphasis:::augment_trees(brts, pars8(l,m,link), as.integer(N), 20L*as.integer(N), 100000L, 1e6, 1L, c(0L,0L,0L), link, rho)
fhat_of <- function(r) { lw <- r$logf - r$logg; m <- max(lw); log(sum(exp(lw - m))) + m - log(length(lw) + r$rejected_zero_weights) }
ess_of <- function(r) { w <- exp(r$logf - r$logg - max(r$logf - r$logg)); sum(w)^2 / sum(w^2) }
for (rho in c(1, 0.6)) {
  a <- t(replicate(10, { r <- aug(0.6, 0.2, 3000L, rho, 0L); c(fhat_of(r), ess_of(r)) }))
  b <- t(replicate(10, { r <- aug(0.6, 0.2, 3000L, rho, 1L); c(fhat_of(r), ess_of(r)) }))
  ref <- nee_rho(0.6, 0.2, rho, brts)
  cat(sprintf("rho %.1f N 3000 x10: link0 gap %.4f (se %.4f, ess %.0f) | link1 gap %.4f (se %.4f, ess %.0f) | z(link1-link0) = %.2f\n",
      rho, mean(a[,1])-ref, sd(a[,1])/sqrt(10), mean(a[,2]), mean(b[,1])-ref, sd(b[,1])/sqrt(10), mean(b[,2]),
      (mean(b[,1])-mean(a[,1]))/sqrt(var(a[,1])/10+var(b[,1])/10)))
}

.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths())); library(emphasis)
f <- emphasis:::.bdi_integral_cr
t1 <- 1; t2 <- 3; n <- 2; k <- 3; tp <- 5; lam <- 0.5
# closed-form critical limit (derived: p = 1/(1+lam*tau), 1-p = lam*tau/(1+lam*tau))
a1 <- 1 + lam*(tp-t1); a2 <- 1 + lam*(tp-t2)
I_lam_lim <- lam*(t2-t1) - log(a1/a2); I_mu_lim <- lam*(t2-t1) + log((tp-t1)/(tp-t2))
crit_lim <- (n+2*k)*I_lam_lim + n*I_mu_lim
# numeric check of the limit formula with integrate at exact criticality
pf <- function(t) 1/(1 + lam*(tp - t))
num_crit <- (n+2*k)*integrate(function(t) lam*(1-pf(t)), t1, t2, rel.tol=1e-12)$value + n*integrate(function(t) lam/(1-pf(t)), t1, t2, rel.tol=1e-12)$value
cat(sprintf("critical limit closed form = %.10f ; numeric integrate at d=0 = %.10f ; code critical branch = %.10f\n", crit_lim, num_crit, f(t1,t2,n,k,lam,lam,tp)))
for (d in c(1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-13, 1e-14, 1e-15, 0)) {
  mu <- lam - d
  code <- f(t1, t2, n, k, lam, mu, tp)
  cat(sprintf("d=%6.0e  code=%.10f  |code - critical limit|=%.2e\n", d, code, abs(code - crit_lim)))
}
# reachability from the sampler with lam == mu exactly
set.seed(4); brts <- sort(c(5, runif(9, 0, 5)), decreasing = TRUE)
aug <- tryCatch(emphasis:::.augment_tree_bdi(brts, pars = c(0.5, 0.5), model_bin = c(0L,0L,0L), sample_size = 50L, link = 0L, rho = 1), error = function(e) e)
if (inherits(aug, "error")) print(aug) else cat("BDI at lam==mu: trees =", length(aug$trees), " mean logg =", mean(aug$logg), " fhat =", aug$fhat, " sd(lw) =", sd(aug$logf-aug$logg), "\n")
aug2 <- emphasis:::.augment_tree_bdi(brts, pars = c(0.5, 0.5 - 1e-9), model_bin = c(0L,0L,0L), sample_size = 50L, link = 0L, rho = 1)
cat("BDI at lam-mu=1e-9: mean logg =", mean(aug2$logg), " fhat =", aug2$fhat, " sd(lw) =", sd(aug2$logf-aug2$logg), "\n")
aug3 <- emphasis:::.augment_tree_bdi(brts, pars = c(0.5, 0.5 - 1e-3), model_bin = c(0L,0L,0L), sample_size = 50L, link = 0L, rho = 1)
cat("BDI at lam-mu=1e-3: mean logg =", mean(aug3$logg), " fhat =", aug3$fhat, " sd(lw) =", sd(aug3$logf-aug3$logg), "\n")

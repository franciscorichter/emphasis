# Replicate the P(no missing) check to see whether the post-fix sampler is
# unbiased or whether a residual bias remains (and how tight 3se really is).
args <- commandArgs(TRUE)
.libPaths(c(args[1], .libPaths()))
suppressMessages(library(emphasis))
augment_trees <- get("augment_trees", envir = asNamespace("emphasis"))
brts7 <- c(6, 4.5, 3.0, 2.0, 1.2, 0.5); T7 <- brts7[1]; s7 <- T7 - brts7[-1]
is_missing_row <- function(df) !(df$t_ext == 0 | df$t_ext == 1e11 | df$t_ext == 5e10)
cr_hazard <- function(lam, mu) {
  knots <- c(0, s7, T7); nvec <- 2 + seq_along(knots) - 1; h <- 0
  for (i in seq_len(length(knots)-1)) { a<-knots[i]; b<-knots[i+1]
    h <- h + nvec[i]*lam*((b-a) - (1/mu)*(exp(-mu*(T7-b)) - exp(-mu*(T7-a)))) }
  h }
lam <- 0.4; mu <- 0.2; n <- 10000L
p0c <- exp(-cr_hazard(lam, mu))
zs <- numeric(6)
for (r in 1:6) {
  raw <- augment_trees(brts=brts7, pars=c(lam,0,0,0,mu,0,0,0), sample_size=n,
                       maxN=50L*n, max_missing=200L, max_lambda=1e6, num_threads=1L,
                       model=c(0L,0L,0L), link=0L, rho=1)
  nm <- vapply(raw$trees, function(df) sum(is_missing_row(df)), numeric(1))
  p0 <- mean(nm == 0); se <- sqrt(p0*(1-p0)/n)
  zs[r] <- (p0 - p0c)/se
  cat(sprintf("rep %d  p0=%.5f  z=%+.2f  mean_missing=%.4f\n", r, p0, zs[r], mean(nm)))
}
cat(sprintf("mean z = %+.2f  (|z|>3 in %d/6)\n", mean(zs), sum(abs(zs)>3)))

# Substantive part of test-thinning-envelope.R with the counter calls removed,
# so it runs on both builds. Usage: Rscript 1.8_substantive.R <lib>
args <- commandArgs(TRUE)
.libPaths(c(args[1], .libPaths()))
suppressMessages(library(emphasis))
augment_trees <- get("augment_trees", envir = asNamespace("emphasis"))

brts7 <- c(6, 4.5, 3.0, 2.0, 1.2, 0.5)
T7 <- brts7[1]
s7 <- T7 - brts7[-1]
t_ext_tip <- 1e11; t_ext_extinct <- 0; t_ext_unsampled <- 5e10
is_missing_row <- function(df) !(df$t_ext == t_ext_extinct | df$t_ext == t_ext_tip | df$t_ext == t_ext_unsampled)

cr_hazard <- function(lam, mu) {
  knots <- c(0, s7, T7); nvec <- 2 + seq_along(knots) - 1; h <- 0
  for (i in seq_len(length(knots) - 1)) {
    a <- knots[i]; b <- knots[i + 1]
    h <- h + nvec[i] * lam * ((b - a) - (1/mu) * (exp(-mu*(T7-b)) - exp(-mu*(T7-a))))
  }
  h
}
draw_cr <- function(lam, mu, n, link = 0L) {
  pars <- if (link == 1L) c(log(lam),0,0,0,log(mu),0,0,0) else c(lam,0,0,0,mu,0,0,0)
  augment_trees(brts = brts7, pars = pars, sample_size = n, maxN = 50L*n,
                max_missing = 200L, max_lambda = 1e6, num_threads = 1L,
                model = c(0L,0L,0L), link = link, rho = 1)
}
is_fhat <- function(logf, logg, n_zero = 0) {
  lw <- logf - logg; m <- max(lw); w <- exp(lw - m)
  c(fhat = log(mean(w)) + m - log(1 + n_zero/length(w)), se = sd(w)/(mean(w)*sqrt(length(w))))
}

set.seed(1)
# --- part 1
n <- 10000L
raw <- draw_cr(0.4, 0.2, n)
n_missing <- vapply(raw$trees, function(df) sum(is_missing_row(df)), numeric(1))
p0 <- mean(n_missing == 0); se <- sqrt(p0*(1-p0)/n); p0c <- exp(-cr_hazard(0.4, 0.2))
cat(sprintf("PART1  p0=%.5f closed=%.5f |diff|=%.5f 3se=%.5f  -> %s\n",
            p0, p0c, abs(p0-p0c), 3*se, ifelse(abs(p0-p0c) < 3*se, "PASS", "FAIL")))
logg_empty <- raw$logg[n_missing == 0]
cat(sprintf("PART1b logg_empty max dev from -H = %.3e -> %s\n",
            max(abs(logg_empty + cr_hazard(0.4,0.2))),
            ifelse(max(abs(logg_empty + cr_hazard(0.4,0.2))) < 1e-8, "PASS", "FAIL")))

# --- part 2
grid <- rbind(c(0.3,0.1), c(0.4,0.2), c(0.5,0.3), c(0.6,0.45), c(0.35,0.3))
n <- 8000L; z <- numeric(nrow(grid))
for (k in seq_len(nrow(grid))) {
  lam <- grid[k,1]; mu <- grid[k,2]
  raw <- draw_cr(lam, mu, n)
  f <- is_fhat(raw$logf, raw$logg, raw$rejected_zero_weights)
  ref <- DDD::bd_loglik(pars1 = c(lam,mu,0,0), pars2 = c(0,0,1,0,2), brts = brts7, missnumspec = 0)
  gap <- unname(f["fhat"] - ref); z[k] <- gap/unname(f["se"])
  cat(sprintf("PART2  lam=%.2f mu=%.2f gap=%+.4f tol=%.4f -> %s\n",
              lam, mu, gap, 0.03 + 3*unname(f["se"]),
              ifelse(abs(gap) < 0.03 + 3*unname(f["se"]), "PASS","FAIL")))
}
cat(sprintf("PART2 pooled sum(z^2)=%.2f  crit=%.2f -> %s\n", sum(z^2), qchisq(0.999,5),
            ifelse(sum(z^2) < qchisq(0.999,5), "PASS","FAIL")))

# --- part 3 (rejected_lambda only; counter unavailable pre-fix)
r1 <- draw_cr(0.4, 0.2, 2000L, link = 1L)
r2 <- augment_trees(brts=brts7, pars=c(0.6,-0.04,0,0,0.2,0,0,0), sample_size=2000L,
                    maxN=100000L, max_missing=200L, max_lambda=1e6, num_threads=1L,
                    model=c(1L,0L,0L), link=0L, rho=1)
r3 <- augment_trees(brts=brts7, pars=c(0.4,0,0.05,0,0.2,0,0.02,0), sample_size=2000L,
                    maxN=100000L, max_missing=200L, max_lambda=1e6, num_threads=1L,
                    model=c(0L,1L,0L), link=0L, rho=1)
cat(sprintf("PART3  rejected_lambda: exp=%d dd=%d M=%d\n",
            r1$rejected_lambda, r2$rejected_lambda, r3$rejected_lambda))

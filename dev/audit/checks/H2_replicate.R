## H2 replication (independent): BDI ignores rho.
## Varies what dev/audit/checks/H2.R did not: different tree, rho = 0.7,
## exponential link, the dd model through BDI, simulate_tree(method = "bdi"),
## and a DDD::bd_loglik cross-check of the closed form at rho = 1.
## Run: Rscript dev/audit/checks/H2_replicate.R
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({ library(emphasis); library(ape) })
options(width = 120)

set.seed(7)
phy  <- ape::rphylo(24, birth = 0.6, death = 0.2)
brts <- sort(ape::branching.times(phy), decreasing = TRUE)
n    <- length(brts) + 1L
cat(sprintf("Tree: n_tips = %d, crown age = %.4f\n", n, brts[1]))
rho <- 0.7

nee_rho <- function(lam, mu, rho, brts) {
  r <- lam - mu
  E <- exp(-r * brts)
  den <- rho * lam + (lam * (1 - rho) - mu) * E
  logp1 <- log(rho) + 2 * log(r) - r * brts - 2 * log(den)
  2 * logp1[1] + sum(log(lam) + logp1[-1])
}

## 0. Closed form vs DDD at rho = 1 (offset must be theta-independent)
if (requireNamespace("DDD", quietly = TRUE)) {
  grid0 <- rbind(c(0.6, 0.2), c(0.4, 0.3), c(0.9, 0.1))
  off <- apply(grid0, 1, function(p)
    nee_rho(p[1], p[2], 1, brts) -
      DDD::bd_loglik(pars1 = c(p[1], p[2], 0, 0), pars2 = c(0, 0, 1, 0, 2),
                     brts = brts, missnumspec = 0))
  cat(sprintf("0. nee(rho=1) - DDD::bd_loglik(cond=0,btorph=1,soc=2): %s (spread %.2e)\n",
              paste(round(off, 6), collapse = " "), diff(range(off))))
}

## 1. CR, linear link, rho = 0.7: BDI draws, unsampled count, shift
lam <- 0.6; mu <- 0.2
b1 <- emphasis:::.augment_tree_bdi(brts, pars = c(lam, mu), model_bin = c(0L,0L,0L),
                                   sample_size = 100L, link = 0L, rho = 1.0)
b7 <- emphasis:::.augment_tree_bdi(brts, pars = c(lam, mu), model_bin = c(0L,0L,0L),
                                   sample_size = 100L, link = 0L, rho = rho)
cat(sprintf("1. CR/linear: unsampled nodes per BDI tree at rho=0.7: max %d; lw sd rho=1: %.2e, rho=0.7: %.2e\n",
            max(sapply(b7$trees, function(tr) sum(tr$t_ext == 5e10))), sd(b1$weights), sd(b7$weights)))
cat(sprintf("   BDI fhat(rho=0.7) - fhat(rho=1) = %.6f ; n*log(rho) = %.6f\n",
            b7$fhat - b1$fhat, n * log(rho)))
cat(sprintf("   BDI fhat(rho=1) - nee(rho=1) = %.2e ; BDI fhat(rho=0.7) - nee(rho=0.7) = %.4f\n",
            b1$fhat - nee_rho(lam, mu, 1, brts), b7$fhat - nee_rho(lam, mu, rho, brts)))

thin_fhat <- function(p, rho, link = 0L, N = 1000L, reps = 4L) {
  sapply(seq_len(reps), function(i) {
    a <- emphasis:::.augment_tree_internal(brts, pars = p, model_bin = c(0L,0L,0L),
                                           sample_size = N, maxN = 50L * N,
                                           link = link, rho = rho, num_threads = 1L)
    lw <- a$logf - a$logg; m <- max(lw); log(mean(exp(lw - m))) + m
  })
}
th <- thin_fhat(c(lam, mu), rho)
cat(sprintf("   thinning fhat(rho=0.7): mean %.4f sd %.4f ; nee(rho=0.7) = %.4f ; BDI = %.4f\n",
            mean(th), sd(th), nee_rho(lam, mu, rho, brts), b7$fhat))

## 2. CR, exponential link (pars = log rates)
lp <- c(log(lam), log(mu))
e1 <- emphasis:::.augment_tree_bdi(brts, pars = lp, model_bin = c(0L,0L,0L),
                                   sample_size = 100L, link = 1L, rho = 1.0)
e7 <- emphasis:::.augment_tree_bdi(brts, pars = lp, model_bin = c(0L,0L,0L),
                                   sample_size = 100L, link = 1L, rho = rho)
the <- thin_fhat(lp, rho, link = 1L)
cat(sprintf("2. CR/exp link: unsampled max %d; BDI fhat(1) - nee(1) = %.2e; BDI fhat(0.7) - fhat(1) = %.6f (n log rho = %.6f)\n",
            max(sapply(e7$trees, function(tr) sum(tr$t_ext == 5e10))),
            e1$fhat - nee_rho(lam, mu, 1, brts), e7$fhat - e1$fhat, n * log(rho)))
cat(sprintf("   thinning fhat(0.7): mean %.4f sd %.4f ; nee(0.7) = %.4f ; BDI(0.7) = %.4f\n",
            mean(the), sd(the), nee_rho(lam, mu, rho, brts), e7$fhat))

## 3. dd model through BDI at rho = 0.7: does it ever emit unsampled nodes?
dd7 <- emphasis:::.augment_tree_bdi(brts, pars = c(0.8, -0.01, 0.2, 0.0), model_bin = c(1L,0L,0L),
                                    sample_size = 50L, link = 0L, rho = rho)
dd1 <- emphasis:::.augment_tree_bdi(brts, pars = c(0.8, -0.01, 0.2, 0.0), model_bin = c(1L,0L,0L),
                                    sample_size = 50L, link = 0L, rho = 1.0)
cat(sprintf("3. dd/BDI rho=0.7: %d trees, unsampled max %d; fhat(0.7) = %.4f, fhat(1) = %.4f, diff = %.4f (n log rho = %.4f)\n",
            length(dd7$trees), max(sapply(dd7$trees, function(tr) sum(tr$t_ext == 5e10))),
            dd7$fhat, dd1$fhat, dd7$fhat - dd1$fhat, n * log(rho)))
cat(sprintf("   .bdi_supported(dd, linear) = %s\n", emphasis:::.bdi_supported(c(1L,0L,0L), 0L)))

## 4. simulate_tree(method = "bdi", rho = 0.7) dispatch: any unsampled lineages in output?
st <- tryCatch(simulate_tree(brts, pars = c(lam, mu), model = "cr", n_trees = 20L,
                             method = "bdi", rho = rho),
               error = function(e) { cat("   simulate_tree error:", conditionMessage(e), "\n"); NULL })
if (!is.null(st)) {
  cat("4. simulate_tree(method='bdi', rho=0.7) returned names:", paste(names(st), collapse = ","), "\n")
  trs <- if (!is.null(st$trees)) st$trees else if (!is.null(st$tas)) st$tas else NULL
  if (is.list(trs) && length(trs) && is.data.frame(trs[[1]]))
    cat(sprintf("   unsampled nodes max over trees: %d\n",
                max(sapply(trs, function(tr) sum(tr$t_ext == 5e10)))))
}

## 5. Does estimate_rates say anything when rho < 1 with bdi? Capture messages.
msgs <- character(0)
f_bdi <- withCallingHandlers(
  estimate_rates(brts, method = "mcem", model = "cr", init_pars = c(0.4, 0.1),
                 control = list(rho = rho, sampling = "bdi", sample_size = 100L,
                                max_iter = 25L, max_time = 60, num_threads = 1L,
                                lower_bound = c(1e-3, 0), upper_bound = c(3, 3), verbose = TRUE)),
  message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
f_thin <- estimate_rates(brts, method = "mcem", model = "cr", init_pars = c(0.4, 0.1),
                 control = list(rho = rho, sampling = "dynamic_fresh", sample_size = 100L,
                                max_iter = 25L, max_time = 60, num_threads = 1L,
                                lower_bound = c(1e-3, 0), upper_bound = c(3, 3)))
mle <- function(rho) {
  o <- optim(c(0.4, 0.1), function(p) if (p[1] <= p[2] || p[2] < 0) 1e10 else -nee_rho(p[1], p[2], rho, brts),
             method = "L-BFGS-B", lower = c(1e-3, 0), upper = c(5, 5)); o$par }
m1 <- mle(1); m7 <- mle(rho)
cat(sprintf("5. closed-form MLE rho=1: (%.4f, %.4f) ; rho=0.7: (%.4f, %.4f)\n", m1[1], m1[2], m7[1], m7[2]))
cat(sprintf("   estimate_rates bdi   (rho=0.7): pars = (%.4f, %.4f) loglik %.4f\n", f_bdi$pars[1], f_bdi$pars[2], f_bdi$loglik))
cat(sprintf("   estimate_rates thin  (rho=0.7): pars = (%.4f, %.4f) loglik %.4f\n", f_thin$pars[1], f_thin$pars[2], f_thin$loglik))
cat("   messages mentioning rho/fallback in verbose bdi run:",
    sum(grepl("rho|fall", msgs, ignore.case = TRUE)), "of", length(msgs), "\n")

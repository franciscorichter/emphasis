## H10 replication — independent variation of dev/audit/checks/H10.R
##  * different tree (seed 5, dd_sim(c(0.6, 0.1, 18), age 8))
##  * the PACKAGE's own .augment_tree_bdi (not a re-implemented loop), with
##    rejections counted by trace() on .bdi_augment_one
##  * link = 1 (exponential) as the CR control, i.e. the other link
##  * smaller theta grid vs DDD::dd_loglik(cond = 0, btorph = 1, soc = 2)
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressMessages(library(emphasis))
suppressMessages(library(DDD))
options(width = 140)

set.seed(5)
sim  <- DDD::dd_sim(pars = c(0.6, 0.1, 18), age = 8)
phy  <- sim$tes
brts <- sort(ape::branching.times(phy), decreasing = TRUE)
cat(sprintf("Tree: %d tips, crown age %.3f\n\n", ape::Ntip(phy), brts[1]))

## count NULL returns of .bdi_augment_one inside the package function
.H10_att <- 0L; .H10_rej <- 0L
trace(".bdi_augment_one", where = asNamespace("emphasis"), print = FALSE,
      exit = quote({ .H10_att <<- .H10_att + 1L
                     if (is.null(returnValue())) .H10_rej <<- .H10_rej + 1L }))

run_pkg <- function(pars8, model_bin, link, sample_size, max_missing = 1e6) {
  .H10_att <<- 0L; .H10_rej <<- 0L
  r <- emphasis:::.augment_tree_bdi(brts, pars8, model_bin, sample_size = sample_size,
                                    max_missing = as.integer(max_missing), link = link, rho = 1)
  acc <- 1 - .H10_rej / .H10_att
  lw  <- r$weights; ok <- is.finite(lw)
  c(n_att = .H10_att, n_rej = .H10_rej, n_trees = length(r$trees), acc = acc,
    nonfin = sum(!ok), fhat_pkg = r$fhat, fhat_corr = r$fhat + log(acc),
    sd_lw = sd(lw[ok]))
}

## ---------- A. CR control on the other link ----------------------------
cat("== A. CR control (lambda 0.6, mu 0.1): cr/link0 exact vs dd/link1 with beta_N = gamma_N = 0 ==\n")
p_cr <- c(0.6, 0, 0, 0, 0.1, 0, 0, 0)
p_dd1 <- c(log(0.6), 0, 0, 0, log(0.1), 0, 0, 0)     # exponential link, same rates
set.seed(1)
cr <- run_pkg(p_cr, c(0L, 0L, 0L), 0L, 200L)
print(round(cr, 4))
set.seed(2)
A <- t(sapply(1:4, function(i) run_pkg(p_dd1, c(1L, 0L, 0L), 1L, 300L)))
print(round(A, 4))
cat(sprintf("fhat_dd(link1) - fhat_cr: mean %.4f ;  -log(acc): mean %.4f ;  corrected - cr: mean %.4f (sd %.4f)\n\n",
            mean(A[, "fhat_pkg"] - cr["fhat_pkg"]), mean(-log(A[, "acc"])),
            mean(A[, "fhat_corr"] - cr["fhat_pkg"]), sd(A[, "fhat_corr"] - cr["fhat_pkg"])))

## ---------- B. theta grid vs DDD::dd_loglik (linear link) ----------------
cat("== B. theta grid vs DDD::dd_loglik ==\n")
grid <- rbind(c(0.6, 0.1, 18), c(0.6, 0.1, 30), c(0.6, 0.1, 100),
              c(0.9, 0.3, 25), c(0.9, 0.3, 60), c(0.5, 0.05, 40), c(1.2, 0.6, 30))
set.seed(3)
B <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  l0 <- grid[i, 1]; m0 <- grid[i, 2]; K <- grid[i, 3]
  pars8 <- c(l0, -(l0 - m0) / K, 0, 0, m0, 0, 0, 0)
  ref <- DDD::dd_loglik(pars1 = c(l0, m0, K), pars2 = c(300, 1, 0, 1, 0, 2),
                        brts = brts, missnumspec = 0)
  reps <- t(sapply(1:2, function(r) run_pkg(pars8, c(1L, 0L, 0L), 0L, 300L)))
  data.frame(l0 = l0, m0 = m0, K = K, dd_loglik = ref,
             acc = mean(reps[, "acc"]), neg_log_acc = -log(mean(reps[, "acc"])),
             fhat_pkg = mean(reps[, "fhat_pkg"]), sd_rep = sd(reps[, "fhat_pkg"]),
             gap_raw = mean(reps[, "fhat_pkg"]) - ref,
             gap_corr = mean(reps[, "fhat_corr"]) - ref,
             nonfin = sum(reps[, "nonfin"]), sd_lw = mean(reps[, "sd_lw"]))
}))
print(B, digits = 4, row.names = FALSE)
Bc <- B[is.finite(B$dd_loglik) & B$nonfin == 0, ]
fit <- lm(gap_raw ~ neg_log_acc, data = Bc)
cat(sprintf("\n%d clean points: sd(gap_raw) = %.4f, sd(gap_corr) = %.4f; lm slope = %.3f, intercept = %.4f, R^2 = %.3f\n\n",
            nrow(Bc), sd(Bc$gap_raw), sd(Bc$gap_corr), coef(fit)[2], coef(fit)[1], summary(fit)$r.squared))

## ---------- C. max_missing channel via the package function ------------
cat("== C. max_missing channel (l0 0.9, mu 0.3, K 25) ==\n")
pars8 <- c(0.9, -(0.9 - 0.3) / 25, 0, 0, 0.3, 0, 0, 0)
set.seed(4)
C <- t(sapply(c(1e6, 8, 5, 3), function(mm) c(max_missing = mm, run_pkg(pars8, c(1L, 0L, 0L), 0L, 200L, max_missing = mm))))
print(round(C, 4))
cat("(n_trees < sample_size when max_tries = 5*sample_size exhausted; no rejection field returned)\n")
untrace(".bdi_augment_one", where = asNamespace("emphasis"))

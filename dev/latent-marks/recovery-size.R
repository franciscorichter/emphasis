# From how many tips can the latent rate be read?  The recovery study's model (cr with the trend-free
# latent rate, lambda0 0.35, mu0 0.07) at three tree sizes: crown age 14, 16.5 and 19 put the median
# clade at about 100, 200 and 400 extant tips, and a tree is kept when it lands within a factor 0.8 to
# 1.25 of that size.  Each tree is fitted by MCEM with the mark-sampling proposal, with the latent rate
# and with cr alone (whose log-likelihood is exact); 20 trees per sigma and size.
#   Rscript dev/latent-marks/recovery-size.R <task 1-120> <out_dir> [lib]
args <- commandArgs(trailingOnly = TRUE)
task <- as.integer(args[1]); out <- args[2]
if (length(args) >= 3 && nzchar(args[3])) .libPaths(c(args[3], .libPaths()))
suppressMessages(library(emphasis)); e <- asNamespace("emphasis")
GRID <- expand.grid(rep = 1:20, size = c(100, 200, 400), sigma = c(0.25, 0.5))
g <- GRID[task, ]; sg <- g$sigma; n <- g$size; rep <- g$rep
age <- c(`100` = 14, `200` = 16.5, `400` = 19)[[as.character(n)]]
mb <- e$.resolve_model("cr"); truth <- c(log(0.35), log(0.07))
mko <- function(x) { s <- exp(x[1]); e$.latent_rate_marks(s, log_alpha = -s^2 / 2, kmax = 20) }
one <- function(x) list(S = 1L, u = 0, log_clam = 0, log_cmu = 0, C = 1, P = 1, pi = 1)
set.seed(3e6 + 1e5 * (sg * 100) + 1e3 * (n / 100) + rep)
tries <- 0L; capped <- 0L
repeat {
  tries <- tries + 1L
  r <- e$.simulate_marks(truth, "cr", age, mko(log(sg)), max_N = 40000, max_tries = 1L)
  if (r$status == "too_large") capped <- capped + 1L
  if (r$status == "done") { nt <- sum(r$Ltable[, 4] == -1); if (nt >= 0.8 * n && nt <= 1.25 * n) break }
}
phy <- DDD::L2phylo(r$Ltable, dropextinct = TRUE)
t0 <- Sys.time()
fm <- e$.mcem_marks(phy, "cr", c(log(0.3), log(0.1)), log(0.3), mko, n_iter = 15L, N = 200L, aware = "sample")
tm <- as.numeric(Sys.time() - t0, units = "secs")
set.seed(11); b <- e$.pm_draws(phy, exp(fm$estimate[1]), exp(fm$estimate[2]), mko(fm$estimate[3]), 500L)
lw <- e$.eval_logf_marks(e$.expand_pars(fm$estimate[1:2], mb), b$trees, mb, mko(fm$estimate[3])) - b$logg
lw <- lw[is.finite(lw)]; ll_m <- max(lw) + log(mean(exp(lw - max(lw))))
fc <- e$.mcem_marks(phy, "cr", c(log(0.3), log(0.1)), numeric(0), one, n_iter = 15L, N = 200L, aware = "sample")
Bp <- e$.pm_branches(phy); l0 <- exp(fc$estimate[1]); m0 <- exp(fc$estimate[2])
ll_c <- e$.pm_prune(Bp, e$.pm_fac_obs(Bp, l0, m0), l0, matrix(0, 1, 1), 0)$root
dir.create(out, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(data.frame(task = task, sigma = sg, size = n, rep = rep, age = age, tips = ape::Ntip(phy),
  tries = tries, capped = capped, lambda0 = exp(fm$estimate[1]), mu0 = exp(fm$estimate[2]), sigma_hat = exp(fm$estimate[3]),
  ess_last = utils::tail(fm$trace[, "ess"], 1), lambda0_cr = l0, mu0_cr = m0, ll_marks = ll_m, ll_cr = ll_c,
  lr = 2 * (ll_m - ll_c), sec = tm, build = as.character(utils::packageVersion("emphasis"))),
  file.path(out, sprintf("task-%03d.csv", task)), row.names = FALSE)

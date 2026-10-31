# Which proposal limits the latent-rate E-step: the effective sample of 500 draws at the generating
# values, with and without the marks, under the thinning proposal (augment_trees, the one the
# prototype .mcem_marks uses) and the surrogate (.augment_tree_bdi, which samples the hidden lineages
# of cr exactly).  The trees are those of the recovery study (recovery-mcem.R: same seeds, same size
# window), so the rows join its per-tree results.
#   Rscript dev/latent-marks/proposal-ess-surrogate.R <task> <out_dir> [lib]
args <- commandArgs(trailingOnly = TRUE)
task <- as.integer(args[1]); out <- args[2]
if (length(args) >= 3 && nzchar(args[3])) .libPaths(c(args[3], .libPaths()))
suppressMessages(library(emphasis))
e <- asNamespace("emphasis")
SIGMAS <- c(0, 0.25, 0.5); REPS <- 30L
sg  <- SIGMAS[(task - 1L) %/% REPS + 1L]; rep <- (task - 1L) %% REPS + 1L
mb  <- e$.resolve_model("cr")
mko <- function(x) { s <- exp(x[1]); e$.latent_rate_marks(s, log_alpha = -s^2 / 2, kmax = 15) }
one <- function(x) list(S = 1L, u = 0, log_clam = 0, log_cmu = 0, C = 1, P = 1, pi = 1)
truth <- c(log(0.35), log(0.07))

set.seed(1e5 * (sg * 100) + rep)
repeat {
  r <- e$.simulate_marks(truth, "cr", 14, if (sg > 0) mko(log(sg)) else one(), max_N = 20000, max_tries = 1L)
  if (r$status == "done") { nt <- sum(r$Ltable[, 4] == -1); if (nt >= 50 && nt <= 200) break }
}
phy  <- DDD::L2phylo(r$Ltable, dropextinct = TRUE)
brts <- e$.extract_brts(phy)
pf   <- e$.expand_pars(truth, mb)
marks <- if (sg > 0) mko(log(sg)) else one()
ess <- function(lw) { lw <- lw[is.finite(lw)]; if (!length(lw)) return(0); w <- exp(lw - max(lw)); sum(w)^2 / sum(w^2) }
lme <- function(lw) { lw <- lw[is.finite(lw)]; m <- max(lw); m + log(mean(exp(lw - m))) }
N <- 500L

t0 <- Sys.time()
a <- e$augment_trees(brts, pf, N, 50L * N, 5000L, 1e6, 1L, model = mb, link = 1L, rho = 1,
                     parent_tip_start = e$.pts(brts), seed = 7L, parent_id = e$.pid(brts))
lm_t <- if (sg > 0) e$.eval_logf_marks(pf, a$trees, mb, marks) else a$logf
sec_t <- as.numeric(Sys.time() - t0, units = "secs")

t0 <- Sys.time()
set.seed(7)
b <- e$.augment_tree_bdi(phy, pars = truth, model_bin = mb, sample_size = N, link = 1L, rho = 1)
lm_s <- if (sg > 0) e$.eval_logf_marks(pf, b$trees, mb, marks) else b$logf
sec_s <- as.numeric(Sys.time() - t0, units = "secs")

dir.create(out, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(data.frame(
  task = task, sigma = sg, rep = rep, tips = ape::Ntip(phy), draws = N,
  ess_thin_nomarks = ess(a$logf - a$logg), ess_thin_marks = ess(lm_t - a$logg),
  ess_surr_nomarks = ess(b$logf - b$logg), ess_surr_marks = ess(lm_s - b$logg),
  surr_acc = b$acc, surr_valid = length(b$trees),
  # log-likelihood estimates with the marks, each proposal's own (the surrogate's with its rejections)
  ll_thin_marks = lme(lm_t - a$logg), ll_surr_marks = lme(lm_s - b$logg) + log(b$acc),
  sec_thin = sec_t, sec_surr = sec_s,
  build = as.character(utils::packageVersion("emphasis"))),
  file.path(out, sprintf("task-%03d.csv", task)), row.names = FALSE)

# The mark-sampling proposal (proposal-marks.R) on the trees of the recovery study: the effective
# sample of N draws at the generating values with the marks, and the log-likelihood estimate.
#   Rscript dev/latent-marks/proposal-ess-marks.R <task> <out_dir> [lib] [N = 200]
args <- commandArgs(trailingOnly = TRUE)
task <- as.integer(args[1]); out <- args[2]
if (length(args) >= 3 && nzchar(args[3])) .libPaths(c(args[3], .libPaths()))
N <- if (length(args) >= 4) as.integer(args[4]) else 200L
suppressMessages(library(emphasis)); e <- asNamespace("emphasis")
for (f in c(".ltable_complete_frame", ".eval_logf_marks", ".latent_rate_marks", ".simulate_marks", ".expand_pars", ".resolve_model", ".pm_draws"))
  assign(f, get(f, e))


SIGMAS <- c(0, 0.25, 0.5); REPS <- 30L
sg <- SIGMAS[(task - 1L) %/% REPS + 1L]; rep <- (task - 1L) %% REPS + 1L
mb <- .resolve_model("cr"); truth <- c(log(0.35), log(0.07))
mko <- function(s) .latent_rate_marks(s, log_alpha = -s^2 / 2, kmax = 15)
one <- list(S = 1L, u = 0, log_clam = 0, log_cmu = 0, C = 1, P = 1, pi = 1)
set.seed(1e5 * (sg * 100) + rep)
repeat { r <- .simulate_marks(truth, "cr", 14, if (sg > 0) mko(sg) else one, max_N = 20000, max_tries = 1L)
  if (r$status == "done") { nt <- sum(r$Ltable[, 4] == -1); if (nt >= 50 && nt <= 200) break } }
phy <- DDD::L2phylo(r$Ltable, dropextinct = TRUE); pf <- .expand_pars(truth, mb)
mk <- if (sg > 0) mko(sg) else one
set.seed(5); t0 <- Sys.time()
d <- .pm_draws(phy, exp(truth[1]), exp(truth[2]), mk, N)
lw <- .eval_logf_marks(pf, d$trees, mb, mk) - d$logg
sec <- as.numeric(Sys.time() - t0, units = "secs")
ok <- is.finite(lw); w <- exp(lw[ok] - max(lw[ok]))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(data.frame(task = task, sigma = sg, rep = rep, tips = ape::Ntip(phy), draws = N,
  ess_marks = sum(w)^2 / sum(w^2), sd_logw = stats::sd(lw[ok]), hidden_median = stats::median(d$n_hidden),
  ll_marks = max(lw[ok]) + log(mean(w)), sec = sec, build = as.character(utils::packageVersion("emphasis"))),
  file.path(out, sprintf("task-%03d.csv", task)), row.names = FALSE)

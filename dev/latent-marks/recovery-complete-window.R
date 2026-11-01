# Is the shortfall of sigma-hat in the recovery study the reconstruction's, the size window's, or the
# single tree's?  Per-tree ML of (lambda0, mu0, sigma) on COMPLETE trees (extinct lineages kept),
#   arm "window":   the recovery study's trees (same seeds, kept at 50-200 extant tips);
#   arm "nowindow": from the same seed, the first tree in which both crown lineages survive, any
#                   size up to the simulator's cap (5000 lineages; larger ones recorded as capped).
#   Rscript dev/latent-marks/recovery-complete-window.R <task> <arm> <out_dir> [lib]
args <- commandArgs(trailingOnly = TRUE)
task <- as.integer(args[1]); arm <- args[2]; out <- args[3]
if (length(args) >= 4 && nzchar(args[4])) .libPaths(c(args[4], .libPaths()))
suppressMessages(library(emphasis)); e <- asNamespace("emphasis")
SIGMAS <- c(0, 0.25, 0.5); REPS <- 30L
sg <- SIGMAS[(task - 1L) %/% REPS + 1L]; rep <- (task - 1L) %% REPS + 1L
mb <- e$.resolve_model("cr"); truth <- c(log(0.35), log(0.07))
mko <- function(s) e$.latent_rate_marks(s, log_alpha = -s^2 / 2, kmax = 15)
one <- list(S = 1L, u = 0, log_clam = 0, log_cmu = 0, C = 1, P = 1, pi = 1)
set.seed(1e5 * (sg * 100) + rep)
tries <- 0L; capped <- 0L
repeat {
  tries <- tries + 1L
  r <- e$.simulate_marks(truth, "cr", 14, if (sg > 0) mko(sg) else one,
                         max_N = if (arm == "window") 20000 else 5000, max_tries = 1L)
  if (r$status == "too_large") capped <- capped + 1L
  if (r$status != "done") next
  ext <- r$Ltable[, 4] == -1; nt <- sum(ext)
  both <- any(r$Ltable[ext, 3] < 0) && any(r$Ltable[ext, 3] > 0)
  if (arm == "window") { if (nt >= 50 && nt <= 200) break } else if (both) break
}
fr <- r$frame
nll <- function(th) {
  s <- exp(th[3]); mk <- if (s < 1e-4) one else mko(s)
  v <- -sum(e$.eval_logf_marks(e$.expand_pars(th[1:2], mb), list(fr), mb, mk))
  if (is.finite(v)) v else 1e10
}
t0 <- Sys.time()
f <- stats::optim(c(log(0.3), log(0.1), log(0.3)), nll, control = list(maxit = 300))
f0 <- stats::optim(c(log(0.3), log(0.1)), function(th) nll(c(th, log(1e-6))), control = list(maxit = 200))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(data.frame(task = task, arm = arm, sigma = sg, rep = rep, tips = nt, lineages = nrow(r$Ltable),
  tries = tries, capped_before = capped, lambda0 = exp(f$par[1]), mu0 = exp(f$par[2]), sigma_hat = exp(f$par[3]),
  lr = 2 * (f0$value - f$value), conv = f$convergence, sec = as.numeric(Sys.time() - t0, units = "secs"),
  build = as.character(utils::packageVersion("emphasis"))),
  file.path(out, sprintf("%s-task-%03d.csv", arm, task)), row.names = FALSE)

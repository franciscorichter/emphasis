# Phase 3 recovery study (dev/latent-marks/design.md), one tree per call:
#   Rscript dev/latent-marks/recovery-mcem.R <task> <out_dir> [lib] [proposal: none | lambda | both]
# Task t (1-based) is generating sigma SIGMAS[(t - 1) %/% REPS + 1], replicate (t - 1) %% REPS + 1.
# A tree: cr with the trend-free latent rate (log_alpha = -sigma^2 / 2), lambda0 0.35, mu0 0.07,
# crown age 14, redrawn until it has 50 to 200 tips (a conditioning on size, reported with the
# results). Two MCEM fits with the prototype (.mcem_marks): cr with the latent rate (sigma free)
# and cr alone; each ends with a log-likelihood estimate from 2000 draws at its estimate. Also
# recorded: the effective sample of 500 draws at the generating values, with and without marks.
# Writes <out_dir>/task-<t>.csv.
args <- commandArgs(trailingOnly = TRUE)
task <- as.integer(args[1]); out <- args[2]
if (length(args) >= 3 && nzchar(args[3])) .libPaths(c(args[3], .libPaths()))
PROP <- if (length(args) >= 4) args[4] else "none"
suppressMessages(library(emphasis))
e <- asNamespace("emphasis")
SIGMAS <- c(0, 0.25, 0.5); REPS <- 30L
sg  <- SIGMAS[(task - 1L) %/% REPS + 1L]; rep <- (task - 1L) %% REPS + 1L
mb  <- e$.resolve_model("cr")
mko <- function(x) { s <- exp(x[1]); e$.latent_rate_marks(s, log_alpha = -s^2 / 2, kmax = 15) }
one <- function(x) list(S = 1L, u = 0, log_clam = 0, log_cmu = 0, C = 1, P = 1, pi = 1)
truth <- c(log(0.35), log(0.07))

set.seed(1e5 * (sg * 100) + rep)
tries <- 0L
repeat {
  tries <- tries + 1L
  r <- e$.simulate_marks(truth, "cr", 14, if (sg > 0) mko(log(sg)) else one(), max_N = 20000, max_tries = 1L)
  if (r$status == "done") { nt <- sum(r$Ltable[, 4] == -1); if (nt >= 50 && nt <= 200) break }
}
phy  <- DDD::L2phylo(r$Ltable, dropextinct = TRUE)
brts <- e$.extract_brts(phy)

draws <- function(th, n, seed) {
  e$augment_trees(brts, e$.expand_pars(th, mb), as.integer(n), as.integer(50 * n), 5000L, 1e6, 1L,
                  model = mb, link = 1L, rho = 1, parent_tip_start = e$.pts(brts), seed = as.integer(seed),
                  parent_id = e$.pid(brts))
}
lme <- function(lw) { lw <- lw[is.finite(lw)]; m <- max(lw); m + log(mean(exp(lw - m))) }
ess <- function(lw) { lw <- lw[is.finite(lw)]; w <- exp(lw - max(lw)); sum(w)^2 / sum(w^2) }

a <- draws(truth, 500, 7)
ess_nomarks <- ess(a$logf - a$logg)
ess_marks   <- if (sg > 0) ess(e$.eval_logf_marks(e$.expand_pars(truth, mb), a$trees, mb, mko(log(sg))) - a$logg) else NA
ess_aware <- function(with_mu) {
  if (sg == 0) return(NA)
  o <- e$.mark_offsets(e$.expand_pars(truth, mb), phy, "cr", mko(log(sg)))
  a <- e$augment_trees(brts, e$.expand_pars(truth, mb), 500L, 25000L, 5000L, 1e6, 1L, model = mb, link = 1L, rho = 1,
                       parent_tip_start = e$.pts(brts), seed = 7L, parent_id = e$.pid(brts), offset_id = o$id,
                       offset_from = o$from, offset_to = o$to, offset = o$off,
                       offset_mu = if (with_mu) o$off_mu else numeric(0))
  ess(e$.eval_logf_marks(e$.expand_pars(truth, mb), a$trees, mb, mko(log(sg))) - a$logg)
}
ess_marks_lambda <- ess_aware(FALSE)
ess_marks_both   <- ess_aware(TRUE)

t0 <- Sys.time()
fm <- e$.mcem_marks(phy, "cr", c(log(0.3), log(0.1)), log(0.3), mko, n_iter = 15L, N = 200L, aware = PROP)
tm <- as.numeric(Sys.time() - t0, units = "secs")
b  <- draws(fm$estimate[1:2], 2000, 11)
ll_m <- lme(e$.eval_logf_marks(e$.expand_pars(fm$estimate[1:2], mb), b$trees, mb, mko(fm$estimate[3])) - b$logg)

t0 <- Sys.time()
fc <- e$.mcem_marks(phy, "cr", c(log(0.3), log(0.1)), numeric(0), one, n_iter = 15L, N = 200L)
tc <- as.numeric(Sys.time() - t0, units = "secs")
b  <- draws(fc$estimate[1:2], 2000, 13)
ll_c <- lme(b$logf - b$logg)

dir.create(out, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(data.frame(
  task = task, sigma = sg, rep = rep, tips = ape::Ntip(phy), lineages_ever = nrow(r$Ltable), draws_to_window = tries,
  ess_truth_nomarks = ess_nomarks, ess_truth_marks = ess_marks,
  ess_truth_marks_lambda = ess_marks_lambda, ess_truth_marks_both = ess_marks_both, proposal = PROP,
  lambda0_marks = exp(fm$estimate[1]), mu0_marks = exp(fm$estimate[2]), sigma_hat = exp(fm$estimate[3]),
  loglik_marks = ll_m, ess_last_marks = utils::tail(fm$trace[, "ess"], 1), sec_marks = tm,
  lambda0_cr = exp(fc$estimate[1]), mu0_cr = exp(fc$estimate[2]),
  loglik_cr = ll_c, ess_last_cr = utils::tail(fc$trace[, "ess"], 1), sec_cr = tc,
  lr = 2 * (ll_m - ll_c), build = as.character(utils::packageVersion("emphasis"))),
  file.path(out, sprintf("task-%03d.csv", task)), row.names = FALSE)

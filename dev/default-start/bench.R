# Default start and stopping rule: 40 cr trees (old default vs new), 10 cr trees from the
# box midpoint under the projected rule, 10 dd trees (old, new, 1000-iteration reference).
# Run on studio: Rscript bench.R <outdir> 9 (load_all path below); results in results.txt.
# Usage: Rscript bench.R <outdir> <workers>
args <- commandArgs(TRUE); out <- args[1]; W <- as.integer(args[2]); dir.create(out, FALSE, TRUE)
suppressMessages(devtools::load_all("~/pancho/emphasis-conv", quiet = TRUE)); library(ape)
cond_cr <- function(Tc) function(pars) { l <- exp(pars$beta_0); m <- exp(pars$gamma_0); d <- l - m
  p <- ifelse(abs(d) < 1e-12, 1 / (1 + l * Tc), d / (l - m * exp(-d * Tc))); p^2 }
# 40 cr trees, as in the earlier comparison
cr <- list(); seed <- 0
while (length(cr) < 40) { seed <- seed + 1; set.seed(1000 + seed)
  s <- simulate_tree(pars = log(c(0.5, 0.2)), max_t = 13, model = "cr")
  if (identical(s$status, "done") && Ntip(s$tes) >= 80 && Ntip(s$tes) <= 150) cr[[length(cr) + 1]] <- s$tes }
dd <- list(); seed <- 0
while (length(dd) < 10) { seed <- seed + 1; set.seed(2000 + seed)
  s <- tryCatch(simulate_tree(pars = c(log(0.6), -0.01, log(0.05), 0), max_t = 10, model = "dd"), error = function(e) NULL)
  if (!is.null(s) && identical(s$status, "done") && Ntip(s$tes) >= 60 && Ntip(s$tes) <= 150) dd[[length(dd) + 1]] <- s$tes }
saveRDS(list(cr = cr, dd = dd), file.path(out, "trees.rds"))
ctl_cr <- list(lower_bound = c(-4, -7), upper_bound = c(1.5, 1.5), sample_size = 200, num_threads = 1, max_time = 1800)
mid_cr <- (ctl_cr$lower_bound + ctl_cr$upper_bound) / 2
tasks <- c(
  lapply(1:40, function(i) list(set = "cr", i = i, v = "old")),
  lapply(1:40, function(i) list(set = "cr", i = i, v = "new")),
  lapply(1:10, function(i) list(set = "cr", i = i, v = "rule_only")),
  lapply(1:10, function(i) list(set = "dd", i = i, v = "old")),
  lapply(1:10, function(i) list(set = "dd", i = i, v = "new")),
  lapply(1:10, function(i) list(set = "dd", i = i, v = "ref")))
run <- function(t) {
  f <- file.path(out, sprintf("%s_%02d_%s.rds", t$set, t$i, t$v)); if (file.exists(f)) return(NULL)
  set.seed(t$i); t0 <- proc.time()[[3]]
  if (t$set == "cr") {
    tr <- cr[[t$i]]; Tc <- max(branching.times(tr))
    ctl <- ctl_cr; init <- NULL
    if (t$v == "old") { ctl$stop_rule <- "rel_change"; init <- mid_cr }
    if (t$v == "rule_only") init <- mid_cr
    fit <- estimate_rates(tr, model = "cr", cond = cond_cr(Tc), init_pars = init, control = ctl)
  } else {
    tr <- dd[[t$i]]
    bx <- auto_bounds(tr, model = "dd", train_surv_gam = FALSE, verbose = FALSE)
    ctl <- list(lower_bound = bx$lower_bound, upper_bound = bx$upper_bound, sample_size = 200, num_threads = 1, max_time = 1800)
    init <- NULL
    if (t$v == "old") { ctl$stop_rule <- "rel_change"; init <- (bx$lower_bound + bx$upper_bound) / 2 }
    if (t$v == "ref") { ctl$stop_rule <- "rel_change"; ctl$max_iter <- 1000; ctl$tol <- 1e-4; ctl$max_time <- 3600 }
    fit <- suppressWarnings(estimate_rates(tr, model = "dd", init_pars = init, control = ctl))
  }
  m <- fit$details$mcem
  r <- list(task = t, pars = fit$pars, loglik = fit$loglik, iter = fit$details$iterations,
            stop = fit$details$stop_reason, secs = proc.time()[[3]] - t0, trace = m)
  saveRDS(r, f); NULL
}
invisible(parallel::mclapply(tasks, function(t) tryCatch(run(t), error = function(e) {
  saveRDS(list(task = t, error = conditionMessage(e)), file.path(out, sprintf("%s_%02d_%s.err.rds", t$set, t$i, t$v))) }),
  mc.cores = W, mc.preschedule = FALSE))
cat("done\n")

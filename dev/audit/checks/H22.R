## H22: E-step recovery in .mcem_dynamic_fresh (R/emphasis.R:64-89)
##   fail    -> maxN <- min(2 maxN, 50000); pars <- 0.8 pars + 0.2 (lb+ub)/2; next
##   success -> maxN <- orig_maxN
## Claims to test
##  A. maxN >= sample_size is documented (inference.R:60-62) but not enforced; the
##     default maxN = 2000 (inference.R:114) is a fixed value, so the
##     max(2000, 10*N) fallback (inference.R:365) never fires for a default control.
##  B. A structural failure cause (sample_size > maxN) makes the loop alternate
##     fail/success; every failure pulls theta 20 % toward the box centre and the
##     M-step starts from the pulled value, so the fixed point of the iteration is
##     displaced toward the centre relative to a run with no failures.
##  C. After 8 consecutive failures the returned `pars` is the perturbed vector
##     0.8^8 theta0 + (1 - 0.8^8) centre, and estimate_rates reports it as `pars`.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({library(emphasis); library(ape)})
options(width = 150)

set.seed(22)
tr20 <- ape::rphylo(20, birth = 1, death = 0.3)
brts20 <- sort(ape::branching.times(tr20), decreasing = TRUE)
cat(sprintf("20-tip tree: crown age %.3f, %d branching times\n", brts20[1], length(brts20)))

lb <- c(0, 0); ub <- c(4, 4); centre <- (lb + ub) / 2
init <- c(1, 0.3)

## ---------------------------------------------------------------------------
cat("\n=== A. enforcement of maxN >= sample_size, and the default maxN ===\n")
ctrl_def <- emphasis:::estimate_rates_control("mcem")
cat(sprintf("default control: sampling=%s  maxN=%s  sample_size=%s\n",
            ctrl_def$sampling, ctrl_def$maxN, ctrl_def$sample_size))
ctrl_big <- emphasis:::.resolve_control_aliases(
  modifyList(ctrl_def, list(num_trees = 3000L)), "mcem", list(num_trees = 3000L))
cat(sprintf("after user num_trees=3000: sample_size=%d  maxN=%d  (maxN < sample_size: %s)\n",
            ctrl_big$sample_size, ctrl_big$maxN, ctrl_big$maxN < ctrl_big$sample_size))
cat(sprintf("is.null(maxN) fallback would give %d, but maxN is not NULL so it is skipped\n",
            max(2000L, 10L * ctrl_big$sample_size)))

# Does estimate_rates refuse sample_size > maxN?  (thinning driver, 2 iterations)
res_A <- tryCatch(
  suppressWarnings(estimate_rates(tr20, model = "cr", method = "mcem", init_pars = init,
                 control = list(sampling = "dynamic_fresh", num_trees = 300L, maxN = 200L,
                                max_iter = 2L, lower_bound = lb, upper_bound = ub))),
  error = function(e) e)
if (inherits(res_A, "error")) {
  cat("estimate_rates ERRORED on maxN < sample_size:", conditionMessage(res_A), "\n")
} else {
  cat(sprintf("estimate_rates ACCEPTED maxN=200 < sample_size=300: %d mcem rows in 2 iterations, stop_reason=%s\n",
              res_A$details$iterations, res_A$details$stop_reason))
}

## ---------------------------------------------------------------------------
cat("\n=== B. structural alternation: sample_size=300 with maxN=200 vs maxN=6000 ===\n")
cat("   thinning driver, cr, lb=(0,0) ub=(4,4) centre=(2,2), init=(1,0.3), num_threads=1\n")

run_trace <- function(maxN, max_iter, label) {
  msgs <- character(0)
  fit <- withCallingHandlers(
    suppressWarnings(estimate_rates(tr20, model = "cr", method = "mcem", init_pars = init,
                   control = list(sampling = "dynamic_fresh", num_trees = 300L, maxN = maxN,
                                  max_iter = max_iter, tol = 1e-3, patience = 3L,
                                  lower_bound = lb, upper_bound = ub,
                                  verbose = FALSE))),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  d <- fit$details
  n_fail <- if (d$stop_reason == "max_iter") max_iter - d$iterations else NA_integer_  # verbose=FALSE: no messages
  list(label = label, pars = fit$pars, loglik = fit$loglik, n_success = d$iterations,
       n_fail = n_fail, stop = d$stop_reason, mcem = d$mcem,
       rej = if (!is.null(d$mcem)) d$mcem$rejected else NA)
}

# The failure messages are only emitted with verbose=TRUE; re-run one alternating
# fit with verbose to show the fail/success pattern explicitly.
cat("\n-- B.0 one alternating run with verbose messages (first 10 iterations) --\n")
msgs <- character(0)
fit0 <- withCallingHandlers(
  suppressWarnings(estimate_rates(tr20, model = "cr", method = "mcem", init_pars = init,
                 control = list(sampling = "dynamic_fresh", num_trees = 300L, maxN = 200L,
                                max_iter = 10L, lower_bound = lb, upper_bound = ub,
                                verbose = TRUE))),
  message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
cat(paste0("  ", grep("Iteration", msgs, value = TRUE)), sep = "\n")
cat(sprintf("  -> %d successful rows, stop_reason=%s\n", fit0$details$iterations, fit0$details$stop_reason))

cat("\n-- B.1 replicated endpoints (3 x each) --\n")
tB <- proc.time()[3]
set.seed(2222)
reps <- 3L
clean <- lapply(seq_len(reps), function(r) run_trace(6000L, 30L, "clean maxN=6000"))
alt   <- lapply(seq_len(reps), function(r) run_trace(200L,  60L, "alt   maxN=200"))
summ <- function(lst) {
  P <- t(sapply(lst, `[[`, "pars"))
  data.frame(label = lst[[1]]$label,
             lambda = P[, 1], mu = P[, 2],
             loglik = sapply(lst, `[[`, "loglik"),
             n_success = sapply(lst, `[[`, "n_success"),
             n_fail = sapply(lst, `[[`, "n_fail"),
             stop = sapply(lst, `[[`, "stop"),
             max_rej = sapply(lst, function(x) max(x$rej)))
}
tab <- rbind(summ(clean), summ(alt))
print(tab, row.names = FALSE)
cat(sprintf("(B.1 elapsed %.0fs)\n", proc.time()[3] - tB))
mc <- colMeans(as.matrix(tab[tab$label == "clean maxN=6000", c("lambda", "mu")]))
ma <- colMeans(as.matrix(tab[tab$label == "alt   maxN=200",  c("lambda", "mu")]))
cat(sprintf("\nmean clean: lambda=%.4f mu=%.4f | mean alt: lambda=%.4f mu=%.4f | centre=(%.1f,%.1f)\n",
            mc[1], mc[2], ma[1], ma[2], centre[1], centre[2]))
cat(sprintf("displacement alt-clean: dlambda=%+.4f  dmu=%+.4f  (toward centre: %s, %s)\n",
            ma[1] - mc[1], ma[2] - mc[2],
            sign(ma[1] - mc[1]) == sign(centre[1] - mc[1]),
            sign(ma[2] - mc[2]) == sign(centre[2] - mc[2])))
# Implied per-success EM contraction r from the alternating fixed point:
#   theta* = [(1-r) thetahat + 0.2 r c] / (1 - 0.8 r)  =>  solve for r per coordinate
r_imp <- sapply(1:2, function(j) {
  f <- function(r) ((1 - r) * mc[j] + 0.2 * r * centre[j]) / (1 - 0.8 * r) - ma[j]
  tryCatch(uniroot(f, c(1e-6, 0.999999))$root, error = function(e) NA_real_)
})
cat(sprintf("implied EM contraction r (lambda, mu) = %.3f, %.3f\n", r_imp[1], r_imp[2]))

cat("\n-- B.2 last 6 successful rows of one clean and one alternating run --\n")
cat("clean:\n"); print(tail(clean[[1]]$mcem[, c("par1", "par5", "fhat", "delta_max", "rejected")], 6), row.names = FALSE)
cat("alt:\n");   print(tail(alt[[1]]$mcem[,   c("par1", "par5", "fhat", "delta_max", "rejected")], 6), row.names = FALSE)

# Reference: unconditioned CR MLE from DDD (same likelihood target as cond=NULL)
ref <- tryCatch({
  o <- utils::capture.output(m <- DDD::bd_ML(brts = brts20, cond = 0, btorph = 0, soc = 2,
                                             initparsopt = c(1, 0.3), idparsopt = 1:2,
                                             parsfix = c(0, 0), idparsfix = 3:4, verbose = FALSE))
  m
}, error = function(e) NULL)
if (!is.null(ref)) cat(sprintf("DDD::bd_ML (cond=0): lambda0=%.4f mu0=%.4f loglik=%.4f\n",
                               ref$lambda0, ref$mu0, ref$loglik))

## ---------------------------------------------------------------------------
cat("\n=== C. eight consecutive failures: returned pars is the perturbed vector ===\n")
set.seed(3)
tr6 <- ape::rphylo(6, birth = 1, death = 0.3)
tC0 <- proc.time()[3]
msgsC <- character(0); warnC <- character(0)
fitC <- withCallingHandlers(
  estimate_rates(tr6, model = "cr", method = "mcem", init_pars = init,
                 control = list(sampling = "dynamic_fresh", num_trees = 50001L, maxN = 50000L,
                                max_iter = 20L, lower_bound = lb, upper_bound = ub,
                                verbose = TRUE)),
  message = function(m) { msgsC <<- c(msgsC, conditionMessage(m)); invokeRestart("muffleMessage") },
  warning = function(w) { warnC <<- c(warnC, conditionMessage(w)); invokeRestart("muffleWarning") })
cat(sprintf("elapsed %.1fs; %d fail messages; stop_reason=%s; iterations=%s\n",
            proc.time()[3] - tC0, sum(grepl("E-step failed", msgsC)),
            fitC$details$stop_reason, fitC$details$iterations))
expected <- 0.8^8 * init + (1 - 0.8^8) * centre
cat(sprintf("returned pars : lambda=%.6f mu=%.6f\n", fitC$pars[1], fitC$pars[2]))
cat(sprintf("0.8^8*init+(1-0.8^8)*centre: lambda=%.6f mu=%.6f  (max abs diff %.2e)\n",
            expected[1], expected[2], max(abs(fitC$pars - expected))))
cat(sprintf("loglik=%s  AIC=%s  loglik_var=%s\n", fitC$loglik, fitC$AIC, fitC$loglik_var))
cat("warning text:", substr(paste(warnC, collapse = " | "), 1, 300), "\n")
cat("last fail message:", tail(grep("E-step failed", msgsC, value = TRUE), 1), "\n")

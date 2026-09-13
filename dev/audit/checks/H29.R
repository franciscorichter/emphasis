## H29 — does control$rho in emphasis_pipeline() reach the gam/cem/mcem stages?
##
## Two independent probes:
##  (1) trace(): record the `rho` argument at every back-end entry point the
##      pipeline can call (eval_logf, augment_trees, em_cpp, emphasis_cem,
##      .mcem_bdi, .mcem_dynamic_fresh, estimate_likelihood_surface,
##      auto_bounds).
##  (2) numeric signature: under CR the BDI proposal has constant log-weights,
##      so .augment_tree_bdi(brts, pars, rho) gives a deterministic fhat at
##      fixed pars and rho only adds n_tips*log(rho) (H2).  Re-evaluating the
##      pipeline's own mcem pars at rho = 1 and rho = 0.5 shows which rho the
##      reported loglik was computed with.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages(library(emphasis))

## 15-tip CR-like tree, fixed branching times (crown age 10)
brts <- c(10, 8.7, 7.9, 6.5, 6.1, 5.2, 4.4, 3.9, 3.1, 2.6, 2.0, 1.5, 0.9, 0.4)
n_tips <- length(brts) + 1L
lb <- c(0.05, 0.00)
ub <- c(1.50, 0.80)

## ---- probe 1: trace the rho argument -----------------------------------
h29log <- list()
mk_tracer <- function(fn) substitute({
  assign("h29log",
         c(get("h29log", envir = .GlobalEnv),
           list(list(fn = FN, rho = if (exists("rho", inherits = FALSE)) rho else NA))),
         envir = .GlobalEnv)
}, list(FN = fn))
targets <- c("eval_logf", "augment_trees", "em_cpp", "emphasis_cem",
             ".mcem_bdi", ".mcem_dynamic_fresh", "estimate_likelihood_surface",
             "auto_bounds")
ns <- asNamespace("emphasis")
for (fn in targets)
  suppressMessages(trace(fn, tracer = mk_tracer(fn), where = ns, print = FALSE))

tab <- function(tag) {
  df <- do.call(rbind, lapply(h29log, function(x) data.frame(fn = x$fn, rho = x$rho)))
  if (is.null(df)) { cat(tag, ": no calls traced\n"); return(invisible()) }
  cat("\n[", tag, "] rho argument seen at each back-end entry point:\n", sep = "")
  print(as.data.frame(table(fn = df$fn, rho = df$rho)))
  invisible(df)
}

small_ctrl <- list(
  num_threads = 1L, max_time = 120,
  gam  = list(n_grid = 30L, grid_points = 6L, sample_size = 10L),
  cem  = list(max_iter = 3L, num_particles = 10L, num_trees = 2L),
  mcem = list(max_iter = 3L, sample_size = 30L)
)

## Run A: rho at the TOP level of control (what README:168 tells the user to do)
h29log <<- list()
ctrl_A <- c(list(rho = 0.5, lower_bound = lb, upper_bound = ub), small_ctrl)
set.seed(1)
fitA <- emphasis_pipeline(brts, model = "cr", link = "linear",
                          stages = c("gam", "cem", "mcem"),
                          control = ctrl_A, verbose = FALSE)
dfA <- tab("Run A: control = list(rho = 0.5, ...) top level, stages gam/cem/mcem")

## Run B: rho nested per stage (the only way it can reach the back-ends)
h29log <<- list()
ctrl_B <- ctrl_A
ctrl_B$rho <- NULL
ctrl_B$gam$rho <- 0.5; ctrl_B$cem$rho <- 0.5; ctrl_B$mcem$rho <- 0.5
set.seed(1)
fitB <- emphasis_pipeline(brts, model = "cr", link = "linear",
                          stages = c("gam", "cem", "mcem"),
                          control = ctrl_B, verbose = FALSE)
dfB <- tab("Run B: rho = 0.5 nested in control$gam/$cem/$mcem")

## Run C: the bounds stage alone -- does auto_bounds get rho = 0.5?
h29log <<- list()
fitC <- tryCatch({
  setTimeLimit(elapsed = 150, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf), add = TRUE)
  emphasis_pipeline(brts, model = "cr", link = "linear",
                    stages = "bounds",
                    control = list(rho = 0.5, num_threads = 1L,
                                   bounds = list(n_test = 1L)),
                    verbose = TRUE)
}, error = function(e) { cat("Run C failed/timeout:", conditionMessage(e), "\n"); NULL })
dfC <- tab("Run C: stages = 'bounds' only, control$rho = 0.5")

for (fn in targets) suppressMessages(untrace(fn, where = ns))

## ---- probe 2: numeric signature on the mcem fit of Run A -----------------
cat("\n[probe 2] CR-BDI fhat is deterministic at fixed pars; rho enters as n_tips*log(rho)\n")
ev <- function(pars_compact, rho) {
  p8 <- emphasis:::.expand_pars(pars_compact, c(0L, 0L, 0L))
  emphasis:::.augment_tree_bdi(brts, p8, model_bin = c(0L, 0L, 0L),
                               sample_size = 20L, link = 0L, rho = rho)
}
pA <- fitA$fits$mcem$pars
rA <- list(rho1 = ev(pA, 1.0), rho05 = ev(pA, 0.5))
cat(sprintf("  pipeline mcem pars (Run A): beta_0=%.5f gamma_0=%.5f\n", pA[1], pA[2]))
cat(sprintf("  sd(lw) at rho=1: %.2e   sd(lw) at rho=0.5: %.2e  (zero-variance check)\n",
            sd(rA$rho1$weights), sd(rA$rho05$weights)))
cat(sprintf("  fhat(pars_A, rho=1)   = %.6f\n", rA$rho1$fhat))
cat(sprintf("  fhat(pars_A, rho=0.5) = %.6f\n", rA$rho05$fhat))
cat(sprintf("  difference            = %.6f   (n_tips*log(0.5) = %.6f)\n",
            rA$rho05$fhat - rA$rho1$fhat, n_tips * log(0.5)))
cat(sprintf("  pipeline-reported fits$mcem$loglik (Run A, rho=0.5 requested) = %.6f\n",
            fitA$fits$mcem$loglik))
cat(sprintf("    |loglik - fhat(rho=1)|   = %.2e\n", abs(fitA$fits$mcem$loglik - rA$rho1$fhat)))
cat(sprintf("    |loglik - fhat(rho=0.5)| = %.2e\n", abs(fitA$fits$mcem$loglik - rA$rho05$fhat)))

pB <- fitB$fits$mcem$pars
rB <- list(rho1 = ev(pB, 1.0), rho05 = ev(pB, 0.5))
cat(sprintf("\n  Run B (nested rho): fits$mcem$loglik = %.6f\n", fitB$fits$mcem$loglik))
cat(sprintf("    |loglik - fhat(rho=1)|   = %.2e\n", abs(fitB$fits$mcem$loglik - rB$rho1$fhat)))
cat(sprintf("    |loglik - fhat(rho=0.5)| = %.2e\n", abs(fitB$fits$mcem$loglik - rB$rho05$fhat)))

## Does the fit object record rho anywhere?
cat("\n[details] names(fits$mcem):", paste(names(fitA$fits$mcem), collapse = ", "), "\n")
cat("[details] names(fits$mcem$details):", paste(names(fitA$fits$mcem$details), collapse = ", "), "\n")
cat("[details] 'rho' appears anywhere in fits$mcem? ",
    any(grepl("rho", capture.output(str(fitA$fits$mcem, max.level = 3)), fixed = TRUE)), "\n")

## Also: estimate_rates() direct call honours control$rho (for contrast)
h29log <<- list()
suppressMessages(trace("eval_logf", tracer = mk_tracer("eval_logf"), where = ns, print = FALSE))
fd <- estimate_rates(brts, method = "mcem", model = "cr",
                     control = list(rho = 0.5, lower_bound = lb, upper_bound = ub,
                                    max_iter = 2L, sample_size = 20L, num_threads = 1L))
suppressMessages(untrace("eval_logf", where = ns))
tab("Direct estimate_rates(method='mcem', control=list(rho=0.5,...))")

## Verdict
top_level_reached <- !is.null(dfA) && any(dfA$rho == 0.5 & dfA$fn != "auto_bounds", na.rm = TRUE)
cat("\nVERDICT: top-level control$rho reached gam/cem/mcem back-ends in Run A? ",
    top_level_reached, "\n")
cat("         nested control$<stage>$rho reached back-ends in Run B? ",
    !is.null(dfB) && any(dfB$rho == 0.5, na.rm = TRUE), "\n")

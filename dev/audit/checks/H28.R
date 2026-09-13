# H28: auto_bounds infeasible-centre fallback is logged as "ok" by the pipeline
# and later stages run unconditioned on the wide box without a log flag.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({library(emphasis); library(ape)})
set.seed(1)

cat("=== Part A: forced infeasibility via a real tree ===\n")
# 40 tips, crown age 1e-3: tip_lo = floor(0.1*40) = 4, but no candidate rate in
# .find_feasible_center (lam <= 5*r_hat/(1-0.5) = 300, r_hat = log(20)/0.1)
# can grow 2 -> >=4 lineages in 1e-3 time units in >=50% of 5 sims.
tr <- rlineage(1, 0, Tmax = 5); tr <- drop.fossil(tr)
while (Ntip(tr) < 30 || Ntip(tr) > 60) { tr <- drop.fossil(rlineage(1, 0, Tmax = 5)) }
tr$edge.length <- tr$edge.length / max(branching.times(tr)) * 1e-3
cat("n_tips =", Ntip(tr), " crown age =", max(branching.times(tr)), "\n")

warn_A <- character()
resA <- withCallingHandlers(
  emphasis_pipeline(tr, model = "cr", link = "linear",
                    stages = c("bounds", "gam"),
                    control = list(num_threads = 1L, max_time = 120,
                                   gam = list(n_grid = 12, sample_size = 5)),
                    verbose = FALSE),
  warning = function(w) { warn_A <<- c(warn_A, conditionMessage(w)); invokeRestart("muffleWarning") })

cat("warnings raised:\n"); print(warn_A)
cat("stage log:\n"); print(resA$log[, c("stage", "status", "loglik")])
cat("result$cond          =", resA$cond, "\n")
cat("bounds$survival_gam  =", class(resA$bounds$survival_gam), "\n")
cat("bounds$center        =", if (is.null(resA$bounds$center)) "NULL" else "non-NULL", "\n")
wide <- emphasis:::.wide_bounds(emphasis:::.resolve_model("cr"), 0L,
                                max(branching.times(tr)), Ntip(tr))
cat("lower_bound == wide$lb :", isTRUE(all.equal(unname(resA$bounds$lower_bound), wide$lb)), "\n")
cat("upper_bound == wide$ub :", isTRUE(all.equal(unname(resA$bounds$upper_bound), wide$ub)), "\n")
cat("fits$gam$cond        =", if (!is.null(resA$fits$gam)) resA$fits$gam$cond else NA, "\n")
cat("log columns:", paste(names(resA$log), collapse = ", "), "\n")
cat("any log column mentions cond/downgrade/infeasible? ",
    any(grepl("cond|downgrad|infeas", unlist(resA$log), ignore.case = TRUE)), "\n")
cat("printed summary:\n"); print(resA)

cat("\n=== Part B: same pipeline on a feasible tree (control) ===\n")
tr2 <- drop.fossil(rlineage(0.4, 0.1, Tmax = 8))
while (Ntip(tr2) < 15 || Ntip(tr2) > 40) tr2 <- drop.fossil(rlineage(0.4, 0.1, Tmax = 8))
cat("n_tips =", Ntip(tr2), " crown age =", round(max(branching.times(tr2)), 3), "\n")
warn_B <- character()
resB <- withCallingHandlers(
  emphasis_pipeline(tr2, model = "cr", link = "linear",
                    stages = c("bounds", "gam"),
                    control = list(num_threads = 1L, max_time = 120,
                                   gam = list(n_grid = 12, sample_size = 5)),
                    verbose = FALSE),
  warning = function(w) { warn_B <<- c(warn_B, conditionMessage(w)); invokeRestart("muffleWarning") })
cat("warnings raised:\n"); print(warn_B)
cat("stage log:\n"); print(resB$log[, c("stage", "status", "loglik")])
cat("result$cond          =", resB$cond, "\n")
cat("bounds$survival_gam  =", class(resB$bounds$survival_gam)[1], "\n")
cat("fits$gam$cond        =", if (!is.null(resB$fits$gam)) resB$fits$gam$cond else NA, "\n")

cat("\n=== Part C: bounds-row of the log is identical between A (infeasible) and B (feasible)? ===\n")
rowA <- resA$log[resA$log$stage == "bounds", c("status", "loglik", "AIC", "pars")]
rowB <- resB$log[resB$log$stage == "bounds", c("status", "loglik", "AIC", "pars")]
print(rowA); print(rowB)
cat("identical bounds-row (status/loglik/AIC/pars):", identical(rowA, rowB), "\n")

cat("\n=== Part D: does the hypothesis' own example (3-tip tree) force infeasibility? ===\n")
tr3 <- read.tree(text = "((a:1,b:1):1,c:2);")
warn_D <- character()
abD <- withCallingHandlers(
  auto_bounds(tr3, model = "cr", link = "linear", num_threads = 1L,
              verbose = FALSE, train_surv_gam = FALSE),
  warning = function(w) { warn_D <<- c(warn_D, conditionMessage(w)); invokeRestart("muffleWarning") })
cat("3-tip warnings:", if (length(warn_D)) warn_D else "<none>", "\n")
cat("3-tip center:", if (is.null(abD$center)) "NULL (infeasible)" else paste(round(abD$center, 4), collapse = ", "), "\n")

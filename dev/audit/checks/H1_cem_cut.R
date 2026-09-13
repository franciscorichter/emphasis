.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({library(emphasis); library(ape)})
lam <- 0.2; mu <- 0.05
set.seed(7); repeat { t <- ape::rlineage(lam, mu, Tmax = 25); t2 <- ape::drop.fossil(t); if (Ntip(t2) >= 35 && Ntip(t2) <= 60) break }
b <- sort(as.numeric(branching.times(t2)), decreasing = TRUE); s <- 1e8
lb <- c(0.001, 0) / s; ub <- c(1, 0.5) / s
ws <- character(0); t0 <- Sys.time()
f <- withCallingHandlers(tryCatch(estimate_rates(b * s, method = "cem", model = "cr", init_pars = c(lam, mu) / s,
        control = list(lower_bound = lb, upper_bound = ub, max_time = 60, num_particles = 4L, max_iter = 1L, maxN = 10L, num_trees = 1L, verbose = FALSE)),
     error = function(e) { cat("cem ERROR:", conditionMessage(e), "\n"); NULL }), warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") })
cat(sprintf("cem at s=1e8 (%.0fs): pars=%s loglik=%s stop=%s warnings: %s\n", as.numeric(Sys.time() - t0, units = "secs"),
    if (is.null(f)) "NULL" else paste(signif(f$pars, 3), collapse = ","), if (is.null(f)) "NULL" else format(f$loglik),
    if (is.null(f)) "" else paste(f$stop_reason, f$cem$stop_reason), paste(substr(ws, 1, 120), collapse = " | ")))
if (!is.null(f)) str(f[setdiff(names(f), c("cem","tree","brts"))], max.level = 1)

## H29 replication with varied inputs: different tree (25 tips from ape::rlineage),
## model "dd", link "exponential", MCEM with thinning sampler (dynamic_fresh),
## each stage run ALONE, top-level rho vs nested rho. Trace rho at back-ends.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({library(emphasis); library(ape)})

set.seed(7)
tr <- NULL
repeat { t0 <- rlineage(0.5, 0.1, Tmax = 6); t0 <- drop.fossil(t0)
         if (!is.null(t0) && Ntip(t0) >= 20 && Ntip(t0) <= 40) { tr <- t0; break } }
brts <- sort(as.numeric(branching.times(tr)), decreasing = TRUE)
cat("n_tips =", Ntip(tr), " crown =", round(max(brts), 3), "\n")
lb <- c(0.05, 0.0, -0.05, -0.05); ub <- c(1.5, 0.6, 0.05, 0.05)   # dd compact: beta0,gamma0,beta_N,gamma_N

h29log <- list()
mk_tracer <- function(fn) substitute({
  assign("h29log", c(get("h29log", envir = .GlobalEnv),
    list(list(fn = FN, rho = if (exists("rho", inherits = FALSE)) rho else NA))), envir = .GlobalEnv)
}, list(FN = fn))
targets <- c("eval_logf", "augment_trees", "em_cpp", "emphasis_cem",
             ".mcem_bdi", ".mcem_dynamic_fresh", "estimate_likelihood_surface", "auto_bounds")
ns <- asNamespace("emphasis")
for (fn in targets) suppressMessages(trace(fn, tracer = mk_tracer(fn), where = ns, print = FALSE))
tab <- function(tag) {
  df <- do.call(rbind, lapply(h29log, function(x) data.frame(fn = x$fn, rho = x$rho)))
  cat("\n[", tag, "]\n", sep = "")
  if (is.null(df)) { cat("  no calls traced\n"); return(invisible(NULL)) }
  print(as.data.frame(table(fn = df$fn, rho = df$rho)))
  invisible(df)
}
base <- list(num_threads = 1L, max_time = 60, lower_bound = lb, upper_bound = ub,
             gam  = list(n_grid = 20L, grid_points = 5L, sample_size = 5L),
             cem  = list(max_iter = 2L, num_particles = 8L, num_trees = 2L),
             mcem = list(max_iter = 2L, sample_size = 10L, sampling = "dynamic_fresh", maxN = 200L))
res <- list()
for (st in c("gam", "cem", "mcem")) {
  h29log <<- list()
  ctl <- c(list(rho = 0.6), base)
  f <- tryCatch(emphasis_pipeline(brts, model = "dd", link = "exponential", stages = st,
                                  control = ctl, verbose = FALSE),
                error = function(e) { cat("ERR", st, conditionMessage(e), "\n"); NULL })
  df <- tab(paste("TOP-LEVEL rho=0.6, stage =", st))
  res[[paste0("top_", st)]] <- if (!is.null(df)) unique(df$rho) else NA
  # is rho recorded anywhere in the fit?
  if (!is.null(f)) cat("  'rho' in str(fit)?", any(grepl("rho", capture.output(str(f, max.level = 4)), fixed = TRUE)), "\n")

  h29log <<- list()
  ctl2 <- base; ctl2[[st]]$rho <- 0.6
  f2 <- tryCatch(emphasis_pipeline(brts, model = "dd", link = "exponential", stages = st,
                                   control = ctl2, verbose = FALSE),
                 error = function(e) { cat("ERR", st, conditionMessage(e), "\n"); NULL })
  df2 <- tab(paste("NESTED rho=0.6 in control$", st, ", stage =", st))
  res[[paste0("nested_", st)]] <- if (!is.null(df2)) unique(df2$rho) else NA
}
for (fn in targets) suppressMessages(untrace(fn, where = ns))
cat("\nSUMMARY (unique rho values seen at back-ends):\n"); print(res)

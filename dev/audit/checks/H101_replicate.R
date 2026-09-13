## H101 replication — independent check of "16 functions have no tests".
## Part 1: DYNAMIC coverage. The verifier only counted name *mentions*. A test
##         can exercise a function without naming it (estimate_rates() calls
##         .resolve_control_aliases before the bounds check that test-inference.R
##         expects to error).  Trace every one of the 16 names in the installed
##         namespace and run the existing test suite against the installed
##         package (no compilation, no load_all).
## Part 2: smoke tests varied from the verifier's (20-tip tree, dd model,
##         exponential link where the verifier used cr/linear; thinning + bdi).
## Part 3: the incidental "conditioned" label on a bounds-only pipeline.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({ library(emphasis); library(ape); library(testthat) })
pkg <- "/Users/pancho/Code/emphasis"

fns <- c("emphasis_pipeline", "auto_bounds", "print.emphasis_pipeline",
         "diagnose_mcem", "diagnose_cem", ".run_mcem", ".run_cem", ".run_gam",
         ".resolve_control_aliases", ".validate_linear_init_pars",
         ".build_cond_fun", "estimate_likelihood_surface",
         ".sim_tree_conditional", ".aug_to_Ltable", "prune_to_extant",
         "select_diversification_model")

## ---------------------------------------------------------------- Part 1
cat("=== Part 1: dynamic coverage of the installed test suite\n")
if (!nzchar(Sys.getenv("H101_SKIP_PART1"))) {
hits <- new.env()
for (f in fns) assign(f, 0L, envir = hits)
ns <- asNamespace("emphasis")
for (f in fns) {
  tracer <- bquote(assign(.(f), get(.(f), envir = hits) + 1L, envir = hits))
  suppressMessages(trace(f, where = ns, tracer = tracer, print = FALSE))
}
res <- testthat::test_dir(file.path(pkg, "tests", "testthat"),
                          package = "emphasis", load_package = "installed",
                          reporter = "summary", stop_on_failure = FALSE)
df <- as.data.frame(res)
cat(sprintf("\nsuite: %d blocks, %d passed, %d failed, %d skipped, %d warnings\n",
            nrow(df), sum(df$passed), sum(df$failed), sum(df$skipped), sum(df$warning)))
for (f in fns) suppressMessages(untrace(f, where = ns))
cov <- data.frame(fn = fns, calls_during_tests = sapply(fns, get, envir = hits),
                  row.names = NULL)
print(cov)
cat("functions executed by at least one (unskipped) test:",
    sum(cov$calls_during_tests > 0), "of", length(fns), "\n")
}

## ---------------------------------------------------------------- Part 2
cat("\n=== Part 2: smoke tests, varied (20-tip tree, dd / exponential)\n")
ok <- character(); bad <- character()
run <- function(name, expr) {
  t0 <- proc.time()[3]
  out <- tryCatch(list(ok = TRUE, val = expr, err = NA_character_),
                  error = function(e) list(ok = FALSE, val = NULL, err = conditionMessage(e)))
  cat(sprintf("[%-34s] %s (%.1fs)%s\n", name, if (out$ok) "OK" else "ERROR",
              proc.time()[3] - t0, if (out$ok) "" else paste0("  -> ", out$err)))
  if (out$ok) ok <<- c(ok, name) else bad <<- c(bad, name)
  invisible(out$val)
}
set.seed(7)
tree20 <- rphylo(20, birth = 0.6, death = 0.2, T0 = 6, fossils = FALSE)
brts20 <- sort(branching.times(tree20), decreasing = TRUE)
cat("tree20: ntip =", Ntip(tree20), " crown age =", round(brts20[1], 3), "\n")

## prune_to_extant on a fresh simulate_tree with extinct lineages
sim <- NULL
for (k in 1:80) {
  s <- tryCatch(simulate_tree(pars = c(0.7, 0.35), max_t = 4, model = "cr"),
                error = function(e) NULL)
  if (!is.null(s) && !is.null(s$tas) && Ntip(s$tas) > Ntip(s$tes) && Ntip(s$tes) >= 5) { sim <- s; break }
}
if (!is.null(sim)) {
  pr <- run("prune_to_extant", emphasis:::prune_to_extant(sim$tas))
  cat("   tas", Ntip(sim$tas), "-> pruned", Ntip(pr), " tes", Ntip(sim$tes),
      " brts equal:", isTRUE(all.equal(sort(branching.times(pr)), sort(branching.times(sim$tes)), tolerance = 1e-8)), "\n")
} else cat("   (no simulate_tree draw with extinct lineages in 80 tries)\n")

## auto_bounds: dd model, EXPONENTIAL link (verifier used cr/linear)
ab <- run("auto_bounds dd/exponential", suppressMessages(
  auto_bounds(tree20, model = "dd", link = "exponential", n_test = 3,
              bisect_steps = 4, num_threads = 1)))
if (!is.null(ab)) {
  cat("   lb:", round(ab$lower_bound, 4), "\n   ub:", round(ab$upper_bound, 4),
      "\n   lb<ub:", all(ab$lower_bound < ab$upper_bound),
      " gam:", class(ab$survival_gam)[1], "\n")
}
## auto_bounds: cr / linear on the 20-tip tree as well (for the pipeline below)
abcr <- run("auto_bounds cr/linear", suppressMessages(
  auto_bounds(tree20, model = "cr", link = "linear", n_test = 3, bisect_steps = 4,
              num_threads = 1)))

## pipeline, bounds only, dd model
pp <- run("emphasis_pipeline dd bounds", suppressMessages(
  emphasis_pipeline(tree20, model = "dd", stages = "bounds",
                    control = list(num_threads = 1, bounds = list(n_test = 3, bisect_steps = 4)),
                    verbose = FALSE)))
if (!is.null(pp)) {
  cat("   class:", paste(class(pp), collapse = "/"), " best:", pp$best_stage,
      " log:", paste(pp$log$stage, pp$log$status), " cond flag:", pp$cond, "\n")
  run("print.emphasis_pipeline", capture.output(print(pp)))
}
run("pipeline cem w/o bounds errors",
    tryCatch({ emphasis_pipeline(tree20, stages = "cem", verbose = FALSE); stop("did NOT error") },
             error = function(e) if (grepl("bounds", conditionMessage(e))) TRUE else stop(e)))

## control alias resolution: gam method (verifier covered mcem/cem)
ra <- run(".resolve_control_aliases gam", {
  c1 <- emphasis:::.resolve_control_aliases(emphasis:::estimate_rates_control("gam"), "gam",
                                            user_ctrl = list(num_points = 12))
  c2 <- emphasis:::.resolve_control_aliases(emphasis:::estimate_rates_control("gam"), "gam",
                                            user_ctrl = list(sample_size = 7))
  list(c1 = c1, c2 = c2)
})
if (!is.null(ra)) cat("   c1 num_points/num_particles:", ra$c1$num_points, ra$c1$num_particles,
                      " c2 sample_size/num_trees:", ra$c2$sample_size, ra$c2$num_trees, "\n")

## .validate_linear_init_pars on a 20-tip tree, dd
vi <- run(".validate_linear_init_pars", emphasis:::.validate_linear_init_pars(
  c(0.5, -0.1, 0.1, 0), c(1L, 0L, 0L), brts20,
  lower_bound = c(0.01, -1, 0, -1), upper_bound = c(3, 1, 3, 1)))
if (!is.null(vi)) cat("   init:", round(vi, 4), " lambda at N=20:", vi[1] + vi[2] * 20, "\n")

## .build_cond_fun on the dd/exponential GAM
if (!is.null(ab) && !is.null(ab$survival_gam)) {
  cf <- run(".build_cond_fun", emphasis:::.build_cond_fun(ab$survival_gam, c(1L, 0L, 0L)))
  if (!is.null(cf)) cat("   logP at centre:", cf((ab$lower_bound + ab$upper_bound) / 2), "\n")
}

## .sim_tree_conditional: thinning and bdi with a dd model
for (smp in c("thinning", "bdi")) {
  st <- run(paste(".sim_tree_conditional", smp), emphasis:::.sim_tree_conditional(
    tree20, pars = c(0.6, -0.005, 0.15, 0), model_bin = c(1L, 0L, 0L), link = 0L,
    n_trees = 3L, method = smp, num_threads = 1L))
  if (!is.null(st)) cat("   n:", length(st$trees), " tips:", sapply(st$trees, Ntip),
                        " log_q:", round(st$log_q, 2), "\n")
}

## estimate_likelihood_surface: 2-point grid, dd
els <- run("estimate_likelihood_surface dd", suppressMessages(
  emphasis:::estimate_likelihood_surface(tree20,
    pars_mat = rbind(c(0.5, -0.005, 0.1, 0), c(0.7, -0.01, 0.2, 0)),
    model = "dd", sample_size = 30, num_threads = 1)))
if (!is.null(els)) print(els)

## .run_mcem via estimate_rates, thinning + bdi, dd; diagnose_mcem
for (smp in c("thinning", "bdi")) {
  fm <- run(paste(".run_mcem", smp), suppressMessages(estimate_rates(
    tree20, method = "mcem", model = "dd",
    control = list(lower_bound = c(0.05, -0.1, 0, -0.1), upper_bound = c(3, 0.1, 2, 0.1),
                   num_trees = 40, max_iter = 4, num_threads = 1, sampler = smp))))
  if (!is.null(fm)) {
    cat("   pars:", round(fm$pars, 3), " loglik:", round(fm$loglik, 2), "\n")
    dm <- run(paste("diagnose_mcem", smp), diagnose_mcem(fm, plot = FALSE))
    if (!is.null(dm)) cat("   diag names:", paste(names(dm), collapse = ","), "\n")
  }
}
## .run_cem + diagnose_cem, dd
fc <- run(".run_cem", suppressMessages(estimate_rates(
  tree20, method = "cem", model = "dd",
  control = list(lower_bound = c(0.05, -0.1, 0, -0.1), upper_bound = c(3, 0.1, 2, 0.1),
                 num_particles = 10, num_trees = 10, max_iter = 3, num_threads = 1))))
if (!is.null(fc)) {
  cat("   loglik:", round(fc$loglik, 2), "\n")
  dc <- run("diagnose_cem", diagnose_cem(fc, plot = FALSE))
  if (!is.null(dc)) cat("   diag names:", paste(names(dc), collapse = ","), "\n")
}
## .run_gam, dd
fg <- run(".run_gam", suppressMessages(estimate_rates(
  tree20, method = "gam", model = "dd",
  control = list(lower_bound = c(0.05, -0.1, 0, -0.1), upper_bound = c(3, 0.1, 2, 0.1),
                 n_grid = 30, sample_size = 20, num_threads = 1))))
if (!is.null(fg)) cat("   loglik:", round(fg$loglik, 2), "\n")

## .aug_to_Ltable via the mcem fit's augmented trees, if exposed
cat("\n=== Part 3: 'conditioned' label on bounds-only pipeline\n")
if (!is.null(pp)) {
  out <- capture.output(print(pp))
  cat("   header line:", out[1], "\n   pp$cond =", pp$cond, " best_stage =", pp$best_stage, "\n")
}

cat("\nsmoke OK:", length(ok), " ERROR:", length(bad), if (length(bad)) paste(":", bad), "\n")

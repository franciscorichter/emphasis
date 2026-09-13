## H101 — "emphasis_pipeline, auto_bounds, print.emphasis_pipeline, diagnose_mcem,
## diagnose_cem, .run_*, .resolve_control_aliases, .validate_linear_init_pars,
## .build_cond_fun, estimate_likelihood_surface, .sim_tree_conditional,
## .aug_to_Ltable, prune_to_extant, select_diversification_model have no tests."
##
## Part A: static census — does any file under tests/ mention each name?
## Part B: smoke tests on a 10-tip tree — do the untested functions run, and
##         what do they return?  (Decides whether "untested" also hides "broken".)
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({ library(emphasis); library(ape) })
pkg <- "/Users/pancho/Code/emphasis"
set.seed(101)

## ---------------------------------------------------------------- Part A
fns <- c("emphasis_pipeline", "auto_bounds", "print.emphasis_pipeline",
         "diagnose_mcem", "diagnose_cem", ".run_mcem", ".run_cem", ".run_gam",
         ".resolve_control_aliases", ".validate_linear_init_pars",
         ".build_cond_fun", "estimate_likelihood_surface",
         ".sim_tree_conditional", ".aug_to_Ltable", "prune_to_extant",
         "select_diversification_model")
test_files <- list.files(file.path(pkg, "tests"), pattern = "\\.R$",
                         recursive = TRUE, full.names = TRUE)
test_txt <- unlist(lapply(test_files, readLines, warn = FALSE))
r_files  <- list.files(file.path(pkg, "R"), pattern = "\\.R$", full.names = TRUE)
r_txt    <- unlist(lapply(r_files, readLines, warn = FALSE))
cat("=== Part A: static census over", length(test_files), "test files,",
    length(test_txt), "lines\n")
census <- data.frame(fn = fns,
  defined_in_R = sapply(fns, function(f) any(grepl(paste0("^\\s*`?", gsub(".", "\\.", f, fixed = TRUE), "`?\\s*<-\\s*function"), r_txt))),
  test_mentions = sapply(fns, function(f) sum(grepl(f, test_txt, fixed = TRUE))),
  exists_in_ns = sapply(fns, function(f) exists(f, envir = asNamespace("emphasis"))),
  row.names = NULL)
print(census)
cat("total test_that blocks:", sum(grepl("test_that(", test_txt, fixed = TRUE)),
    " skip() calls:", sum(grepl("^\\s*skip\\(", test_txt)), "\n")

## ---------------------------------------------------------------- Part B
cat("\n=== Part B: smoke tests on a 10-tip tree\n")
res <- list()
run <- function(name, expr) {
  t0 <- proc.time()[3]
  out <- tryCatch(list(ok = TRUE, val = expr, err = NA_character_),
                  error = function(e) list(ok = FALSE, val = NULL, err = conditionMessage(e)))
  out$secs <- round(proc.time()[3] - t0, 1)
  res[[name]] <<- out
  cat(sprintf("[%-32s] %s (%.1fs)%s\n", name, if (out$ok) "OK" else "ERROR",
              out$secs, if (out$ok) "" else paste0("  -> ", out$err)))
  invisible(out)
}

## 10-tip CR tree from ape (reproducible), plus a simulate_tree() object for tas
tree10 <- rphylo(10, birth = 0.5, death = 0.1, T0 = 5, fossils = FALSE)
brts10 <- sort(branching.times(tree10), decreasing = TRUE)
cat("tree10: ntip =", Ntip(tree10), " crown age =", round(brts10[1], 3), "\n")

## 1. prune_to_extant --------------------------------------------------------
sim <- NULL
for (k in 1:50) {
  s <- tryCatch(simulate_tree(pars = c(0.6, 0.3), max_t = 4, model = "cr"),
                error = function(e) NULL)
  if (!is.null(s) && !is.null(s$tas) && Ntip(s$tas) > Ntip(s$tes) &&
      Ntip(s$tes) >= 5 && Ntip(s$tes) <= 30) { sim <- s; break }
}
run("prune_to_extant", {
  stopifnot(!is.null(sim))
  pr <- emphasis:::prune_to_extant(sim$tas)
  h  <- node.depth.edgelength(pr)[seq_len(Ntip(pr))]
  c(ntip_tas = Ntip(sim$tas), ntip_tes = Ntip(sim$tes), ntip_pruned = Ntip(pr),
    ultrametric_spread = max(h) - min(h),
    brts_match = isTRUE(all.equal(sort(branching.times(pr)), sort(branching.times(sim$tes)),
                                  check.attributes = FALSE)))
})
print(res$prune_to_extant$val)

## 2. auto_bounds -------------------------------------------------------------
run("auto_bounds(cr,linear)", {
  ab <- auto_bounds(tree10, model = "cr", link = "linear", n_test = 3L,
                    bisect_steps = 4L, num_threads = 1L, verbose = FALSE)
  list(lb = ab$lower_bound, ub = ab$upper_bound, center = ab$center,
       has_gam = !is.null(ab$survival_gam), names = names(ab))
})
ab_cr <- res[["auto_bounds(cr,linear)"]]$val
if (!is.null(ab_cr)) { str(ab_cr[c("lb", "ub", "center", "has_gam")]) }

## 3. emphasis_pipeline stages = "bounds" + print method ---------------------
run("emphasis_pipeline(bounds only)", {
  pp <- emphasis_pipeline(tree10, model = "cr", link = "linear",
                          stages = "bounds",
                          control = list(num_threads = 1L, bounds = list(n_test = 3L)),
                          verbose = FALSE)
  list(class = class(pp), best_stage = pp$best_stage, pars = pp$pars,
       loglik = pp$loglik, log = pp$log[, c("stage", "status", "elapsed")],
       lb = pp$bounds$lower_bound, ub = pp$bounds$upper_bound)
})
pp_val <- res[["emphasis_pipeline(bounds only)"]]$val
if (!is.null(pp_val)) { cat("class:", pp_val$class, " best_stage:", pp_val$best_stage,
                            " loglik:", pp_val$loglik, "\n"); print(pp_val$log) }
run("print.emphasis_pipeline(bounds)", {
  pp <- emphasis_pipeline(tree10, model = "cr", stages = "bounds",
                          control = list(num_threads = 1L, bounds = list(n_test = 3L)),
                          verbose = FALSE)
  txt <- capture.output(print(pp))
  cat(txt, sep = "\n"); txt
})
## degenerate: bounds skipped without control bounds must stop()
run("emphasis_pipeline(no bounds->stop)", {
  e <- tryCatch(emphasis_pipeline(tree10, stages = "cem", verbose = FALSE),
                error = function(e) conditionMessage(e))
  stopifnot(is.character(e)); e
})

## 4. .resolve_control_aliases ------------------------------------------------
run(".resolve_control_aliases", {
  rca <- emphasis:::.resolve_control_aliases
  base <- emphasis:::estimate_rates_control("cem")
  a <- rca(utils::modifyList(base, list(num_points = 7L)), "cem", list(num_points = 7L))
  b <- rca(utils::modifyList(base, list(num_particles = 9L)), "cem", list(num_particles = 9L))
  d <- rca(utils::modifyList(base, list(num_particles = 9L, num_points = 7L)), "cem",
           list(num_particles = 9L, num_points = 7L))
  m <- rca(utils::modifyList(emphasis:::estimate_rates_control("mcem"), list(sample_size = 33L)),
           "mcem", list(sample_size = 33L))
  out <- c(old_only_propagates = a$num_particles == 7L && a$num_points == 7L,
           new_only_propagates = b$num_points == 9L && b$num_particles == 9L,
           both_new_wins       = d$num_points == 9L,
           mcem_sample_to_trees = m$num_trees == 33L)
  print(out); stopifnot(all(out)); out
})

## 5. .validate_linear_init_pars ---------------------------------------------
run(".validate_linear_init_pars", {
  v  <- emphasis:::.validate_linear_init_pars
  lb <- c(0, -1, 0, -1); ub <- c(3, 1, 1, 1)
  bad  <- v(c(0.5, -0.2, 0.1, -0.05), c(1L, 0L, 0L), brts10, lb, ub)   # lambda<0 at N=10
  good <- v(c(0.5, -0.01, 0.1, 0.0), c(1L, 0L, 0L), brts10, lb, ub)
  out <- c(clamped_slope = bad[2], lambda_at_N10 = bad[1] + bad[2] * 10,
           good_unchanged = isTRUE(all.equal(good, c(0.5, -0.01, 0.1, 0.0))))
  print(out)
  stopifnot(out["lambda_at_N10"] > 0, isTRUE(as.logical(out["good_unchanged"])))
  out
})

## 6. .build_cond_fun ---------------------------------------------------------
run(".build_cond_fun", {
  ab <- auto_bounds(tree10, model = "cr", link = "linear", n_test = 3L,
                    bisect_steps = 4L, num_threads = 1L, verbose = FALSE,
                    train_surv_gam = TRUE)
  stopifnot(!is.null(ab$survival_gam))
  cf <- emphasis:::.build_cond_fun(ab$survival_gam, c(0L, 0L, 0L))
  p8 <- emphasis:::.expand_pars(ab$center, c(0L, 0L, 0L))
  v  <- cf(p8)
  out <- c(logP_center = v, finite = is.finite(v), nonpositive = v <= 0)
  print(out); stopifnot(is.finite(v), v <= 0); out
})

## 7. .sim_tree_conditional / .aug_to_Ltable ----------------------------------
run(".sim_tree_conditional(thinning)", {
  r <- emphasis:::.sim_tree_conditional(tree10, pars = c(0.5, 0.2),
         model_bin = c(0L, 0L, 0L), n_trees = 3L, method = "thinning",
         num_threads = 1L)
  out <- list(n = length(r$trees),
              classes = sapply(r$trees, function(t) class(t)[1]),
              ntips = sapply(r$trees, function(t) if (inherits(t, "phylo")) Ntip(t) else NA),
              log_q = r$log_q)
  str(out); out
})
run(".sim_tree_conditional(bdi)", {
  r <- emphasis:::.sim_tree_conditional(tree10, pars = c(0.5, 0.2),
         model_bin = c(0L, 0L, 0L), n_trees = 3L, method = "bdi")
  out <- list(n = length(r$trees),
              ntips = sapply(r$trees, function(t) if (inherits(t, "phylo")) Ntip(t) else NA),
              log_q = r$log_q)
  str(out); out
})
run(".aug_to_Ltable(direct)", {
  aug <- emphasis:::.augment_tree_internal(tree10, pars = c(0.5, 0.4),
           model_bin = c(0L, 0L, 0L), sample_size = 20L, num_threads = 1L, link = 0L)
  L0  <- DDD::phylo2L(tree10)
  nm  <- sapply(aug$trees, function(df) sum(df$t_ext != 0 & df$parent_id != -1L & df$t_ext < 1e11))
  j   <- which(nm > 0)[1]; stopifnot(!is.na(j))
  L   <- emphasis:::.aug_to_Ltable(aug$trees[[j]], brts10[1], brts10, L0)
  n_missing <- nm[j]
  out <- c(rows_extant = nrow(L0), rows_aug = nrow(L), n_missing_in_df = n_missing,
           extinct_rows = sum(L[, 4] > 0), ok_L2phylo = inherits(
             tryCatch(DDD::L2phylo(L, dropextinct = FALSE), error = function(e) NULL), "phylo"))
  print(out); stopifnot(out["rows_aug"] == out["rows_extant"] + out["n_missing_in_df"]); out
})

## 8. estimate_likelihood_surface --------------------------------------------
run("estimate_likelihood_surface", {
  grid <- cbind(beta_0 = c(0.3, 0.5, 0.8), gamma_0 = c(0.05, 0.1, 0.2))
  s <- emphasis:::estimate_likelihood_surface(tree10, grid, model = "cr", sample_size = 50L,
                                   num_threads = 1L, verbose = FALSE)
  print(s); stopifnot(nrow(s) == 3L, all(is.finite(s$fhat))); s
})

## 9. .run_mcem / .run_cem / .run_gam via estimate_rates + diagnose_* -------
lb <- c(0.05, 0.0); ub <- c(2.0, 1.0)
run(".run_mcem(bdi)+diagnose_mcem", {
  fit <- estimate_rates(tree10, method = "mcem", model = "cr",
           control = list(lower_bound = lb, upper_bound = ub, num_trees = 50L,
                          max_iter = 4L, num_threads = 1L, verbose = FALSE))
  d <- diagnose_mcem(fit, plot = FALSE, lower_bound = lb, upper_bound = ub)
  out <- list(pars = fit$pars, loglik = fit$loglik, iters = nrow(fit$details$mcem),
              diag_names = names(d), diag_class = class(d))
  str(out); out
})
run(".run_mcem(thinning)+diagnose_mcem", {
  fit <- estimate_rates(tree10, method = "mcem", model = "cr",
           control = list(lower_bound = lb, upper_bound = ub, num_trees = 50L,
                          max_iter = 4L, num_threads = 1L, verbose = FALSE,
                          sampling = "dynamic_fresh"))
  d <- diagnose_mcem(fit, plot = FALSE, lower_bound = lb, upper_bound = ub)
  print(capture.output(print(d))[1:6])
  list(pars = fit$pars, loglik = fit$loglik, iters = nrow(fit$details$mcem))
})
run(".run_cem+diagnose_cem", {
  fit <- estimate_rates(tree10, method = "cem", model = "cr",
           control = list(lower_bound = lb, upper_bound = ub, num_particles = 10L,
                          num_trees = 10L, max_iter = 3L, num_threads = 1L,
                          verbose = FALSE))
  d <- diagnose_cem(fit, plot = FALSE, lower_bound = lb, upper_bound = ub)
  out <- list(pars = fit$pars, loglik = fit$loglik, diag_names = names(d),
              ess_frac = d$IS_quality$ESS_fraction)
  str(out); out
})
run(".run_gam", {
  fit <- estimate_rates(tree10, method = "gam", model = "cr",
           control = list(lower_bound = lb, upper_bound = ub, n_grid = 30L,
                          sample_size = 20L, num_threads = 1L, verbose = FALSE))
  out <- list(pars = fit$pars, loglik = fit$loglik, method = fit$method)
  str(out); out
})

## 10. select_diversification_model -------------------------------------------
run("select_diversification_model", {
  sel <- emphasis:::select_diversification_model(tree10,
           lower_bound = c(0.05, -0.2, 0, -0.2), upper_bound = c(2, 0.2, 1, 0.2),
           control = list(cem  = list(max_iter = 2L, num_particles = 8L, num_trees = 5L,
                                      num_threads = 1L),
                          mcem = list(num_trees = 30L, max_iter = 3L, num_threads = 1L)),
           verbose = FALSE)
  print(sel$summary)
  list(best = sel$best_model, models = names(sel$fits), class = class(sel))
})

## ---------------------------------------------------------------- summary
cat("\n=== SUMMARY\n")
tab <- data.frame(check = names(res),
                  ok = sapply(res, `[[`, "ok"),
                  secs = sapply(res, `[[`, "secs"),
                  err = sapply(res, function(r) substr(ifelse(is.na(r$err), "", r$err), 1, 90)),
                  row.names = NULL)
print(tab, right = FALSE)
cat("functions with zero test mentions:", sum(census$test_mentions == 0), "of", nrow(census), "\n")
cat("smoke failures:", sum(!tab$ok), "of", nrow(tab), "\n")

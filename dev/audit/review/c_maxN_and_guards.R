# (c) maxN semantics across the three paths, and (d) the 8-slot guards.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib-wave1", .libPaths()))
library(emphasis)
set.seed(3)
brts <- sort(ape::branching.times(ape::rcoal(15)), decreasing = TRUE)
lb <- c(0, 0); ub <- c(3, 3)

er <- function(...) tryCatch(estimate_rates(brts, method = "mcem", model = "cr",
                                            control = utils::modifyList(
                                              list(lower_bound = lb, upper_bound = ub,
                                                   num_trees = 30L, max_iter = 2L,
                                                   num_threads = 1L), list(...))),
                             error = function(e) paste("ERROR:", conditionMessage(e)))

cat("## (c1) maxN < num_trees, BDI sampler (BDI never uses maxN)\n")
r <- er(sampling = "bdi", maxN = 10L)
cat(if (is.character(r)) substr(r, 1, 160) else "no error", "\n\n")

cat("## (c2) maxN < num_trees, thinning sampler\n")
r <- er(sampling = "dynamic_fresh", maxN = 10L)
cat(if (is.character(r)) substr(r, 1, 120) else "no error", "\n\n")

cat("## (c3) does maxN change anything on the BDI path?\n")
set.seed(5); a <- er(sampling = "bdi", maxN = 2000L)
set.seed(5); b <- er(sampling = "bdi", maxN = 500000L)
cat("  identical pars:", isTRUE(all.equal(a$pars, b$pars)), "\n")
cat("  .mcem_bdi formals:", paste(names(formals(emphasis:::.mcem_bdi)), collapse = ","), "\n")
cat("  'maxN' among them:", "maxN" %in% names(formals(emphasis:::.mcem_bdi)), "\n\n")

cat("## (c4) thinning ratchet: details$maxN after a run with forced failures\n")
set.seed(9)
res <- tryCatch(emphasis:::.mcem_dynamic_fresh(
  brts = brts, pars = c(3, 0, 0, 0, 0.001, 0, 0, 0),
  sample_size = 20L, maxN = 20L, max_missing = 5L,
  lower_bound = c(0,0,0,0,0,0,0,0), upper_bound = c(5,0,0,0,5,0,0,0),
  max_iter = 4L, xtol = 1e-3, patience = 3L, num_threads = 1L,
  model = c(0L,0L,0L), link = 0L),
  error = function(e) paste("ERROR:", conditionMessage(e)))
if (is.list(res)) {
  cat("  start maxN = 20, end maxN =", res$maxN,
      " n_failed =", res$n_failed, " stop =", res$stop_reason, "\n")
} else cat(" ", substr(res, 1, 120), "\n")
cat("  BDI driver returns a maxN field:",
    "maxN" %in% names(tryCatch(er(sampling = "bdi")$details, error = function(e) list())), "\n")

cat("\n## (c5) documentation of maxN in the control help\n")
rd <- readLines("/Users/pancho/Code/emphasis/man/estimate_rates_control.Rd")
i <- grep("maxN", rd)
cat(paste(" ", rd[sort(unique(c(i, i+1, i+2, i+3)))], collapse = "\n"), "\n")
cat("\n sampling item:\n")
j <- grep("sampling", rd)
cat(paste(" ", rd[sort(unique(c(j, j+1, j+2)))], collapse = "\n"), "\n")

cat("\n## (d) 8-slot guards\n")
p8 <- c(0.5, 0, 0, 0, 0.1, 0, 0, 0)
e <- emphasis:::.augment_tree_bdi(tree = brts, pars = p8, model_bin = c(0L,0L,0L),
                                  sample_size = 5L, max_missing = 1e4L, link = 0L, rho = 1)
es <- list(trees = e$trees, weights = rep(1, length(e$trees)), rejected = 0L,
           rejected_overruns = 0L, rejected_lambda = 0L,
           rejected_zero_weights = 0L, time = 0, fhat = 0)
t1 <- tryCatch(emphasis:::m_cpp(es, init_pars = c(0.5, 0.1), plugin = "rpd1",
                                lower_bound = c(0,0), upper_bound = c(3,3),
                                xtol_rel = 1e-3, num_threads = 1L),
               error = function(e) conditionMessage(e))
cat("  m_cpp compact pars ->", substr(t1, 1, 110), "\n")
t2 <- tryCatch(emphasis:::m_cpp(utils::modifyList(es, list(weights = rep(1, length(e$trees) - 1L))),
                                init_pars = p8, plugin = "rpd1",
                                lower_bound = rep(0, 8), upper_bound = rep(3, 8),
                                xtol_rel = 1e-3, num_threads = 1L),
               error = function(e) conditionMessage(e))
cat("  m_cpp short weights ->", substr(t2, 1, 110), "\n")
t3 <- tryCatch(emphasis:::em_cpp(brts = brts, init_pars = c(0.5, 0.1), sample_size = 5L,
                                 maxN = 500L, max_missing = 1e4, max_lambda = 1e6,
                                 lower_bound = c(0,0), upper_bound = c(3,3),
                                 xtol_rel = 1e-3, num_threads = 1L, copy_trees = FALSE),
               error = function(e) conditionMessage(e))
cat("  em_cpp compact pars ->", substr(t3, 1, 110), "\n")

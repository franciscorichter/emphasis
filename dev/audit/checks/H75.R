# H75: numeric branching-time vector as `tree` in simulate_tree() conditional path.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages(library(emphasis))

brts <- c(5, 3, 1)
cat("== 1. Documented call: simulate_tree(tree = c(5,3,1), pars = c(0.5,0.1)) [method default bdi]\n")
r1 <- tryCatch(simulate_tree(tree = brts, pars = c(0.5, 0.1)), error = function(e) e)
print(if (inherits(r1, "error")) paste("ERROR:", conditionMessage(r1)) else str(r1))

cat("\n== 2. Same with method = 'thinning'\n")
r2 <- tryCatch(simulate_tree(tree = brts, pars = c(0.5, 0.1), method = "thinning"), error = function(e) e)
print(if (inherits(r2, "error")) paste("ERROR:", conditionMessage(r2)) else str(r2))

cat("\n== 3. Same with useDDD = FALSE (no L-table needed downstream)\n")
r3 <- tryCatch(simulate_tree(tree = brts, pars = c(0.5, 0.1), useDDD = FALSE), error = function(e) e)
print(if (inherits(r3, "error")) paste("ERROR:", conditionMessage(r3)) else str(r3))

cat("\n== 4. Same with n_trees = 5\n")
r4 <- tryCatch(simulate_tree(tree = brts, pars = c(0.5, 0.1), n_trees = 5L), error = function(e) e)
print(if (inherits(r4, "error")) paste("ERROR:", conditionMessage(r4)) else str(r4))

cat("\n== 5. Direct: .extract_Ltable(c(5,3,1)) and .extract_brts(c(5,3,1))\n")
print(tryCatch(emphasis:::.extract_Ltable(brts), error = function(e) paste("ERROR:", conditionMessage(e))))
print(emphasis:::.extract_brts(brts))

cat("\n== 6. Control: the same numeric vector works in the internal samplers and estimate_rates\n")
a1 <- tryCatch(emphasis:::.augment_tree_internal(brts, pars = c(0.5, 0.1), sample_size = 3L), error = function(e) e)
cat("thinning internal: ", if (inherits(a1, "error")) conditionMessage(a1) else sprintf("%d trees, logg=%s", length(a1$trees), paste(round(a1$logg, 3), collapse = ",")), "\n")
a2 <- tryCatch(emphasis:::.augment_tree_bdi(brts, pars = c(0.5, 0.1), model_bin = c(0L,0L,0L), sample_size = 3L), error = function(e) e)
cat("bdi internal:      ", if (inherits(a2, "error")) conditionMessage(a2) else sprintf("%d trees, logg=%s", length(a2$trees), paste(round(a2$logg, 3), collapse = ",")), "\n")

brts20 <- sort(ape::branching.times(ape::rlineage(0.5, 0.1, Tmax = 6) |> (\(t) tryCatch(ape::drop.fossil(t), error = function(e) NULL))()), decreasing = TRUE)
if (length(brts20) < 5) brts20 <- sort(runif(19, 0, 6), decreasing = TRUE)
set.seed(1)
fit <- tryCatch(estimate_rates(tree = brts20, model = "cr", init_pars = c(0.5, 0.1),
                               control = list(max_iter = 3L, max_time = 30, num_threads = 1L,
                                              lower_bound = c(0.01, 0.001), upper_bound = c(2, 1))),
                error = function(e) e)
cat("estimate_rates(numeric brts): ", if (inherits(fit, "error")) paste("ERROR:", conditionMessage(fit)) else paste("OK, pars =", paste(round(fit$pars, 3), collapse = ",")), "\n")

cat("\n== 7. Control: phylo input on the same path works\n")
phy <- ape::rphylo(10, 0.5, 0.1)
r7 <- tryCatch(simulate_tree(tree = phy, pars = c(0.5, 0.1)), error = function(e) e)
cat(if (inherits(r7, "error")) paste("ERROR:", conditionMessage(r7)) else sprintf("OK: tas has %d tips, log_q = %.3f", ape::Ntip(r7$tas), r7$log_q), "\n")

.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths())); library(emphasis)
set.seed(1); tr <- ape::rlineage(0.4, 0.1, Tmax = 8); tr <- ape::drop.fossil(tr)
while (ape::Ntip(tr) < 12 || ape::Ntip(tr) > 40) { tr <- ape::drop.fossil(ape::rlineage(0.4,0.1,Tmax=8)) }
brts <- sort(ape::branching.times(tr), decreasing = TRUE)
cat("tips:", ape::Ntip(tr), "\n")
model <- c(1L,0L,0L)  # ~N (dd)
lb <- emphasis:::.expand_pars(c(0, -0.1, 0, -0.1), model); ub <- emphasis:::.expand_pars(c(2, 0.1, 1, 0.1), model)
ip <- emphasis:::.expand_pars(c(0.5, -0.01, 0.1, 0), model)
run <- function(xtol, reps) {
  codes <- character(reps)
  for (r in seq_len(reps)) {
    res <- tryCatch(emphasis:::em_cpp(brts, ip, 50L, 2000L, 1e4, 1e6, lb, ub, xtol, 1L, FALSE, model, 0L, 1.0, NULL),
                    error = function(e) conditionMessage(e))
    codes[r] <- if (is.character(res)) paste0("ERR:", res) else paste0("nlopt=", res$nlopt)
  }
  print(table(codes))
}
cat("--- xtol_rel = 1e-3\n"); run(1e-3, 30)
cat("--- xtol_rel = 1e-14 (provoke roundoff)\n"); run(1e-14, 30)
cat("--- xtol_rel = 0 (provoke roundoff)\n"); run(0, 10)

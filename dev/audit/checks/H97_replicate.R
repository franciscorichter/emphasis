# H97 replication: (1) does the unused-argument error hold for every method and
# every block? (2) is the verifier's "bounds into control" fix a deterministic
# pass, or does it depend on the clock-seeded simulate_tree() giving a usable
# tree? (3) does the documented interface work on a fixed brts vector?
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages(library(emphasis))

cat("== (1) unused-argument error is independent of method / model / link ==\n")
brts <- c(4.5, 3.2, 2.7, 1.9, 1.4, 1.1, 0.8, 0.5, 0.3, 0.1)
for (m in c("mcem", "cem", "gam")) {
  r <- tryCatch(estimate_rates(brts, method = m, model = "cr",
                               lower_bound = c(0, 0), upper_bound = c(2, 1),
                               control = list(sample_size = 5, max_iter = 1)),
                error = function(e) conditionMessage(e))
  cat(sprintf("  method=%-4s -> %s\n", m, if (is.character(r)) r else "NO ERROR"))
}
r <- tryCatch(estimate_rates(brts, method = "mcem", model = "ep", link = "exponential",
                             lower_bound = c(-5,-5,-5,-5), upper_bound = c(2,2,2,2),
                             control = list(sample_size = 5)),
              error = function(e) conditionMessage(e))
cat(sprintf("  ep/exponential (block 5) -> %s\n", if (is.character(r)) r else "NO ERROR"))
cat("  match.call check: is the error raised before any body code? ->",
    tryCatch({ estimate_rates(NULL, lower_bound = 1); "no" },
             error = function(e) grepl("unused argument", conditionMessage(e))), "\n")

cat("\n== (2) do the test blocks' simulate_tree() calls give a tree with $tes? (20 reps each) ==\n")
sim_specs <- list(
  CR = quote(simulate_tree(pars = c(0.5, 0.1), max_t = 5, model = "cr")),
  DD = quote(simulate_tree(pars = c(0.5, -0.005, 0.1, 0), max_t = 8, model = "dd")),
  EP = quote(simulate_tree(pars = c(0.5, 0.05, 0.1, 0.01), max_t = 5, model = "ep"))
)
for (nm in names(sim_specs)) {
  st <- character(20); has_tes <- logical(20); ntip <- integer(20)
  for (i in 1:20) {
    tr <- eval(sim_specs[[nm]])
    st[i] <- tr$status
    has_tes[i] <- !is.null(tr$tes) && inherits(tr$tes, "phylo")
    ntip[i] <- if (has_tes[i]) ape::Ntip(tr$tes) else NA_integer_
  }
  cat(sprintf("  %s: status table = %s | has $tes phylo: %d/20 | ntip range %s\n",
              nm, paste(names(table(st)), table(st), sep = ":", collapse = " "),
              sum(has_tes), paste(range(ntip, na.rm = TRUE), collapse = "-")))
}

cat("\n== (3) documented interface (bounds in control) on a fixed brts vector, 3 reps ==\n")
for (i in 1:3) {
  fit <- tryCatch(estimate_rates(brts, method = "mcem", model = "cr",
                                 control = list(lower_bound = c(0, 0), upper_bound = c(2, 1),
                                                sample_size = 20, max_iter = 3, max_time = 60,
                                                num_threads = 1)),
                  error = function(e) e)
  if (inherits(fit, "error")) cat("  rep", i, "ERROR:", conditionMessage(fit), "\n") else
    cat(sprintf("  rep %d: pars=%s loglik=%.3f\n", i,
                paste(names(fit$pars), round(fit$pars, 3), sep = "=", collapse = ","), fit$loglik))
}

cat("\n== (4) the 'fixed' CR block as the verifier wrote it (simulate + fit), 10 reps ==\n")
out <- character(10)
for (i in 1:10) {
  out[i] <- tryCatch({
    tr <- simulate_tree(pars = c(0.5, 0.1), max_t = 5, model = "cr")
    fit <- estimate_rates(tr, method = "mcem", model = "cr",
                          control = list(lower_bound = c(0, 0), upper_bound = c(2, 1),
                                         sample_size = 20, max_iter = 3, max_time = 60,
                                         num_threads = 1))
    sprintf("ok (ntip=%d)", ape::Ntip(tr$tes))
  }, error = function(e) paste("ERROR:", substr(conditionMessage(e), 1, 90)))
}
print(table(out))

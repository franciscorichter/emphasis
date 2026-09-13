.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths())); library(emphasis)
set.seed(7)
P <- rbind(c(0.6,0.1), c(0.5,0.1), c(NA,0.1), c(Inf,0.1), c(-1,-1))
for (i in seq_len(nrow(P))) cat("row", i, paste(P[i,],collapse=","), "->", tryCatch({r <- simulate_tree(pars = P[i,], max_t=3, model="cr"); paste("ok status", r$status)}, error=function(e) paste("ERROR:", conditionMessage(e))), "\n")
# batch under mclapply with a bad row that raises an R error: use max_tries as a per-row trigger? no; inject via pars of wrong type
cat("batch (num_threads=2) ->", tryCatch({b <- simulate_tree(pars = P, max_t=3, model="cr", num_threads=2L); paste("ok surv", b$survival_prob)}, error=function(e) paste("ERROR:", conditionMessage(e))), "\n")
# reproduce the vapply failure mode on a try-error element as mclapply would return it
sims <- list(list(survival_prob=1), structure("boom", class="try-error"))
cat("vapply on try-error row ->", tryCatch(mean(vapply(sims, `[[`, 0.0, "survival_prob")), error=function(e) paste("ERROR:", conditionMessage(e))), "\n")
# does mclapply return try-error on a child error? confirm
x <- parallel::mclapply(1:2, function(i) if (i==2) stop("child fail") else 1, mc.cores=2)
cat("mclapply child error class:", class(x[[2]]), "\n")
# conditional batch shape
sim <- NULL; while (is.null(sim$tes) || ape::Ntip(sim$tes) < 10 || ape::Ntip(sim$tes) > 40) sim <- simulate_tree(pars = c(0.6,0.1), max_t=5, model="cr")
cb <- simulate_tree(tree=sim, pars=P[1:2,], model="cr")
cat("conditional batch: class", class(cb), "len", length(cb), "names[[1]]", paste(names(cb[[1]]),collapse=","), "\n")

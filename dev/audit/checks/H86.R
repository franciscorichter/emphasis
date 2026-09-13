.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages(library(emphasis))
set.seed(1)
# Small CEM fit so the details object has the real structure
phy <- ape::rlineage(0.5, 0.1, Tmax = 4)
tes <- ape::drop.fossil(phy)
while (ape::Ntip(tes) < 8 || ape::Ntip(tes) > 40) {
  phy <- ape::rlineage(0.5, 0.1, Tmax = 4); tes <- ape::drop.fossil(phy)
}
cat("tips:", ape::Ntip(tes), "\n")
fit <- estimate_rates(tes, method = "cem", model = "cr",
  control = list(lower_bound = c(0.05, 0), upper_bound = c(2, 1),
                 max_iter = 2, num_particles = 6, num_trees = 3,
                 max_time = 60, num_threads = 1, verbose = FALSE))
cat("class(fit$details$best_IS$ESS):", class(fit$details$best_IS$ESS), "\n")

# (a) visible vs invisible return
r <- withVisible(diagnose_cem(fit, plot = FALSE))
cat("diagnose_cem visible:", r$visible, "\n")
sink_out <- capture.output(print(diagnose_cem(fit, plot = FALSE)))
cat("print has 'B=200' line:", any(grepl("B=200", sink_out)), "\n")

# (b) fit lacking ESS
fit2 <- fit
fit2$details$best_IS$ESS <- NULL
res <- tryCatch({ diagnose_cem(fit2, plot = FALSE); "OK" },
                error = function(e) paste("ERROR:", conditionMessage(e)))
cat("diagnose_cem without ESS ->", res, "\n")

# (c) ESS from fit vs recomputed from lw
lw <- fit$details$best_IS$lw
ess_re <- emphasis:::.ess_from_lw(lw)
cat("ESS in fit:", fit$details$best_IS$ESS, " recomputed:", ess_re, "\n")

# (d) diagnose_mcem par_names vs par[0-9]+ columns
fitm <- estimate_rates(tes, method = "mcem", model = "cr",
  control = list(lower_bound = c(0.05, 0), upper_bound = c(2, 1),
                 max_iter = 2, sample_size = 20, max_time = 60,
                 num_threads = 1, verbose = FALSE))
pc <- grep("^par[0-9]+$", names(fitm$details$mcem), value = TRUE)
cat("mcem par cols:", length(pc), " names(x$pars):", length(fitm$pars), "\n")
rm <- withVisible(diagnose_mcem(fitm, plot = FALSE))
cat("diagnose_mcem visible:", rm$visible, "\n")
fitm2 <- fitm; fitm2$details$final_IS$ESS <- NULL
res2 <- tryCatch({ diagnose_mcem(fitm2, plot = FALSE); "OK" },
                 error = function(e) paste("ERROR:", conditionMessage(e)))
cat("diagnose_mcem without ESS ->", res2, "\n")

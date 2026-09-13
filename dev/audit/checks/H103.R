.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths())); library(emphasis)
cat("mgcv installed here:", requireNamespace("mgcv", quietly = TRUE), "\n")
for (f in c("train_GAM", "train_likelihood_GAM", "predict_survival", ".run_gam",
            "estimate_rates", "diagnose_gam", "emphasis_pipeline", "auto_bounds")) {
  d <- deparse(get(f, envir = asNamespace("emphasis")))
  cat(sprintf("%-22s calls mgcv:: %-5s  requireNamespace('mgcv') guard: %s\n", f,
              any(grepl("mgcv::", d)), any(grepl("requireNamespace\\(\"mgcv\"", d))))
}
# Emulate a machine without mgcv: block the namespace loader for mgcv only.
unloadNamespace("mgcv")
trace(loadNamespace, quote(if (identical(as.character(package), "mgcv")) stop("no package called 'mgcv'")), print = FALSE)
cat("requireNamespace('mgcv') now:", requireNamespace("mgcv", quietly = TRUE), "\n")
n <- 50
pars_mat <- cbind(beta_0 = runif(n, 0.2, 1.0), gamma_0 = runif(n, 0.05, 0.3))
sims <- lapply(rbinom(n, 1, 0.6), function(s) list(status = if (s == 1) "done" else "extinct"))
r <- tryCatch(train_GAM(sims, pars_mat, model = "cr"), error = function(e) conditionMessage(e))
cat("train_GAM without mgcv ->", if (is.character(r)) r else class(r)[1], "\n")

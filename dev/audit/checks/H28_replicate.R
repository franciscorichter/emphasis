# H28 replication: vary link (exponential), model (dd), stages (bounds+cem),
# verbose = TRUE (capture what the pipeline itself prints on the fallback path).
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({library(emphasis); library(ape)})
set.seed(7)
tr <- drop.fossil(rlineage(1, 0, Tmax = 4))
while (Ntip(tr) < 40 || Ntip(tr) > 60) tr <- drop.fossil(rlineage(1, 0, Tmax = 4))
tr$edge.length <- tr$edge.length / max(branching.times(tr)) * 1e-3
cat("n_tips =", Ntip(tr), " crown age =", max(branching.times(tr)), "\n")

run <- function(model, link, stages, verbose) {
  warns <- character()
  out <- capture.output(
    res <- withCallingHandlers(
      emphasis_pipeline(tr, model = model, link = link, stages = stages,
        control = list(num_threads = 1L, max_time = 60,
                       gam = list(n_grid = 10, sample_size = 5),
                       cem = list(max_iter = 2, num_particles = 8, num_trees = 2,
                                  sample_size = 1, maxN = 5)),
        verbose = verbose),
      warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") }))
  list(res = res, warns = warns, out = out)
}

cat("\n=== V1: dd / exponential / bounds+gam, verbose = FALSE ===\n")
v1 <- run("dd", "exponential", c("bounds", "gam"), FALSE)
print(v1$warns); print(v1$res$log[, c("stage", "status", "loglik")])
cat("cond =", v1$res$cond, " survival_gam NULL:", is.null(v1$res$bounds$survival_gam),
    " center NULL:", is.null(v1$res$bounds$center), "\n")
wide <- emphasis:::.wide_bounds(emphasis:::.resolve_model("dd"), 1L, max(branching.times(tr)), Ntip(tr))
cat("bounds == wide box:", isTRUE(all.equal(unname(v1$res$bounds$lower_bound), wide$lb)),
    isTRUE(all.equal(unname(v1$res$bounds$upper_bound), wide$ub)), "\n")

cat("\n=== V2: cr / linear / bounds+cem, verbose = TRUE (what does the pipeline print?) ===\n")
v2 <- run("cr", "linear", c("bounds", "cem"), TRUE)
print(v2$warns); print(v2$res$log[, c("stage", "status", "loglik")])
cat("cond =", v2$res$cond, " fits$cem$cond =", if (!is.null(v2$res$fits$cem)) v2$res$fits$cem$cond else NA, "\n")
cat("--- verbose output lines mentioning infeasib/feasible/center/bounds/cond/warn: ---\n")
print(grep("infeas|feasib|center|bounds\\]|cond|warn|survival", v2$out, value = TRUE, ignore.case = TRUE))
cat("--- verbose lines between '[Stage 1]' and '[Stage' (the whole bounds stage transcript): ---\n")
i1 <- grep("\\[Stage 1\\]", v2$out); i2 <- grep("\\[Stage [2-4]\\]", v2$out)
print(v2$out[i1:(if (length(i2)) i2[1] else length(v2$out))])

cat("\n=== V3: is result$cond == FALSE a unique fingerprint of the fallback when 'bounds' ran? ===\n")
cat("pipeline never passes train_surv_gam to auto_bounds:",
    !grepl("train_surv_gam", paste(deparse(emphasis_pipeline), collapse = "")), "\n")
cat("auto_bounds success path sets surv_gam only under train_surv_gam (default TRUE) => cond NULL iff fallback\n")

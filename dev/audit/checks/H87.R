.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages(library(emphasis))
mk <- function(ll, np, var = 0.05, model = c(0L, 0L, 0L)) {
  structure(list(pars = c(0.5, 0.1), loglik = ll, loglik_var = var,
                 n_pars = np, AIC = -2 * ll + 2 * np, method = "mcem",
                 model = model, cond = FALSE, details = NULL),
            class = "emphasis_fit")
}
try_it <- function(label, expr) {
  res <- tryCatch({ x <- expr; paste("OK; nrow =", nrow(x)) },
                  error = function(e) paste("ERROR:", conditionMessage(e)),
                  warning = function(w) paste("WARNING:", conditionMessage(w)))
  cat(sprintf("%-45s -> %s\n", label, res))
}
# 1. two unnamed CR fits, both with loglik_var  (duplicate label "CR")
try_it("two CR fits, integer n_pars, with var",
       emphasis:::compare_models(mk(-10, 2L), mk(-12, 2L)))
# 2. same but no variance (pairwise block skipped)
try_it("two CR fits, integer n_pars, no var",
       emphasis:::compare_models(mk(-10, 2L, NA), mk(-12, 2L, NA)))
# 3. double n_pars (what pipeline / a user-built fit might carry)
try_it("two fits, double n_pars",
       emphasis:::compare_models(mk(-10, 2, NA, c(0L,0L,0L)), mk(-12, 4, NA, c(1L,0L,0L))))
# 4. via list(a, b) as the hypothesis phrases it
try_it("emphasis:::compare_models(list(a, b))",
       emphasis:::compare_models(list(mk(-10, 2L), mk(-12, 2L))))
# 5. what estimate_rates actually stores for n_pars
cat("estimate_rates stores n_pars as length(pars) -> typeof:",
    typeof(length(c(a = 1, b = 2))), "\n")

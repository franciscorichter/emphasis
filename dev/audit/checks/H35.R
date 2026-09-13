# H35: does n_pars count parameters fixed via lb == ub?
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths())); library(emphasis)
set.seed(5)
tr <- ape::rcoal(20); brts <- sort(ape::branching.times(tr), decreasing = TRUE); brts <- brts/max(brts)*5
fit <- suppressWarnings(estimate_rates(brts, model = "dd", link = "linear", method = "mcem",
        control = list(lower_bound = c(0, -0.1, 0, 0), upper_bound = c(1, 0, 0.5, 0),   # gammaN fixed at 0
                       sample_size = 20L, max_iter = 2L, max_time = 60)))
print(fit$pars); cat("n_pars:", fit$n_pars, " free (lb<ub):", sum(c(0,-0.1,0,0) != c(1,0,0.5,0)),
    " AIC:", fit$AIC, " -2ll+2*3:", -2*fit$loglik + 6, "\n")

.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths())); library(emphasis)
# Probe the lam+mu == 0 path in the forward simulator under linear link + D covariate (all rates clipped to 0 at new t)
set.seed(7)
st <- table(unlist(lapply(1:200, function(i) {
  s <- tryCatch(simulate_tree(pars = c(0.8, -0.6, 0.0, 0.0), max_t = 3, model = "d", max_lin = 2000, max_tries = 1), error = function(e) "error")
  if (identical(s, "error")) "error" else s$status
})))
print(st)
st2 <- table(unlist(lapply(1:200, function(i) {
  s <- tryCatch(simulate_tree(pars = c(0.8, -0.6, 0.5, -0.6), max_t = 3, model = "d", max_lin = 2000, max_tries = 1), error = function(e) "error")
  if (identical(s, "error")) "error" else s$status
})))
print(st2)

# Captures the cr/dd numbers that must not move when the pendant-age covariates
# are given the observed topology.  Run before and after the change and diff.
suppressWarnings(suppressMessages(devtools::load_all(".", quiet = TRUE)))
out <- commandArgs(trailingOnly = TRUE)[1]

set.seed(42)
phy  <- ape::rphylo(12, 0.5, 0)
brts <- sort(ape::branching.times(phy), decreasing = TRUE)

# Fixed augmented trees, captured once and reused by both runs.
tf <- "dev/pendant_baseline_trees.rds"
if (!file.exists(tf)) {
  a_cr <- augment_trees(as.numeric(brts), c(0.6, 0, 0, 0, 0.2, 0, 0, 0), 5L, 2000L,
                        1000L, 1e6, 1L, model = c(0L, 0L, 0L), link = 0L, rho = 1)
  a_dd <- augment_trees(as.numeric(brts), c(0.8, -0.02, 0, 0, 0.2, 0, 0, 0), 5L, 2000L,
                        1000L, 1e6, 1L, model = c(1L, 0L, 0L), link = 0L, rho = 1)
  saveRDS(list(cr = a_cr$trees, dd = a_dd$trees), tf)
}
trs <- readRDS(tf)

res <- list()
res$logf_cr     <- eval_logf(c(0.6, 0, 0, 0, 0.2, 0, 0, 0), trs$cr,
                             model = c(0L, 0L, 0L), link = 0L, rho = 1)
res$logf_dd     <- eval_logf(c(0.8, -0.02, 0, 0, 0.2, 0, 0, 0), trs$dd,
                             model = c(1L, 0L, 0L), link = 0L, rho = 1)
res$logf_cr_exp <- eval_logf(c(-0.5, 0, 0, 0, -1.6, 0, 0, 0), trs$cr,
                             model = c(0L, 0L, 0L), link = 1L, rho = 1)
res$logf_dd_exp <- eval_logf(c(-0.2, -0.01, 0, 0, -1.6, 0, 0, 0), trs$dd,
                             model = c(1L, 0L, 0L), link = 1L, rho = 1)

fit <- function(seed, model, init, lb, ub) {
  set.seed(seed)
  f <- estimate_rates(phy, method = "mcem", model = model, init_pars = init,
                      control = list(lower_bound = lb, upper_bound = ub,
                                     num_trees = 30L, max_iter = 5L, verbose = FALSE))
  list(pars = unname(f$pars), loglik = f$loglik, iterations = f$iterations,
       stop_reason = f$stop_reason)
}
res$fit_cr <- fit(7, "cr", c(0.5, 0.1), c(0, 0), c(2, 1))
res$fit_dd <- fit(9, "dd", c(0.8, -0.02, 0.2, 0),
                  c(0.1, -0.5, 0, -0.01), c(3, 0.01, 1, 0.01))

es <- list(trees = trs$dd, weights = rep(1, length(trs$dd)), rejected = 0L,
           rejected_overruns = 0L, rejected_lambda = 0L,
           rejected_zero_weights = 0L, time = 0, fhat = 0)
res$m_dd <- m_cpp(e_step = es, init_pars = c(0.8, -0.02, 0, 0, 0.2, 0, 0, 0),
                  plugin = "rpd1",
                  lower_bound = c(0.1, -0.5, 0, 0, 0, -0.01, 0, 0),
                  upper_bound = c(3, 0.01, 0, 0, 1, 0.01, 0, 0),
                  xtol_rel = 1e-4, num_threads = 1L,
                  model = c(1L, 0L, 0L), link = 0L, rho = 1)$estimates

saveRDS(res, out)
cat("wrote", out, "\n")

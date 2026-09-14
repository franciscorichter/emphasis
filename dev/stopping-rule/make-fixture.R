# Freeze the scoring corpus the stopping-rule test rescores against.
# Run once, from the package root, on a build whose scorers are believed
# correct; the test then fails if eval_logf or the M-step ever move.
devtools::load_all(".", quiet = TRUE)

mk_corpus <- function(model, link, pars, brts, n = 60L, seed = 1L) {
  augment_trees(brts = brts, pars = pars, sample_size = n, maxN = 40L * n,
                max_missing = 1e4, max_lambda = 1e6, num_threads = 1L,
                model = model, link = link, rho = 1,
                parent_tip_start = numeric(0), seed = seed)
}

set.seed(1)
sim_cr <- simulate_tree(pars = c(0.6, 0.15), max_t = 8, model = "cr")
brts_cr <- sort(ape::branching.times(sim_cr$tes), decreasing = TRUE)
set.seed(2)
sim_dd <- simulate_tree(pars = c(0.8, -0.0175, 0.1, 0), max_t = 10, model = "dd")
brts_dd <- sort(ape::branching.times(sim_dd$tes), decreasing = TRUE)
cat("cr tips:", length(brts_cr) + 1, " dd tips:", length(brts_dd) + 1, "\n")

defs <- list(
  cr_lin = list(model = c(0L, 0L, 0L), link = 0L, brts = brts_cr,
                pars = .expand_pars(c(0.5, 0.1), c(0L, 0L, 0L)),
                lower = .expand_pars(c(0, 0), c(0L, 0L, 0L)),
                upper = .expand_pars(c(3, 3), c(0L, 0L, 0L))),
  cr_exp = list(model = c(0L, 0L, 0L), link = 1L, brts = brts_cr,
                pars = .expand_pars(c(log(0.5), log(0.1)), c(0L, 0L, 0L)),
                lower = .expand_pars(c(-5, -5), c(0L, 0L, 0L)),
                upper = .expand_pars(c(2, 2), c(0L, 0L, 0L))),
  dd_lin = list(model = c(1L, 0L, 0L), link = 0L, brts = brts_dd,
                pars = .expand_pars(c(0.8, -0.0175, 0.1, 0), c(1L, 0L, 0L)),
                lower = .expand_pars(c(1e-3, -0.0875, 0, 0), c(1L, 0L, 0L)),
                upper = .expand_pars(c(4, 0, 4, 0), c(1L, 0L, 0L))),
  dd_exp = list(model = c(1L, 0L, 0L), link = 1L, brts = brts_dd,
                pars = .expand_pars(c(log(0.8), -0.01, log(0.1), 0), c(1L, 0L, 0L)),
                lower = .expand_pars(c(-5, -0.1, -5, 0), c(1L, 0L, 0L)),
                upper = .expand_pars(c(2, 0, 2, 0), c(1L, 0L, 0L))))

trees <- list(); cases <- list(); msteps <- list()
for (nm in names(defs)) {
  d <- defs[[nm]]
  cp <- mk_corpus(d$model, d$link, d$pars, d$brts, n = 60L, seed = 20250914L)
  trees[[nm]] <- cp$trees
  cat(nm, ": ", length(cp$trees), " trees\n", sep = "")
  set.seed(99)
  for (k in 1:5) {
    th <- d$pars
    free <- which(d$lower != d$upper)
    th[free] <- th[free] + (k - 3) * 0.05 * pmax(abs(th[free]), 0.05)
    th <- pmin(pmax(th, d$lower), d$upper)
    ev <- eval_logf(th, cp$trees, model = d$model, link = d$link, rho = 1)
    cases[[length(cases) + 1L]] <- list(corpus = nm, theta = th,
                                        model = d$model, link = d$link, rho = 1,
                                        logf = ev$logf, logg = ev$logg)
    lw <- cp$logf - cp$logg
    w <- exp(lw - max(lw)) * (1 + 0.1 * ((seq_along(lw) + k) %% 3))
    m <- m_cpp(e_step = list(trees = cp$trees, weights = w, rejected = 0L,
                             rejected_overruns = 0L, rejected_lambda = 0L,
                             rejected_zero_weights = 0L, time = 0, fhat = 0),
               init_pars = th, plugin = "rpd1", lower_bound = d$lower,
               upper_bound = d$upper, xtol_rel = 1e-3, num_threads = 1L,
               model = d$model, link = d$link, rho = 1, rconditional = NULL)
    msteps[[length(msteps) + 1L]] <- list(corpus = nm, weights = w, init = th,
                                          lower = d$lower, upper = d$upper,
                                          model = d$model, link = d$link, rho = 1,
                                          estimates = as.numeric(m$estimates))
  }
}
dir.create("tests/testthat/fixtures", showWarnings = FALSE, recursive = TRUE)
saveRDS(list(trees = trees, cases = cases, msteps = msteps,
             built = list(r = as.character(getRversion()), when = Sys.Date())),
        "tests/testthat/fixtures/scoring-corpus.rds", version = 2)
cat("fixture written:", length(cases), "eval_logf cases,", length(msteps),
    "M-step cases\n")

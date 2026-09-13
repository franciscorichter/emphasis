## H47 replicate, part 2 -- elite re-evaluation under frozen forked streams,
## traced through the real resampler and the real CEM driver. ~1 min.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({ library(emphasis); library(ape) })
same <- function(a, b) isTRUE(all.equal(a, b))

set.seed(11)
phy <- ape::drop.fossil(ape::rlineage(0.5, 0.15, Tmax = 5))
while (ape::Ntip(phy) < 10 || ape::Ntip(phy) > 40)
  phy <- ape::drop.fossil(ape::rlineage(0.5, 0.15, Tmax = 5))
brts <- sort(as.numeric(ape::branching.times(phy)), decreasing = TRUE)

P <- rbind(c(0.7, -0.006, 0, 0, 0.20, 0, 0, 0),   # best of the four last time
           c(0.6, -0.005, 0, 0, 0.15, 0, 0, 0),
           c(0.5, -0.004, 0, 0, 0.10, 0, 0, 0),
           c(0.4, -0.002, 0, 0, 0.05, 0, 0, 0))
mkpop <- function(P) list(pars = as.data.frame(P), fhat = rep(NA_real_, nrow(P)),
                          log_q = vector("list", nrow(P)), trees = vector("list", nrow(P)))
input <- list(brts = brts, sample_size = 10L, maxN = 2000L, max_missing = 1e4,
              max_lambda = 1e6, lower_bound = P[1, ] - 1, upper_bound = P[1, ] + 1,
              shared_trees = FALSE, bias_correct = FALSE,
              model = c(1L, 0L, 0L), link = 0L, rho = 1)
ev <- function(pop) emphasis:::.eval_independent(pop, input, 2L)$pop$fhat
lb8 <- c(0.05, -0.05, 0, 0, 0.01, 0, 0, 0); ub8 <- c(2, 0, 0, 0, 1, 0, 0, 0)
sdv <- (ub8 - lb8) / 4

trace_gens <- function(ngen, label) {
  pop <- mkpop(P); pop$fhat <- ev(pop); prev <- pop
  for (g in 2:ngen) {
    set.seed(100 + g)   # resampler's own R draws, held fixed
    pop <- emphasis:::.resample_particles(prev, 4L, 0.5, lb8, ub8, sdv)
    pop$fhat <- ev(pop)
    row1_same_pars <- same(as.numeric(pop$pars[1, ]), as.numeric(prev$pars[1, ]))
    cat(sprintf("%s gen %d: argmax(prev)=%d; row1 pars unchanged from prev row1: %s; fhat[1] == prev fhat[1]: %s; fhat[1] == prev max: %s\n",
                label, g, which.max(prev$fhat), row1_same_pars,
                same(pop$fhat[1], prev$fhat[1]), same(pop$fhat[1], max(prev$fhat))))
    prev <- pop
  }
}

## ---- unprimed (this process has not touched the engines yet) --------------
trace_gens(4, "U  unprimed")

## ---- prime, then repeat --------------------------------------------------
invisible(emphasis:::augment_trees(brts, P[1, ], 5L, 500L, 1e4L, 1e6, 1L, c(1L,0L,0L), 0L, 1))
trace_gens(4, "P  primed  ")

## ---- BDI sampler, dd model (4 compact pars), set.seed reach ----------------
bdi_dd <- replicate(2, { set.seed(1)
  emphasis:::.augment_tree_bdi(phy, pars = c(0.6, -0.005, 0.15, 0), model_bin = c(1L,0L,0L),
                               sample_size = 10L, link = 0L, rho = 1)$logf }, simplify = FALSE)
cat(sprintf("D  BDI dd (4 pars), set.seed(1) x2 identical logf: %s\n", same(bdi_dd[[1]], bdi_dd[[2]])))

## ---- real CEM driver in this primed session: exact ties of best_loglik -----
cemfit <- function(nt) { set.seed(5)
  estimate_rates(phy, model = "cr", method = "cem",
                 control = list(num_particles = 8L, num_trees = 10L, max_iter = 6L,
                                max_time = 200, num_threads = nt, patience = 10L,
                                lower_bound = c(0.01, 0), upper_bound = c(3, 2))) }
for (nt in c(2L, 1L)) {
  f <- cemfit(nt); bl <- f$details$best_loglik; bp <- f$details$best_pars
  ties <- vapply(2:length(bl), function(k) same(bl[k], bl[k - 1]), NA)
  same_best <- vapply(2:nrow(bp), function(k) same(bp[k, ], bp[k - 1, ]), NA)
  cat(sprintf("R  cem threads=%d: best_loglik = %s; best particle unchanged: %s; exact tie: %s; stop=%s\n",
              nt, paste(round(bl, 4), collapse = " "),
              paste(same_best, collapse = ","), paste(ties, collapse = ","), f$details$converged))
}

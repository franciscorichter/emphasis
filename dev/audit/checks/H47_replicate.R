## H47 replicate -- independent re-run with a different tree, the dd model,
## both links, the REAL CEM dispatcher (.eval_independent / .resample_particles)
## and the real estimate_rates(method = "cem") driver. ~1-2 min.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({ library(emphasis); library(ape) })
same <- function(a, b) isTRUE(all.equal(a, b))

set.seed(11)
phy <- ape::drop.fossil(ape::rlineage(0.5, 0.15, Tmax = 5))
while (ape::Ntip(phy) < 10 || ape::Ntip(phy) > 40)
  phy <- ape::drop.fossil(ape::rlineage(0.5, 0.15, Tmax = 5))
brts <- sort(as.numeric(ape::branching.times(phy)), decreasing = TRUE)
cat(sprintf("Tree: %d tips, crown age %.3f\n", ape::Ntip(phy), brts[1]))

## dd model (lambda = b0 + bN*N), linear link; 4 distinct particles in pars8 layout
P <- rbind(c(0.6, -0.005, 0, 0, 0.15, 0, 0, 0),
           c(0.5, -0.004, 0, 0, 0.10, 0, 0, 0),
           c(0.7, -0.006, 0, 0, 0.20, 0, 0, 0),
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

## ---- U. UNPRIMED parent (this process has not touched the C++ engines) ----
## Everything here forks before any serial augmentation in the parent.
U1 <- ev(mkpop(P)); U2 <- ev(mkpop(P))
cat(sprintf("U  unprimed: two forked generations identical fhat: %s (fhat1 = %s)\n",
            same(U1, U2), paste(round(U1, 3), collapse = " ")))
## identical-pars rows 1&2 (children's first tasks) and 3&4 (second tasks)
Pdup <- P[c(1, 1, 3, 3), ]
Ud <- ev(mkpop(Pdup))
cat(sprintf("U  unprimed, duplicate pars: fhat[1]==fhat[2]: %s; fhat[3]==fhat[4]: %s\n",
            same(Ud[1], Ud[2]), same(Ud[3], Ud[4])))

## ---- prime the parent with ONE serial thinning augmentation ---------------
invisible(emphasis:::augment_trees(brts, P[1, ], 5L, 500L, 1e4L, 1e6, 1L, c(1L,0L,0L), 0L, 1))

## ---- P. PRIMED parent ----------------------------------------------------
P1 <- ev(mkpop(P)); P2 <- ev(mkpop(P))
cat(sprintf("P  primed: two forked generations identical fhat: %s (fhat1 = %s)\n",
            same(P1, P2), paste(round(P1, 3), collapse = " ")))
Pd <- ev(mkpop(Pdup))
cat(sprintf("P  primed, duplicate pars: fhat[1]==fhat[2]: %s; fhat[3]==fhat[4]: %s; fhat[1]==fhat[3]: %s\n",
            same(Pd[1], Pd[2]), same(Pd[3], Pd[4]), same(Pd[1], Pd[3])))
## exponential link variant of the same thing
input_exp <- input; input_exp$link <- 1L
Pexp <- rbind(c(log(0.6), -0.005, 0, 0, log(0.15), 0, 0, 0))[c(1, 1, 1, 1), ]
Pe <- emphasis:::.eval_independent(mkpop(Pexp), input_exp, 2L)$pop$fhat
cat(sprintf("P  primed, exp link, 4 identical particles: fhat = %s; [1]==[2]: %s; [1]==[3]: %s\n",
            paste(round(Pe, 4), collapse = " "), same(Pe[1], Pe[2]), same(Pe[1], Pe[3])))

## ---- Elite re-evaluation defeat, via the real resampler --------------------
gen <- function(pop) { pop$fhat <- ev(pop); pop }
run_gens <- function(ngen = 4) {
  pop <- gen(mkpop(P)); out <- list(pop$fhat)
  for (g in 2:ngen) {
    pop <- emphasis:::.resample_particles(pop, 4L, 0.5, lb8, ub8, sdv)
    pop <- gen(pop); out[[g]] <- pop$fhat
  }
  out
}
H <- run_gens(4)
top_frozen <- vapply(2:4, function(g) same(H[[g]][1], max(H[[g - 1]])), NA)
cat(sprintf("P  primed, 4 generations via .resample_particles: top elite fhat identical to previous gen max: %s\n",
            paste(top_frozen, collapse = ",")))

## ---- D. set.seed reach, dd model, both links (thinning) --------------------
augdd <- function(link, p) { set.seed(1)
  emphasis:::augment_trees(brts, p, 10L, 2000L, 1e4L, 1e6, 1L, c(1L,0L,0L), link, 1)$logf }
d_lin <- replicate(3, augdd(0L, P[1, ]), simplify = FALSE)
d_exp <- replicate(3, augdd(1L, Pexp[1, ]), simplify = FALSE)
cat(sprintf("D  dd thinning, set.seed(1) x3: linear identical %d/2; exp identical %d/2\n",
            same(d_lin[[1]], d_lin[[2]]) + same(d_lin[[1]], d_lin[[3]]),
            same(d_exp[[1]], d_exp[[2]]) + same(d_exp[[1]], d_exp[[3]])))
## BDI sampler with dd model under set.seed (R RNG): should be reproducible
bdi_dd <- tryCatch(replicate(2, { set.seed(1)
  emphasis:::.augment_tree_bdi(phy, pars = c(0.6, -0.005, 0.15), model_bin = c(1L,0L,0L),
                               sample_size = 10L, link = 0L, rho = 1)$logf }, simplify = FALSE),
  error = function(e) conditionMessage(e))
cat(sprintf("D  BDI dd, set.seed(1) x2 identical: %s\n",
            if (is.list(bdi_dd)) same(bdi_dd[[1]], bdi_dd[[2]]) else bdi_dd))

## ---- R. real driver: estimate_rates(method = "cem") in this primed session ---
## Per generation the resampler puts the best elite at row 1 with fhat = NA
## ("always re-evaluate"). With frozen child streams, its re-evaluation should
## reproduce the previous generation's value exactly when threads > 1.
cemfit <- function(nt) { set.seed(5)
  estimate_rates(phy, model = "cr", method = "cem",
                 control = list(num_particles = 8L, num_trees = 10L, max_iter = 4L,
                                max_time = 200, num_threads = nt,
                                lower_bound = c(0.01, 0), upper_bound = c(3, 2))) }
frozen_count <- function(fit) {
  h <- fit$details$hist_fhat_all; h <- h[!vapply(h, is.null, NA)]
  if (length(h) < 2) return(NA)
  sum(vapply(2:length(h), function(g) same(h[[g]][1], max(h[[g - 1]], na.rm = TRUE)), NA))
}
f2 <- cemfit(2L); f1 <- cemfit(1L)
cat(sprintf("R  estimate_rates(cem, 8 particles x 10 trees, 4 iters): generations where re-evaluated top elite == previous max: threads=2: %s of %d; threads=1: %s of %d\n",
            frozen_count(f2), length(Filter(Negate(is.null), f2$details$hist_fhat_all)) - 1,
            frozen_count(f1), length(Filter(Negate(is.null), f1$details$hist_fhat_all)) - 1))
cat(sprintf("R  pars threads=2: %s | threads=1: %s\n",
            paste(signif(f2$pars, 5), collapse = ","), paste(signif(f1$pars, 5), collapse = ",")))

## ---- N. default CEM control never forks (num_trees = 1 < 10) ---------------
ctl <- emphasis:::estimate_rates_control("cem")
cat(sprintf("N  CEM defaults: num_trees = %s, num_threads = %s -> R-level fork requires sample_size >= 10\n",
            ctl$num_trees, ctl$num_threads))

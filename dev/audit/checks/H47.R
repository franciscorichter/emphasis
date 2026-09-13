## H47 — RNG seeding: set.seed() reach, engine identity, fork inheritance.
## Self-contained; ~1 min. num_threads = 1 everywhere except the fork test,
## which is the hypothesis' own subject (mclapply children, not TBB threads).
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({ library(emphasis); library(ape) })

set.seed(2024)
phy  <- ape::rlineage(0.6, 0.2, Tmax = 6)          # birth-death, extinct lineages kept
phy  <- ape::drop.fossil(phy)
while (ape::Ntip(phy) < 10 || ape::Ntip(phy) > 40) {
  phy <- ape::drop.fossil(ape::rlineage(0.6, 0.2, Tmax = 6))
}
brts <- sort(as.numeric(ape::branching.times(phy)), decreasing = TRUE)
cat(sprintf("Tree: %d tips, crown age %.3f\n", ape::Ntip(phy), brts[1]))
pars8 <- c(0.6, 0, 0, 0, 0.2, 0, 0, 0)

same <- function(a, b) isTRUE(all.equal(a, b))

## ---- A. forward simulate_tree() under set.seed(1) -------------------------
fwd <- function() { set.seed(1); simulate_tree(pars = c(0.6, 0.1), max_t = 3, max_tries = 50,
                                               model = "cr", useDDD = FALSE)$L }
fwdL <- replicate(5, fwd(), simplify = FALSE)
a_ident <- sum(vapply(fwdL[-1], function(L) same(L, fwdL[[1]]), NA))
cat(sprintf("A  forward simulate_tree, 5 runs set.seed(1): rows = %s; identical to run 1: %d/4\n",
            paste(vapply(fwdL, NROW, 1L), collapse = ","), a_ident))

## ---- B. emphasis:::augment_trees() (C++ thinning) under set.seed(1) ------------------
aug <- function() { set.seed(1)
  emphasis:::augment_trees(brts, pars8, sample_size = 20L, maxN = 2000L, max_missing = 1e4L,
                max_lambda = 500, num_threads = 1L, model = c(0L,0L,0L), link = 0L, rho = 1) }
augR <- replicate(5, aug(), simplify = FALSE)
b_ident <- sum(vapply(augR[-1], function(r) same(r$logf, augR[[1]]$logf), NA))
cat(sprintf("B  augment_trees (thinning), 5 runs set.seed(1): identical logf to run 1: %d/4; run1 logf[1:3] = %s; run2 = %s\n",
            b_ident, paste(round(augR[[1]]$logf[1:3], 3), collapse = " "),
            paste(round(augR[[2]]$logf[1:3], 3), collapse = " ")))

## ---- C. BDI sampler (R RNG) under set.seed(1) -----------------------------
bdi <- function() { set.seed(1)
  emphasis:::.augment_tree_bdi(phy, pars = c(0.6, 0.2), model_bin = c(0L,0L,0L),
                               sample_size = 20L, link = 0L, rho = 1) }
bdiR <- replicate(5, bdi(), simplify = FALSE)
c_ident <- sum(vapply(bdiR[-1], function(r) same(r$logg, bdiR[[1]]$logg) && same(r$logf, bdiR[[1]]$logf), NA))
cat(sprintf("C  BDI sampler, 5 runs set.seed(1): identical (logg,logf) to run 1: %d/4\n", c_ident))

## ---- D. full emphasis() cr fit, default (bdi) vs thinning, set.seed(1) ---
fit <- function(sampling) { set.seed(1)
  estimate_rates(phy, model = "cr", init_pars = c(0.5, 0.1), method = "mcem",
           control = list(sample_size = 30L, max_iter = 3L, max_time = 120,
                          lower_bound = c(0.01, 0.0), upper_bound = c(3, 2),
                          sampling = sampling, num_threads = 1L)) }
fit_bdi <- replicate(2, fit("bdi"), simplify = FALSE)
fit_thn <- replicate(2, fit("dynamic_fresh"), simplify = FALSE)
cat(sprintf("D  estimate_rates(cr, mcem, bdi): pars run1 = %s | run2 = %s | identical = %s\n",
            paste(signif(fit_bdi[[1]]$pars, 6), collapse = ","),
            paste(signif(fit_bdi[[2]]$pars, 6), collapse = ","),
            same(fit_bdi[[1]]$pars, fit_bdi[[2]]$pars)))
cat(sprintf("D  estimate_rates(cr, mcem, thinning): pars run1 = %s | run2 = %s | identical = %s\n",
            paste(signif(fit_thn[[1]]$pars, 6), collapse = ","),
            paste(signif(fit_thn[[2]]$pars, 6), collapse = ","),
            same(fit_thn[[1]]$pars, fit_thn[[2]]$pars)))

## ---- E. fork inheritance of the thread_local engines ----------------------
## Parent has already called augment_trees (B above) so both thread_local
## engines exist with a definite state. mclapply forks 2 children which each
## receive a byte copy of that state. Prescheduled: child 1 gets tasks 1,3;
## child 2 gets tasks 2,4.
aug_one <- function(i) emphasis:::augment_trees(brts, pars8, 20L, 2000L, 1e4L, 500,
                                     1L, c(0L,0L,0L), 0L, 1)$logf
forked <- parallel::mclapply(1:4, aug_one, mc.cores = 2L, mc.preschedule = TRUE)
cat(sprintf("E  primed parent, mclapply(2 cores): task1==task2 (different children): %s; task1==task3 (same child): %s; task3==task4: %s\n",
            same(forked[[1]], forked[[2]]), same(forked[[1]], forked[[3]]), same(forked[[3]], forked[[4]])))
## Same with the exact CEM worker (.simulate_particle) at two different particles:
sp <- function(i) emphasis:::.simulate_particle(brts, pars8, c(0L,0L,0L), 0L, 20L, 2000L, 1e4, 500, 1L, 1)
spf <- parallel::mclapply(1:2, sp, mc.cores = 2L)
cat(sprintf("E  CEM worker .simulate_particle in 2 forked children, same pars: identical logg: %s\n",
            same(spf[[1]]$logg, spf[[2]]$logg)))

## Un-primed parent: a fresh R process that forks BEFORE ever touching the engines.
code <- '
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages(library(emphasis))
brts <- c(%s); pars8 <- c(0.6,0,0,0,0.2,0,0,0)
f <- function(i) emphasis:::augment_trees(brts, pars8, 20L, 2000L, 1e4L, 500, 1L, c(0L,0L,0L), 0L, 1)$logf
r <- parallel::mclapply(1:2, f, mc.cores = 2L)
cat(isTRUE(all.equal(r[[1]], r[[2]])))
'
tmp <- tempfile(fileext = ".R")
writeLines(sprintf(code, paste(brts, collapse = ",")), tmp)
unprimed <- system2("Rscript", tmp, stdout = TRUE)
cat(sprintf("E  un-primed parent (fresh process), 2 children identical: %s\n", tail(unprimed, 1)))

## ---- F. engine identity on this platform -----------------------------------
sdk <- tryCatch(system2("xcrun", "--show-sdk-path", stdout = TRUE), error = function(e) "")
hdr <- file.path(sdk, "usr/include/c++/v1/__random/default_random_engine.h")
if (file.exists(hdr)) cat("F  libc++:", grep("typedef", readLines(hdr), value = TRUE), "\n")

## ---- E2. cross-call reuse: parent never draws, so every later mclapply
## starts its children from the SAME frozen state (the CEM/GAM loop shape).
forked2 <- parallel::mclapply(1:4, aug_one, mc.cores = 2L, mc.preschedule = TRUE)
cat(sprintf("E2 second mclapply in primed parent: call1.task1==call2.task1: %s; call1.task3==call2.task3: %s\n",
            same(forked[[1]], forked2[[1]]), same(forked[[3]], forked2[[3]])))
## after the parent draws once more (any serial augment_trees), children move on:
invisible(aug_one(1))
forked3 <- parallel::mclapply(1:2, aug_one, mc.cores = 2L)
cat(sprintf("E2 after one serial call in parent: call3.task1==call1.task1: %s; call3.task1==call3.task2: %s\n",
            same(forked3[[1]], forked[[1]]), same(forked3[[1]], forked3[[2]])))

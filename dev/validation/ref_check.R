# NOTE: DDD 5.2.4 pars2 layouts are
#   bd_loglik: c(tdmodel, cond, btorph, verbose, soc)
#   dd_loglik: c(lx, ddmodel, cond, btorph, verbose, soc)
# An earlier version of this file used a 5-slot layout ending (verbose, soc) =
# (2, 0), i.e. it computed stem-age references, which is why the dd - bd offset
# appeared to vary with theta.  See dev/validation/00-design.md section 4.
# Reference-side checks for the validation study: do the exact CR references agree,
# under which conditioning, and how does DDD truncate the DD rates?
suppressPackageStartupMessages(suppressWarnings({ library(ape); library(DDD); library(TreeSim) }))
set.seed(7)
lam <- 0.6; mu <- 0.2; n <- 30
tr <- sim.bd.taxa(n = n, numbsim = 1, lambda = lam, mu = mu, complete = FALSE)[[1]]
brts <- sort(as.numeric(branching.times(tr)), decreasing = TRUE)
cat("crown age", round(brts[1], 3), " n tips", Ntip(tr), "\n\n")

# ape::birthdeath: Nee et al. (1994) likelihood, params (d/b, b-d)
bd <- birthdeath(tr)
r <- bd$para["b-d"]; a <- bd$para["d/b"]
ape_mle <- c(lambda = unname(r/(1 - a)), mu = unname(a*r/(1 - a)))
cat("ape::birthdeath MLE      :", round(ape_mle, 4), "  dev =", round(bd$dev, 4), "\n")

# DDD::bd_ML under the four (cond, btorph) combinations, soc = 2 (crown)
for (cond in 0:1) for (btorph in 0:1) {
  junk <- capture.output(fit <- suppressMessages(bd_ML(brts = brts, cond = cond, btorph = btorph, soc = 2,
                                initparsopt = c(lam, mu), idparsopt = 1:2, tdmodel = 0, verbose = FALSE)))
  cat(sprintf("DDD::bd_ML cond=%d btorph=%d: lambda=%.4f mu=%.4f loglik=%.4f conv=%d\n",
              cond, btorph, fit$lambda0, fit$mu0, fit$loglik, fit$conv))
}
# Which DDD surface does ape maximise? evaluate bd_loglik at ape's MLE for each cond
for (cond in 0:1) cat(sprintf("bd_loglik at ape MLE, cond=%d btorph=1: %.4f\n", cond,
  bd_loglik(pars1 = c(ape_mle[1], ape_mle[2], 0, 0), pars2 = c(0, cond, 1, 0, 2), brts = brts, missnumspec = 0)))
# fixed-parameter loglik surface: bd_loglik vs dd_loglik with huge K (should coincide)
ll_bd <- bd_loglik(pars1 = c(lam, mu, 0, 0), pars2 = c(0, 0, 1, 0, 2), brts = brts, missnumspec = 0)
ll_dd <- dd_loglik(pars1 = c(lam, mu, 1e6), pars2 = c(500, 1, 0, 1, 0, 2), brts = brts, missnumspec = 0)
cat(sprintf("\nbd_loglik(cond0,btorph1) = %.4f   dd_loglik(K=1e6,ddmodel1,cond0,btorph1) = %.4f\n", ll_bd, ll_dd))

# DDD's rate function for ddmodel = 1: does it truncate lambda at 0?
cat("\nDDD:::lambdamu (ddmodel=1) source:\n")
f <- tryCatch(get("lambdamu", asNamespace("DDD")), error = function(e) NULL)
if (!is.null(f)) print(head(deparse(f), 40)) else cat("lambdamu not found; searching...\n")
src <- grep("ddmodel == 1", capture.output(print(DDD:::dd_loglik_rhs)), value = TRUE)
if (length(src)) print(src)

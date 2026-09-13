# Adversarial probes of the H7 envelope fix. Run on POST-FIX (counter needed).
args <- commandArgs(TRUE)
.libPaths(c(args[1], .libPaths()))
suppressMessages(library(emphasis))
aug <- get("augment_trees", envir = asNamespace("emphasis"))
viol <- get("thinning_envelope_violations", envir = asNamespace("emphasis"))

run <- function(label, brts, pars, model, link, rho = 1, n = 2000L,
                max_lambda = 1e6, threads = 1L, max_missing = 500L) {
  viol(reset = TRUE)
  r <- tryCatch(aug(brts = brts, pars = pars, sample_size = n, maxN = 200L * n,
                    max_missing = max_missing, max_lambda = max_lambda,
                    num_threads = threads, model = model, link = link, rho = rho),
                error = function(e) e)
  v <- viol(reset = TRUE)
  if (inherits(r, "error")) {
    cat(sprintf("%-42s ERROR: %s\n", label, conditionMessage(r))); return(invisible(NULL))
  }
  nm <- vapply(r$trees, nrow, numeric(1))
  cat(sprintf("%-42s viol=%-8g ntrees=%-6d rej_lam=%-5d rej_ovr=%-5d rej0w=%-5d mean_rows=%.1f finite_logf=%d\n",
              label, v, length(r$trees), r$rejected_lambda, r$rejected_overruns,
              r$rejected_zero_weights, mean(nm), sum(is.finite(r$logf))))
  invisible(r)
}

brts7  <- c(6, 4.5, 3.0, 2.0, 1.2, 0.5)
brts2  <- c(6)                      # crown only, 2 tips
set.seed(2)
brts30 <- sort(c(10, runif(28, 0.1, 9.9)), decreasing = TRUE)

cat("--- CR boundary values (constant-rate branch) ---\n")
run("cr mu<lam",        brts7, c(0.4,0,0,0,0.2,0,0,0), c(0L,0L,0L), 0L)
run("cr mu==lam",       brts7, c(0.4,0,0,0,0.4,0,0,0), c(0L,0L,0L), 0L)
run("cr mu>lam",        brts7, c(0.3,0,0,0,0.9,0,0,0), c(0L,0L,0L), 0L)
run("cr mu=0",          brts7, c(0.4,0,0,0,0.0,0,0,0), c(0L,0L,0L), 0L)
run("cr lam=0",         brts7, c(0.0,0,0,0,0.2,0,0,0), c(0L,0L,0L), 0L)
run("cr lam=1e-8",      brts7, c(1e-8,0,0,0,0.2,0,0,0), c(0L,0L,0L), 0L)
run("cr rho=0.5",       brts7, c(0.4,0,0,0,0.2,0,0,0), c(0L,0L,0L), 0L, rho = 0.5)
run("cr rho=0.2",       brts7, c(0.4,0,0,0,0.2,0,0,0), c(0L,0L,0L), 0L, rho = 0.2)
run("cr 2-tip tree",    brts2, c(0.4,0,0,0,0.2,0,0,0), c(0L,0L,0L), 0L)
run("cr 30-tip tree",   brts30, c(0.4,0,0,0,0.2,0,0,0), c(0L,0L,0L), 0L, n = 500L)
run("cr sample_size=1", brts7, c(0.4,0,0,0,0.2,0,0,0), c(0L,0L,0L), 0L, n = 1L)
run("cr threads=4",     brts7, c(0.4,0,0,0,0.2,0,0,0), c(0L,0L,0L), 0L, threads = 4L)
run("cr exp link",      brts7, c(log(0.4),0,0,0,log(0.2),0,0,0), c(0L,0L,0L), 1L)
run("cr gauss link",    brts7, c(0.4,0,0,0,0.2,0,0,0), c(0L,0L,0L), 2L)

cat("\n--- dd (N only, constant-rate branch) ---\n")
run("dd lin",  brts7, c(0.6,-0.04,0,0,0.2,0,0,0), c(1L,0L,0L), 0L)
run("dd lin steep", brts7, c(1.2,-0.15,0,0,0.2,0,0,0), c(1L,0L,0L), 0L)
run("dd exp",  brts7, c(log(0.6),-0.05,0,0,log(0.2),0,0,0), c(1L,0L,0L), 1L)
run("dd gauss",brts7, c(0.8,0.1,0,0,0.3,0.05,0,0), c(1L,0L,0L), 2L)

cat("\n--- d / nd (D active -> safety-factor branch) ---\n")
for (bd in c(-2, -1, -0.5, -0.2, 0.2, 0.5, 1, 2)) {
  run(sprintf("nd lin  beta_D=%+.1f", bd), brts7,
      c(0.4,-0.01,0,bd,0.2,0,0,0), c(1L,0L,1L), 0L)
}
for (bd in c(-2, -1, -0.5, 0.5, 1, 2)) {
  run(sprintf("nd exp  beta_D=%+.1f", bd), brts7,
      c(log(0.4),-0.01,0,bd,log(0.2),0,0,0), c(1L,0L,1L), 1L)
}
for (gd in c(-1, -0.5, 0.5, 1, 2)) {
  run(sprintf("nd exp  gamma_D=%+.1f", gd), brts7,
      c(log(0.4),-0.01,0,0,log(0.2),0,0,gd), c(1L,0L,1L), 1L)
}
for (bd in c(-1, -0.5, 0.5, 1)) {
  run(sprintf("d gauss beta_D=%+.1f", bd), brts7,
      c(0.5,0,0,bd,0.2,0,0,0), c(0L,0L,1L), 2L)
}
run("nd exp big tree bD=-1", brts30, c(log(0.4),-0.01,0,-1,log(0.2),0,0,0), c(1L,0L,1L), 1L, n = 300L)
run("nd exp bD=-1 maxlam=500", brts7, c(log(0.4),-0.01,0,-1,log(0.2),0,0,0), c(1L,0L,1L), 1L, max_lambda = 500)

cat("\n--- M active via raw vector (H77 path) ---\n")
for (bm in c(-1, -0.3, 0.3, 1)) {
  run(sprintf("M-active lin beta_M=%+.1f", bm), brts7,
      c(0.4,0,bm,0,0.2,0,0,0), c(0L,1L,0L), 0L)
}
for (bm in c(-1, 1)) {
  run(sprintf("M-active gauss beta_M=%+.1f", bm), brts7,
      c(0.5,0,bm,0,0.2,0,0,0), c(0L,1L,0L), 2L)
}

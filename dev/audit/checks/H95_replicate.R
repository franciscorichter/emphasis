# H95 replicate: vary what the verifier did not.
#  R1: default estimate_rates(method="mcem") dispatch really reaches .mcem_bdi / .augment_tree_bdi
#  R2: the skip() message's own advice (filter='em') still skips -> tests unreachable anywhere
#  R3: CR exactness on a DIFFERENT tree (30 tips, lambda 0.4 / mu 0.2), N = 50 and N = 500,
#      LINEAR and EXPONENTIAL link, vs DDD::bd_loglik(btorph = 1, cond = 0, soc = 2)
#  R4: small tree (12 tips) + near-critical theta (lambda - mu small but > 0)
#  R5: dd on a second tree: BDI vs thinning vs DDD::dd_loglik - lfactorial(n-1)
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressMessages({ library(emphasis); library(ape); library(DDD); library(testthat) })
ns  <- asNamespace("emphasis")
pkg <- "/Users/pancho/Code/emphasis"
set.seed(9595)
ess <- function(lw) { w <- exp(lw - max(lw)); sum(w)^2 / sum(w^2) }

cat("=== R1: does the DEFAULT mcem control reach BDI? ===\n")
cnt <- new.env(); for (f in c(".mcem_bdi", ".augment_tree_bdi", ".mcem_dynamic_fresh")) assign(f, 0L, envir = cnt)
for (f in ls(cnt)) suppressMessages(trace(f, where = ns, print = FALSE,
  tracer = substitute({ assign(FN, get(FN, envir = cnt) + 1L, envir = cnt) }, list(FN = f))))
phy0  <- ape::rphylo(15, 0.5, 0.1)
brts0 <- sort(as.numeric(ape::branching.times(phy0)), decreasing = TRUE)
fit <- tryCatch(estimate_rates(brts0, method = "mcem", model = "cr",
                 control = list(lower_bound = c(0.05, 0.001), upper_bound = c(2, 0.9),
                                max_iter = 3, sample_size = 20, num_threads = 1)),
                error = function(e) { cat("estimate_rates error:", conditionMessage(e), "\n"); NULL })
for (f in ls(cnt)) suppressMessages(untrace(f, where = ns))
print(sapply(ls(cnt), function(f) get(f, envir = cnt)))
if (!is.null(fit)) cat("fit pars:", paste(round(fit$pars, 4), collapse = " "), " loglik:", round(fit$loglik, 4), "\n")

cat("\n=== R2: skip() under the filter the message recommends ===\n")
res <- testthat::test_dir(file.path(pkg, "tests/testthat"), package = "emphasis", load_package = "installed",
                          filter = "em", reporter = "silent", stop_on_failure = FALSE)
df <- as.data.frame(res)
cat(sprintf("filter='em': %d tests, %d skipped, %d passed-expectations, %d failed\n",
            nrow(df), sum(df$skipped), sum(df$passed), sum(df$failed)))
res <- testthat::test_dir(file.path(pkg, "tests/testthat"), package = "emphasis", load_package = "installed",
                          filter = "inference", reporter = "silent", stop_on_failure = FALSE)
df <- as.data.frame(res)
cat(sprintf("filter='inference': %d tests, %d skipped, %d passed-expectations, %d failed\n",
            nrow(df), sum(df$skipped), sum(df$passed), sum(df$failed)))

bdi_fhat <- function(brts, pars, model_bin, N = 200L, rho = 1, link = 0L) {
  r <- ns$.augment_tree_bdi(brts, pars, model_bin = model_bin, sample_size = N, link = link, rho = rho)
  c(fhat = r$fhat, sd_lw = sd(r$weights), ess = ess(r$weights), n = length(r$weights),
    n_nonfinite = sum(!is.finite(r$weights)))
}

cat("\n=== R3: CR exactness on a different tree, both links, N = 50 / 500 ===\n")
phy  <- ape::rphylo(30, 0.4, 0.2)
brts <- sort(as.numeric(ape::branching.times(phy)), decreasing = TRUE)
n <- length(brts) + 1
cat(sprintf("tree: %d tips, crown %.3f\n", n, brts[1]))
grid <- rbind(c(0.4, 0.2), c(0.7, 0.05), c(0.3, 0.28), c(1.5, 1.2), c(0.6, 0))
for (N in c(50L, 500L)) for (link in c(0L, 1L)) {
  out <- t(apply(grid, 1, function(th) {
    pars <- if (link == 0L) th else log(pmax(th, 1e-300))   # exp link: eta = log(rate)
    b <- tryCatch(bdi_fhat(brts, pars, c(0L, 0L, 0L), N = N, link = link),
                  error = function(e) c(fhat = NA, sd_lw = NA, ess = NA, n = NA, n_nonfinite = NA))
    d1 <- DDD::bd_loglik(pars1 = c(th[1], th[2], 0, 0), pars2 = c(0, 0, 1, 0, 2), brts = brts, missnumspec = 0)
    c(lambda = th[1], mu = th[2], fhat = b[["fhat"]], sd_lw = b[["sd_lw"]], ess = b[["ess"]],
      ddd_phylo = d1, gap = b[["fhat"]] - d1)
  }))
  cat(sprintf("-- N = %d, link = %d (%s)\n", N, link, if (link == 0L) "linear" else "exponential"))
  print(round(out, 8))
}
cat(sprintf("lfactorial(n-1) = %.4f for n = %d\n", lfactorial(n - 1), n))

cat("\n=== R4: small tree, near-critical theta (lambda > mu, small gap) ===\n")
phy2  <- ape::rphylo(12, 1, 0.5)
brts2 <- sort(as.numeric(ape::branching.times(phy2)), decreasing = TRUE)
for (th in list(c(1, 0.5), c(0.5, 0.499), c(0.5, 0.4999), c(0.5, 0.49999))) {
  b <- tryCatch(bdi_fhat(brts2, th, c(0L, 0L, 0L), N = 100L),
                error = function(e) paste("ERROR:", conditionMessage(e)))
  d1 <- DDD::bd_loglik(pars1 = c(th[1], th[2], 0, 0), pars2 = c(0, 0, 1, 0, 2), brts = brts2, missnumspec = 0)
  cat(sprintf("theta=(%g,%g): BDI -> %s | DDD phylo %.6f\n", th[1], th[2],
              if (is.character(b)) b else sprintf("fhat %.6f sd_lw %.1e ess %.1f gap %.2e", b[["fhat"]], b[["sd_lw"]], b[["ess"]], b[["fhat"]] - d1), d1))
}

cat("\n=== R5: dd on a second tree: BDI vs thinning vs DDD::dd_loglik - lfactorial(n-1) ===\n")
thin_fhat <- function(brts, pars, model_bin, N = 300L, rho = 1) {
  p8 <- ns$.expand_pars(pars, model_bin)
  r <- ns$augment_trees(brts, p8, sample_size = N, maxN = 20L * N, max_missing = 1e4, max_lambda = 500,
                        num_threads = 1L, model = as.integer(model_bin), link = 0L, rho = rho)
  lw <- r$logf - r$logg; lw <- lw[is.finite(lw)]; m <- max(lw)
  fh <- log(sum(exp(lw - m)) / (length(lw) + r$rejected_zero_weights)) + m
  c(fhat = fh, ess = ess(lw), n = length(lw))
}
dd_true <- c(0.7, -0.015, 0.15, 0)
sim <- NULL
for (try in 1:30) { s <- tryCatch(simulate_tree(pars = dd_true, max_t = 5, model = "dd", max_lin = 200), error = function(e) NULL)
  if (!is.null(s) && inherits(s$tes, "phylo") && Ntip(s$tes) >= 10 && Ntip(s$tes) <= 40) { sim <- s; break } }
if (is.null(sim)) cat("no dd tree in range; skipping R5\n") else {
  brts_dd <- ns$.extract_brts(sim); n_dd <- length(brts_dd) + 1
  cat(sprintf("dd tree: %d tips, crown %.3f\n", n_dd, brts_dd[1]))
  for (th in list(c(0.7, -0.015, 0.15, 0), c(0.9, -0.02, 0.1, 0))) {
    K <- (th[1] - th[3]) / (-th[2])
    b  <- t(replicate(3, bdi_fhat(brts_dd, th, c(1L, 0L, 0L), N = 200L)))
    tt <- t(replicate(3, thin_fhat(brts_dd, th, c(1L, 0L, 0L), N = 300L)))
    ddd <- DDD::dd_loglik(pars1 = c(th[1], th[3], K), pars2 = c(200, 1, 0, 0, 0, 2), brts = brts_dd, missnumspec = 0)
    cat(sprintf("theta=(%g,%g,%g,0) K=%.1f | BDI %s (ess %s, nonfinite %s) | thin %s (ess %s) | DDD-lfact(n-1) %.3f\n",
                th[1], th[2], th[3], K, paste(round(b[, "fhat"], 3), collapse = "/"), paste(round(b[, "ess"]), collapse = "/"),
                paste(b[, "n_nonfinite"], collapse = "/"), paste(round(tt[, "fhat"], 3), collapse = "/"),
                paste(round(tt[, "ess"]), collapse = "/"), ddd - lfactorial(n_dd - 1)))
    cat(sprintf("   BDI - (DDD-lfact) = %s ; thin - (DDD-lfact) = %s ; mean BDI - mean thin = %.3f\n",
                paste(round(b[, "fhat"] - ddd + lfactorial(n_dd - 1), 3), collapse = "/"),
                paste(round(tt[, "fhat"] - ddd + lfactorial(n_dd - 1), 3), collapse = "/"),
                mean(b[, "fhat"]) - mean(tt[, "fhat"])))
  }
}
cat("\nDONE\n")

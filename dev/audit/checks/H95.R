# H95: "No test touches any .bdi_* function though BDI is the default proposal
#       and README claims exactness."
# Part A: static — grep the shipped tests.
# Part B: dynamic — run the shipped test suite against the installed build with
#         call counters traced onto every .bdi_* / .augment_tree_bdi / .mcem_bdi.
# Part C: the four pins the hypothesis proposes, each against a reference:
#   (i)   CR, lambda > mu: sd(lw) and ESS == N
#   (ii)  fhat == DDD::bd_loglik(cond = 0) + const across a theta grid
#   (iii) dd: BDI fhat vs thinning fhat vs DDD::dd_loglik(ddmodel = 1, cond = 0)
#   (iv)  mu > lambda, and rho < 1
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressMessages({ library(emphasis); library(ape); library(DDD); library(testthat) })
ns  <- asNamespace("emphasis")
pkg <- "/Users/pancho/Code/emphasis"
set.seed(95)

cat("=================== PART A: static grep of tests/ ===================\n")
tfiles <- list.files(file.path(pkg, "tests"), pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
hits <- unlist(lapply(tfiles, function(f) { l <- readLines(f); w <- grep("bdi", l, ignore.case = TRUE); if (length(w)) paste0(basename(f), ":", w, ": ", l[w]) }))
cat("test files:", length(tfiles), "  lines mentioning 'bdi':", length(hits), "\n"); if (length(hits)) cat(hits, sep = "\n")
skips <- unlist(lapply(tfiles, function(f) { l <- readLines(f); w <- grep("^\\s*skip\\(", l); if (length(w)) paste0(basename(f), ":", w) }))
cat("unconditional skip() calls:", length(skips), "\n"); cat(skips, sep = "  "); cat("\n")

cat("\n=================== PART B: call counters during the shipped test suite ===================\n")
bdi_funs <- c(grep("^\\.bdi_", ls(ns, all.names = TRUE), value = TRUE), ".augment_tree_bdi", ".mcem_bdi")
.H95_counts <- new.env(); for (f in bdi_funs) assign(f, 0L, envir = .H95_counts)
for (f in bdi_funs) {
  suppressMessages(trace(f, where = ns, print = FALSE,
        tracer = substitute({ assign(FN, get(FN, envir = .H95_counts) + 1L, envir = .H95_counts) }, list(FN = f))))
}
# also count the two dispatchers that would route to BDI
for (f in c(".run_mcem", ".mcem_dynamic_fresh")) { assign(f, 0L, envir = .H95_counts)
  suppressMessages(trace(f, where = ns, print = FALSE,
        tracer = substitute({ assign(FN, get(FN, envir = .H95_counts) + 1L, envir = .H95_counts) }, list(FN = f)))) }
res <- tryCatch(
  testthat::test_dir(file.path(pkg, "tests/testthat"), package = "emphasis", load_package = "installed",
                     reporter = "silent", stop_on_failure = FALSE, stop_on_warning = FALSE),
  error = function(e) { cat("test_dir error:", conditionMessage(e), "\n"); NULL })
if (!is.null(res)) { df <- as.data.frame(res)
  cat(sprintf("test suite: %d tests, %d passed-expectations, %d failed, %d skipped\n",
              nrow(df), sum(df$passed), sum(df$failed), sum(df$skipped))) }
for (f in c(bdi_funs, ".run_mcem", ".mcem_dynamic_fresh")) suppressMessages(untrace(f, where = ns))
cnt <- sapply(c(bdi_funs, ".run_mcem", ".mcem_dynamic_fresh"), function(f) get(f, envir = .H95_counts))
print(cnt)
cat("total calls into any .bdi_* / .augment_tree_bdi / .mcem_bdi during the test suite:", sum(cnt[bdi_funs]), "\n")

cat("\n=================== PART C: the proposed pins ===================\n")
ess <- function(lw) { w <- exp(lw - max(lw)); sum(w)^2 / sum(w^2) }
bdi_fhat <- function(brts, pars, model_bin, N = 200L, rho = 1) {
  r <- ns$.augment_tree_bdi(brts, pars, model_bin = model_bin, sample_size = N, link = 0L, rho = rho)
  c(fhat = r$fhat, sd_lw = sd(r$weights), ess = ess(r$weights), n = length(r$weights),
    n_nonfinite = sum(!is.finite(r$weights)))
}
thin_fhat <- function(brts, pars, model_bin, N = 300L, rho = 1) {
  p8 <- ns$.expand_pars(pars, model_bin)
  r <- ns$augment_trees(brts, p8, sample_size = N, maxN = 20L * N, max_missing = 1e4, max_lambda = 500,
                     num_threads = 1L, model = as.integer(model_bin), link = 0L, rho = rho)
  lw <- r$logf - r$logg; lw <- lw[is.finite(lw)]; m <- max(lw)
  # same convention as src/E_step.cpp:143-144: denominator = N + rejected_zero_weights
  fh <- log(sum(exp(lw - m)) / (length(lw) + r$rejected_zero_weights)) + m
  c(fhat = fh, sd_lw = sd(lw), ess = ess(lw), n = length(lw), rej0 = r$rejected_zero_weights)
}

## ---- a fixed CR tree (reconstructed, 20 tips) ----
phy <- ape::rphylo(20, birth = 0.5, death = 0.1)
brts <- sort(as.numeric(ape::branching.times(phy)), decreasing = TRUE)
cat(sprintf("CR tree: %d tips, crown age %.3f\n", Ntip(phy), brts[1]))

cat("\n--- (i) CR, lambda > mu: sd(lw), ESS (N = 200, 3 replicates) ---\n")
for (k in 1:3) print(round(bdi_fhat(brts, c(0.5, 0.1), c(0L, 0L, 0L)), 10))

cat("\n--- (ii) fhat vs DDD::bd_loglik(cond = 0, crown, soc = 2) across a theta grid ---\n")
grid <- rbind(c(0.5, 0.1), c(0.8, 0.1), c(0.3, 0.05), c(0.5, 0.4), c(0.5, 0.49), c(1.2, 0.9), c(0.5, 0.0))
out <- t(apply(grid, 1, function(th) {
  b <- bdi_fhat(brts, th, c(0L, 0L, 0L))
  d0 <- DDD::bd_loglik(pars1 = c(th[1], th[2], 0, 0), pars2 = c(0, 0, 0, 0, 2), brts = brts, missnumspec = 0)  # branching times
  d1 <- DDD::bd_loglik(pars1 = c(th[1], th[2], 0, 0), pars2 = c(0, 0, 1, 0, 2), brts = brts, missnumspec = 0)  # phylogeny
  c(lambda = th[1], mu = th[2], fhat = b[["fhat"]], sd_lw = b[["sd_lw"]], ess = b[["ess"]],
    ddd_brts = d0, ddd_phylo = d1, gap_brts = b[["fhat"]] - d0, gap_phylo = b[["fhat"]] - d1)
}))
print(round(out, 8))
cat(sprintf("range of (fhat - ddd_brts) over grid: %.3e ; range of (fhat - ddd_phylo): %.3e\n",
            diff(range(out[, "gap_brts"])), diff(range(out[, "gap_phylo"]))))
n <- length(brts) + 1
cat(sprintf("lfactorial(n-1) = %.6f  (n-1)*log(2) = %.6f  (labelled/unlabelled constants for n = %d)\n",
            lfactorial(n - 1), (n - 1) * log(2), n))

cat("\n--- (iii) dd (linear, gamma_N = 0): BDI vs thinning vs DDD::dd_loglik(ddmodel = 1, cond = 0) ---\n")
# generate a dd tree with emphasis's own simulator, keep it modest
dd_true <- c(0.8, -0.02, 0.1, 0)      # lambda(N) = 0.8 - 0.02 N  -> K = (0.8-0.1)/0.02 = 35
sim <- NULL
for (try in 1:20) { s <- tryCatch(simulate_tree(pars = dd_true, max_t = 6, model = "dd", max_lin = 200), error = function(e) NULL)
  if (!is.null(s) && !is.null(s$tes) && inherits(s$tes, "phylo") && Ntip(s$tes) >= 10 && Ntip(s$tes) <= 40) { sim <- s; break } }
if (is.null(sim)) { cat("could not draw a dd tree with 10-40 tips in 20 tries; skipping (iii)\n") } else {
  brts_dd <- ns$.extract_brts(sim)
  cat(sprintf("dd tree: %d tips, crown age %.3f\n", length(brts_dd) + 1, brts_dd[1]))
  th_dd <- rbind(c(0.8, -0.02, 0.1, 0), c(0.6, -0.01, 0.1, 0), c(1.0, -0.03, 0.2, 0))
  for (j in seq_len(nrow(th_dd))) {
    th <- th_dd[j, ]; K <- (th[1] - th[3]) / (-th[2])
    b <- t(replicate(3, bdi_fhat(brts_dd, th, c(1L, 0L, 0L), N = 200L)))
    tt <- t(replicate(3, thin_fhat(brts_dd, th, c(1L, 0L, 0L), N = 300L)))
    ddd <- DDD::dd_loglik(pars1 = c(th[1], th[3], K), pars2 = c(200, 1, 0, 0, 0, 2), brts = brts_dd, missnumspec = 0)
    cat(sprintf("theta=(%.2f,%.3f,%.2f,0) K=%.1f | BDI fhat %s (ess %s, nonfinite %s) | thinning fhat %s (ess %s) | DDD %.4f\n",
                th[1], th[2], th[3], K,
                paste(round(b[, "fhat"], 3), collapse = "/"), paste(round(b[, "ess"], 1), collapse = "/"),
                paste(b[, "n_nonfinite"], collapse = "/"),
                paste(round(tt[, "fhat"], 3), collapse = "/"), paste(round(tt[, "ess"], 1), collapse = "/"), ddd))
    cat(sprintf("   BDI - DDD = %s ; thinning - DDD = %s ; BDI - thinning = %s\n",
                paste(round(b[, "fhat"] - ddd, 3), collapse = "/"), paste(round(tt[, "fhat"] - ddd, 3), collapse = "/"),
                paste(round(mean(b[, "fhat"]) - mean(tt[, "fhat"]), 3), collapse = "/")))
  }
}

cat("\n--- (iv-a) mu > lambda under CR ---\n")
for (th in list(c(0.3, 0.5), c(0.3, 0.3), c(0.3, 0.31))) {
  r <- tryCatch(bdi_fhat(brts, th, c(0L, 0L, 0L), N = 50L), error = function(e) paste("ERROR:", conditionMessage(e)))
  d <- DDD::bd_loglik(pars1 = c(th[1], th[2], 0, 0), pars2 = c(0, 0, 0, 0, 2), brts = brts, missnumspec = 0)
  cat(sprintf("theta=(%.2f,%.2f): BDI -> %s | DDD bd_loglik = %.4f\n", th[1], th[2],
              if (is.character(r)) r else paste(names(r), round(r, 4), collapse = " "), d))
}

cat("\n--- (iv-b) rho < 1 under CR: BDI vs thinning fhat at the same theta ---\n")
for (rho in c(1, 0.7, 0.5)) {
  b <- t(replicate(3, bdi_fhat(brts, c(0.5, 0.1), c(0L, 0L, 0L), N = 200L, rho = rho)))
  tt <- t(replicate(3, thin_fhat(brts, c(0.5, 0.1), c(0L, 0L, 0L), N = 300L, rho = rho)))
  cat(sprintf("rho=%.1f: BDI fhat %s (sd_lw %s) | thinning fhat %s (ess %s) | BDI-thin mean gap %.3f | n_tips*log(rho) = %.3f\n",
              rho, paste(round(b[, "fhat"], 3), collapse = "/"), paste(signif(b[, "sd_lw"], 2), collapse = "/"),
              paste(round(tt[, "fhat"], 3), collapse = "/"), paste(round(tt[, "ess"], 1), collapse = "/"),
              mean(b[, "fhat"]) - mean(tt[, "fhat"]), Ntip(phy) * log(rho)))
}
cat("\nDONE\n")

## H77 — the M slot (model_bin[2]) : reachable? semantics in C++? affects cr/dd?
## Self-contained; ~1 min.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressMessages(library(emphasis))
options(warn = 1)
set.seed(7)

cat("=== A. Reachability of the M slot through the public API ===\n")
tr <- ape::drop.fossil(ape::rlineage(0.5, 0.0, 4))
while (ape::Ntip(tr) < 8 || ape::Ntip(tr) > 25) tr <- ape::drop.fossil(ape::rlineage(0.5, 0.0, 4))
brts <- sort(ape::branching.times(tr), decreasing = TRUE)
cat("ntips =", ape::Ntip(tr), "\n")

cat("resolve c(1,1,0): ", paste(emphasis:::.resolve_model(c(1L, 1L, 0L)), collapse = ","), "\n")
cat("resolve c(0,1,0): ", paste(emphasis:::.resolve_model(c(0L, 1L, 0L)), collapse = ","), "\n")
cat("formula ~ N + M : ",
    tryCatch({emphasis:::.resolve_model(~ N + M); "ACCEPTED"},
             error = function(e) paste("error:", conditionMessage(e))), "\n")
cat("par_names c(1,1,0): ", paste(emphasis:::.par_names(c(1L, 1L, 0L)), collapse = " "), "\n")
cat("model_label c(1,1,0): ", emphasis:::.model_label(c(1L, 1L, 0L)), "\n")
cat("bdi_supported c(1,1,0), link 0: ", emphasis:::.bdi_supported(c(1L, 1L, 0L), 0L), "\n")

fit_M <- tryCatch(
  estimate_rates(tr, method = "mcem", model = c(1L, 1L, 0L), link = "exponential",
                 control = list(lower_bound = c(-3, -0.5, -1, -5, -0.5, -1),
                                upper_bound = c( 1,  0.5,  1,  0,  0.5,  1),
                                sample_size = 10, max_iter = 3, burnin = 1,
                                num_threads = 1, verbose = FALSE)),
  error = function(e) e)
if (inherits(fit_M, "error")) {
  cat("estimate_rates(model = c(1,1,0)) ERROR:", conditionMessage(fit_M), "\n")
} else {
  cat("estimate_rates(model = c(1,1,0)) RAN. pars:\n"); print(round(fit_M$pars, 4))
  cat("class:", class(fit_M)[1], " loglik:", fit_M$loglik, "\n")
}
fit_M2 <- tryCatch(
  estimate_rates(tr, method = "mcem", model = c(1L, 1L, 0L), link = "exponential",
                 control = list(lower_bound = c(-3, -0.5, -1, -5, -0.5, -1),
                                upper_bound = c( 1,  0.5,  1,  0,  0.5,  1),
                                sample_size = 10, max_iter = 2, burnin = 1,
                                num_threads = 1, verbose = TRUE, sampling = "bdi")),
  error = function(e) e)
cat("with sampling='bdi' + M-active: ",
    if (inherits(fit_M2, "error")) paste("ERROR", conditionMessage(fit_M2)) else "ran (see message above)", "\n")

sim_M <- tryCatch(simulate_tree(pars = c(0.6, 0.0, 0.05, 0.05, 0.0, 0.0), model = c(1L, 1L, 0L),
                                max_t = 3, link = "linear"),
                  error = function(e) e)
cat("simulate_tree(model = c(1,1,0)): ",
    if (inherits(sim_M, "error")) paste("ERROR", conditionMessage(sim_M)) else
      paste("ran, ntips =", if (!is.null(sim_M$tes)) ape::Ntip(sim_M$tes) else NA), "\n")

cat("\n=== B. C++ semantics of an M-active model ===\n")
## B1: C++ ignores model_bin[2]; pars8[3] (beta_M) is always in the predictor.
p8_noM <- c(0.5, -0.01, 0.00, 0, 0.10, 0, 0.00, 0)
p8_M   <- c(0.5, -0.01, 0.05, 0, 0.10, 0, 0.02, 0)
aug <- emphasis:::augment_trees(brts, p8_M, 20L, 400L, 1e4, 500, 1L, c(1L, 1L, 0L), 0L, 1.0)
cat("augmented trees:", length(aug$trees), " with missing lineages: ",
    sum(vapply(aug$trees, function(d) any(d$t_ext == 0), logical(1))), "\n")
lf_100 <- emphasis:::eval_logf(p8_M, aug$trees, c(1L, 0L, 0L), 0L, 1.0)$logf
lf_110 <- emphasis:::eval_logf(p8_M, aug$trees, c(1L, 1L, 0L), 0L, 1.0)$logf
lf_010 <- emphasis:::eval_logf(p8_M, aug$trees, c(0L, 1L, 0L), 0L, 1.0)$logf
cat("max |logf(model=100) - logf(model=110)| at beta_M=0.05: ", max(abs(lf_100 - lf_110)), "\n")
cat("max |logf(model=010) - logf(model=110)| at beta_M=0.05: ", max(abs(lf_010 - lf_110)), "\n")
lf_noM <- emphasis:::eval_logf(p8_noM, aug$trees, c(1L, 1L, 0L), 0L, 1.0)$logf
cat("max |logf(beta_M=0.05) - logf(beta_M=0)| (same trees): ", max(abs(lf_110 - lf_noM)),
    "  -> slot is live, not inert\n")

## B2: hand recomputation of loglik with beta_M != 0, linear link, use_D = 0:
##   lambda_i = max(0, b0 + bN n_i + bM pd_i/n_i)  (M evaluated on the node ending the segment)
##   mu_i     = max(0, g0 + gN n_i + gM pd_i/n_i)
##   loglik   = -sum_i dt_i n_i (lambda_i + mu_i) + sum_{spec nodes, i<last} log lambda_i
##              + sum_{ext nodes} log mu_i
hand_loglik <- function(d, p) {
  n <- d$n; pd <- d$pd; M <- ifelse(n > 0, pd / n, 0)
  lam <- pmax(0, p[1] + p[2] * n + p[3] * M)
  mu  <- pmax(0, p[5] + p[6] * n + p[7] * M)
  dt  <- diff(c(0, d$brts))
  ext <- d$t_ext == 0
  k   <- nrow(d)
  -sum(dt * n * (lam + mu)) +
    sum(log(lam[!ext & seq_len(k) != k])) +
    sum(log(pmax(mu[ext], 1e-300)))
}
hand <- vapply(aug$trees, hand_loglik, numeric(1), p = p8_M)
cat("max |eval_logf - hand| on", length(hand), "trees (beta_M=.05, gamma_M=.02): ",
    max(abs(hand - lf_110)), "\n")
## and exponential link
lf_exp <- emphasis:::eval_logf(p8_M, aug$trees, c(1L, 1L, 0L), 1L, 1.0)$logf
hand_exp <- vapply(aug$trees, function(d, p) {
  n <- d$n; M <- ifelse(n > 0, d$pd / n, 0)
  lam <- exp(p[1] + p[2] * n + p[3] * M); mu <- exp(p[5] + p[6] * n + p[7] * M)
  dt <- diff(c(0, d$brts)); ext <- d$t_ext == 0; k <- nrow(d)
  -sum(dt * n * (lam + mu)) + sum(log(lam[!ext & seq_len(k) != k])) + sum(log(mu[ext]))
}, numeric(1), p = p8_M)
cat("max |eval_logf - hand| exponential link: ", max(abs(hand_exp - lf_exp)), "\n")

## B3: what M actually is on observed nodes (tip_start = 0 convention, cf. H9)
d1 <- aug$trees[[1]]
obs <- d1[d1$t_ext > 1e10, ]
cat("observed nodes: pd/(rank*brts) = ", paste(round(obs$pd / (seq_len(nrow(obs)) * obs$brts), 6), collapse = " "),
    "\n  -> M_obs(t_i) = i*t_i/n_i, not the true mean pendant age (H9 territory)\n")

cat("\n=== C. Does use_M = 0 leave cr / dd fits untouched? ===\n")
cat(".expand_pars(c(b0,g0), cr)  = ", paste(emphasis:::.expand_pars(c(0.5, 0.1), c(0L, 0L, 0L)), collapse = " "), "\n")
cat(".expand_pars(dd compact)    = ", paste(emphasis:::.expand_pars(c(0.5, -0.01, 0.1, 0.0), c(1L, 0L, 0L)), collapse = " "), "\n")
cat("bounds for dd -> lb8/ub8 slots 3,7 : ",
    paste(emphasis:::.expand_pars(c(0, -0.1, 0, -0.1), c(1L, 0L, 0L))[c(3, 7)], collapse = ","), "/",
    paste(emphasis:::.expand_pars(c(2, 0.1, 1, 0.1), c(1L, 0L, 0L))[c(3, 7)], collapse = ","), "\n")
## logf of a dd parameter (beta_M = 0) is the same under model 100 and 110 on the same trees
aug_dd <- emphasis:::augment_trees(brts, p8_noM, 20L, 400L, 1e4, 500, 1L, c(1L, 0L, 0L), 0L, 1.0)
a1 <- emphasis:::eval_logf(p8_noM, aug_dd$trees, c(1L, 0L, 0L), 0L, 1.0)
a2 <- emphasis:::eval_logf(p8_noM, aug_dd$trees, c(1L, 1L, 0L), 0L, 1.0)
cat("dd trees, beta_M=0: max|logf 100 - 110| =", max(abs(a1$logf - a2$logf)),
    " max|logg| diff =", max(abs(a1$logg - a2$logg)), "\n")
## and cr against the closed form (rho = 1, cond = NULL): eval_logf on the bare tree at beta_M = 0
d0 <- data.frame(brts = rev(brts[1] - brts)[-1], n = 2:length(brts), t_ext = 1e11, pd = 0)
d0 <- rbind(d0, data.frame(brts = brts[1], n = length(brts) + 1, t_ext = 1e11, pd = 0))
d0$pd <- seq_len(nrow(d0)) * d0$brts
cr8 <- c(0.5, 0, 0, 0, 0.0, 0, 0, 0)
lf_cr <- emphasis:::eval_logf(cr8, list(d0), c(0L, 0L, 0L), 0L, 1.0)$logf
lf_cr_Mflag <- emphasis:::eval_logf(cr8, list(d0), c(0L, 1L, 0L), 0L, 1.0)$logf
## pure-birth closed form for the observed tree: sum_{i<last} log(lambda) - lambda * sum_i n_i dt_i
pb <- (nrow(d0) - 1) * log(0.5) - 0.5 * sum(d0$n * diff(c(0, d0$brts)))
cat("cr, mu=0, no missing: eval_logf =", lf_cr, " pure-birth closed form =", pb,
    " use_M flag =", lf_cr_Mflag, "\n")

cat("\n=== D. Documentation / test tally ===\n")
cat("README 'not a user covariate':", length(grep("not a user covariate", readLines("/Users/pancho/Code/emphasis/README.md"))), "\n")
cat("estimate_rates.Rd 'use_M is internal only and should be 0':",
    length(grep("use_M.*internal only", readLines("/Users/pancho/Code/emphasis/man/estimate_rates.Rd"))), "\n")
cat("emphasis_cem.Rd mentions use_M with no restriction:",
    length(grep("use_M", readLines("/Users/pancho/Code/emphasis/man/emphasis_cem.Rd"))), "\n")
cat("augment_trees roxygen model doc says use_P/use_E (stale):",
    length(grep("use_P, use_E", readLines("/Users/pancho/Code/emphasis/R/RcppExports.R"))), "\n")
tst <- list.files("/Users/pancho/Code/emphasis/tests/testthat", full.names = TRUE, pattern = "\\.R$")
hits <- unlist(lapply(tst, function(f) { l <- readLines(f); if (any(grepl("estimate_rates\\(.*c\\(1L?, 1L?", l) | grepl("model *= *c\\([01]L?, *1L?, *[01]L?\\)", l))) basename(f) }))
cat("test files fitting/using an M-active model:", if (length(hits)) paste(hits, collapse = ", ") else "none (only test-em.R:23, which is skip()ped)", "\n")
cat("done\n")

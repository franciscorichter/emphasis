## H77 replication — vary tree (with extinction), link (linear), method (cem),
## model c(0,1,0), silent BDI fallback, and test whether use_N is any different
## from use_M on the C++ side.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressMessages(library(emphasis))
options(warn = 1)
set.seed(11)
mk <- function() { repeat { tr <- ape::drop.fossil(ape::rlineage(0.7, 0.25, 4)); if (ape::Ntip(tr) >= 15 && ape::Ntip(tr) <= 35) return(tr) } }
tr <- mk(); brts <- sort(ape::branching.times(tr), decreasing = TRUE)
cat("ntips =", ape::Ntip(tr), " crown age =", round(brts[1], 3), "\n")

cat("\n=== 1. estimate_rates M-active, LINEAR link, mcem, verbose=FALSE, sampling=bdi (silent fallback?) ===\n")
out <- capture.output(msgs <- capture.output(type = "message",
  fit1 <- tryCatch(estimate_rates(tr, method = "mcem", model = c(1L, 1L, 0L), link = "linear",
                  control = list(lower_bound = c(0, -0.1, -0.5, 0, -0.1, -0.5),
                                 upper_bound = c(2,  0.1,  0.5, 1,  0.1,  0.5),
                                 sample_size = 10, max_iter = 3, burnin = 1,
                                 num_threads = 1, verbose = FALSE, sampling = "bdi")),
                  error = function(e) e)))
cat("stdout lines:", length(out), " message lines:", length(msgs), "\n")
if (inherits(fit1, "error")) cat("ERROR:", conditionMessage(fit1), "\n") else {
  cat("ran; class", class(fit1)[1], "; model stored:", paste(fit1$model, collapse = ","), "\n")
  print(round(fit1$pars, 4))
  cat("print() header:\n"); print(fit1)
}

cat("\n=== 2. model c(0,1,0) (M only), exponential, mcem ===\n")
fit2 <- tryCatch(estimate_rates(tr, method = "mcem", model = c(0L, 1L, 0L), link = "exponential",
                  control = list(lower_bound = c(-3, -1, -5, -1), upper_bound = c(1, 1, 0, 1),
                                 sample_size = 10, max_iter = 3, burnin = 1, num_threads = 1, verbose = FALSE)),
                 error = function(e) e)
if (inherits(fit2, "error")) cat("ERROR:", conditionMessage(fit2), "\n") else { cat("ran:\n"); print(round(fit2$pars, 4)); cat("label:", emphasis:::.model_label(fit2$model), "\n") }

cat("\n=== 3. method = 'cem' with c(1,1,0) ===\n")
fit3 <- tryCatch(estimate_rates(tr, method = "cem", model = c(1L, 1L, 0L), link = "exponential",
                  control = list(lower_bound = c(-3, -0.5, -1, -5, -0.5, -1),
                                 upper_bound = c( 1,  0.5,  1,  0,  0.5,  1),
                                 sample_size = 10, max_iter = 2, num_threads = 1, verbose = FALSE)),
                 error = function(e) e)
if (inherits(fit3, "error")) cat("ERROR:", conditionMessage(fit3), "\n") else { cat("ran:\n"); print(round(fit3$pars, 4)) }

cat("\n=== 4. simulate_tree: c(0,1,0) and formula ~ N + M ===\n")
s1 <- tryCatch(simulate_tree(pars = c(0.6, 0.1, 0.05, 0.0), model = c(0L, 1L, 0L), max_t = 3, link = "linear"), error = function(e) e)
cat("simulate_tree c(0,1,0):", if (inherits(s1, "error")) conditionMessage(s1) else paste("ran, ntips", ape::Ntip(s1$tes)), "\n")
s2 <- tryCatch(simulate_tree(pars = c(0.6, 0.0, 0.05, 0.05, 0.0, 0.0), model = ~ N + M, max_t = 3), error = function(e) e)
cat("simulate_tree ~ N + M:", if (inherits(s2, "error")) conditionMessage(s2) else "ACCEPTED", "\n")
e2 <- tryCatch(estimate_rates(tr, method = "mcem", model = ~ N + M, link = "exponential",
                 control = list(lower_bound = rep(-1, 6), upper_bound = rep(1, 6), sample_size = 5, max_iter = 1)),
               error = function(e) e)
cat("estimate_rates ~ N + M:", if (inherits(e2, "error")) conditionMessage(e2) else "ACCEPTED", "\n")

cat("\n=== 5. Is use_N any different from use_M in C++? (only model_bin_[2] is read) ===\n")
p8 <- c(0.6, -0.02, 0.05, 0, 0.15, 0.01, 0.02, 0)
aug <- emphasis:::augment_trees(brts, p8, 15L, 400L, 1e4, 500, 1L, c(1L, 1L, 0L), 0L, 1.0)
trees <- aug$trees
lf <- function(mb, link = 0L, p = p8) emphasis:::eval_logf(p, trees, mb, link, 1.0)$logf
cat("beta_N=-0.02 live: max|logf(000) - logf(100)| =", max(abs(lf(c(0L,0L,0L)) - lf(c(1L,0L,0L)))),
    " (0 => C++ ignores use_N exactly as it ignores use_M)\n")
cat("max|logf(000) - logf(110)| =", max(abs(lf(c(0L,0L,0L)) - lf(c(1L,1L,0L)))), "\n")
cat("exp link: max|logf(000) - logf(110)| =", max(abs(lf(c(0L,0L,0L), 1L) - lf(c(1L,1L,0L), 1L))), "\n")
p8z <- p8; p8z[c(2, 6)] <- 0
cat("zeroing beta_N/gamma_N changes logf: max diff =", max(abs(lf(c(1L,1L,0L)) - lf(c(1L,1L,0L), p = p8z))), "\n")
cat("logg unaffected by flags (010 vs 000):", max(abs(emphasis:::eval_logf(p8, trees, c(0L,1L,0L), 0L, 1.0)$logg -
                                                     emphasis:::eval_logf(p8, trees, c(0L,0L,0L), 0L, 1.0)$logg)), "\n")

cat("\n=== 6. expand/contract for c(0,1,0) puts values in slots 3/7 ===\n")
cat(paste(emphasis:::.expand_pars(c(0.6, 0.05, 0.1, 0.02), c(0L, 1L, 0L)), collapse = " "), "\n")
cat("done\n")

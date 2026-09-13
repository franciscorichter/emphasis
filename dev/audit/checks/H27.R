## H27 — are the pipeline's stage logliks comparable, and does the
## "CEM if finite else GAM" init rule ignore their values?
##
## Test: run emphasis_pipeline() on a fixed 25-tip CR tree (3 replicates),
## then re-evaluate every stage's `pars` with ONE common estimator:
##   (a) 2000 thinning trees   (.simulate_particle + .is_fhat, C++)
##   (b) 300 BDI trees         (.augment_tree_bdi; zero-variance under CR, so exact)
##   (c) DDD::bd_loglik        (closed form, cond = 1 on crown survival)
## all with the same survival-GAM correction the pipeline applied.
## Compare with the logged per-stage loglik; record the init choice and
## whether it picked the point with the lower exact likelihood.
## Runs in ~2-4 min single-threaded.

.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressMessages({ library(emphasis); library(ape) })

set.seed(27)                            # reaches ape only; C++ RNG is clock-seeded
tree <- ape::rphylo(25, birth = 0.5, death = 0.2)
brts <- emphasis:::.extract_brts(tree)
model_bin <- c(0L, 0L, 0L); link <- 0L
cat(sprintf("tree: %d tips, crown age %.3f\n", Ntip(tree), brts[1]))

# is the installed pipeline the source's "cem first, break" rule?
body_txt <- paste(deparse(emphasis::emphasis_pipeline), collapse = "\n")
cat("installed pipeline has c(\"cem\",\"gam\") init loop:",
    grepl('c\\("cem", "gam"\\)', body_txt), "\n")

exact_ll <- function(p, cond = 1L) {
  # DDD: pars1 = (lambda, mu), pars2 = (tdmodel=0, cond, btorph=0, print=0, soc=2)
  DDD::bd_loglik(pars1 = c(p[1], p[2]), pars2 = c(0, cond, 0, 0, 2),
                 brts = brts, missnumspec = 0)
}

common_eval <- function(pars_compact, cond_fun, n = 2000L) {
  p8 <- emphasis:::.expand_pars(as.numeric(pars_compact), model_bin)
  thin <- emphasis:::.simulate_particle(brts, p8, model_bin, link,
                                        sample_size = n, maxN = 20L * n,
                                        max_missing = 1e4, max_lambda = 500,
                                        num_threads = 1L, rho = 1)
  f_thin <- if (is.null(thin)) NA_real_ else
    emphasis:::.is_fhat(thin$logf, thin$logg,
                        n_zero_weight = emphasis:::.n0(thin$rejected_zero_weights))
  ess_thin <- if (is.null(thin)) NA_real_ else
    emphasis:::.ess_from_lw(thin$logf - thin$logg)
  # BDI under CR is zero-variance IS (lw constant), so 300 trees is exact;
  # sd_lw_bdi is recorded to prove it.
  bdi <- emphasis:::.augment_tree_bdi(brts, p8, model_bin, sample_size = 300L,
                                      link = link, rho = 1)
  corr <- if (is.null(cond_fun)) 0 else cond_fun(p8)
  c(thin2000 = f_thin - corr, ess_thin = ess_thin,
    bdi2000 = bdi$fhat - corr, sd_lw_bdi = stats::sd(bdi$logf - bdi$logg),
    cond_corr = corr,
    exact_c1 = exact_ll(pars_compact, 1L), exact_c0 = exact_ll(pars_compact, 0L))
}

n_rep <- 3L
rows <- list(); init_rows <- list(); lag_rows <- list()
t_all <- proc.time()[3]
for (r in seq_len(n_rep)) {
  t0 <- proc.time()[3]
  pl <- emphasis_pipeline(tree, model = "cr", link = "linear",
                          control = list(num_threads = 1L, max_time = 150,
                                         mcem = list(max_iter = 25L)),
                          verbose = FALSE)
  cat(sprintf("\n[rep %d] pipeline done in %.0fs; best_stage=%s\n",
              r, proc.time()[3] - t0, pl$best_stage))
  print(pl$log[, c("stage", "status", "loglik", "AIC", "elapsed", "pars")])
  cond_fun <- if (!is.null(pl$bounds$survival_gam))
    emphasis:::.build_cond_fun(pl$bounds$survival_gam, model_bin) else NULL

  for (st in c("gam", "cem", "mcem")) {
    f <- pl$fits[[st]]
    if (is.null(f) || !is.finite(f$loglik)) next
    ce <- common_eval(f$pars, cond_fun)
    rows[[length(rows) + 1L]] <- data.frame(
      rep = r, stage = st, lambda = f$pars[1], mu = f$pars[2],
      logged = f$loglik, t(ce),
      gap_thin = f$loglik - ce[["thin2000"]],
      gap_bdi  = f$loglik - ce[["bdi2000"]])
  }

  # --- init rule: which stage was used, and was it the better point? -----
  ll_cem <- pl$fits$cem$loglik; ll_gam <- pl$fits$gam$loglik
  ex_cem <- exact_ll(pl$fits$cem$pars); ex_gam <- exact_ll(pl$fits$gam$pars)
  init_src <- if (is.finite(ll_cem)) "cem" else if (is.finite(ll_gam)) "gam" else "midpoint"
  init_rows[[r]] <- data.frame(
    rep = r, init_src = init_src, logged_cem = ll_cem, logged_gam = ll_gam,
    exact_cem = ex_cem, exact_gam = ex_gam,
    cem_worse_logged = ll_cem < ll_gam, cem_worse_exact = ex_cem < ex_gam)

  # --- MCEM lag: reported loglik = fhat(theta_{K-1}); pars = theta_K -----
  tr <- pl$mcem_trace
  if (!is.null(tr) && nrow(tr) >= 2L) {
    K <- nrow(tr)
    th_K1c <- emphasis:::.contract_pars(as.numeric(tr[K-1, grep("^par", names(tr))]), model_bin)
    th_Kc  <- emphasis:::.contract_pars(as.numeric(tr[K,   grep("^par", names(tr))]), model_bin)
    e_K1 <- common_eval(th_K1c, cond_fun, n = 1000L)
    e_K  <- common_eval(th_Kc,  cond_fun, n = 1000L)
    lag_rows[[r]] <- data.frame(
      rep = r, K = K, stop = pl$fits$mcem$details$stop_reason,
      reported = pl$fits$mcem$loglik,
      bdi_at_thetaKm1 = e_K1[["bdi2000"]], bdi_at_thetaK = e_K[["bdi2000"]],
      exact_Km1 = e_K1[["exact_c1"]], exact_K = e_K[["exact_c1"]],
      dpar_max = max(abs(th_Kc - th_K1c)))
  }
}
cat(sprintf("\nall replicates: %.0fs\n", proc.time()[3] - t_all))

res <- do.call(rbind, rows); rownames(res) <- NULL
cat("\n=== per-stage: logged loglik vs common re-evaluation of the SAME pars ===\n")
print(res, digits = 5)

cat("\n=== spread of logged vs common values across stages within a replicate ===\n")
for (r in seq_len(n_rep)) {
  s <- res[res$rep == r, ]
  cat(sprintf("rep %d: range(logged)=%.2f  range(thin2000)=%.2f  range(bdi2000)=%.2f  range(exact_c1)=%.2f\n",
              r, diff(range(s$logged)), diff(range(s$thin2000)),
              diff(range(s$bdi2000)), diff(range(s$exact_c1))))
  cat(sprintf("        stage ranking by logged : %s\n", paste(s$stage[order(-s$logged)], collapse = " > ")))
  cat(sprintf("        stage ranking by bdi2000: %s\n", paste(s$stage[order(-s$bdi2000)], collapse = " > ")))
}

cat("\n=== init rule ===\n")
print(do.call(rbind, init_rows), digits = 5)
cat("\n=== MCEM reported loglik lag (theta_{K-1} vs theta_K) ===\n")
if (length(lag_rows)) print(do.call(rbind, lag_rows), digits = 5)

# --- does the init choice change the MCEM endpoint? (rep-1 fits reused) ---
cat("\n=== MCEM from gam-init vs cem-init (same control), 2 runs each ===\n")
ab <- pl$bounds
run_from <- function(init) {
  f <- estimate_rates(tree, method = "mcem", model = "cr", init_pars = init,
                      cond = ab$survival_gam,
                      control = list(lower_bound = ab$lower_bound,
                                     upper_bound = ab$upper_bound,
                                     num_threads = 1L, max_iter = 25L,
                                     max_time = 100, sample_size = 200L))
  c(f$pars, loglik = f$loglik, iters = nrow(f$details$mcem),
    exact = exact_ll(f$pars))
}
out <- rbind(
  gam_init_1 = run_from(pl$fits$gam$pars), gam_init_2 = run_from(pl$fits$gam$pars),
  cem_init_1 = run_from(pl$fits$cem$pars), cem_init_2 = run_from(pl$fits$cem$pars))
print(out, digits = 5)

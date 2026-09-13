## H24 -- CEM adaptive-SD stopping rule cannot fire within the pipeline's max_iter = 20.
## Self-contained. Run: Rscript dev/audit/checks/H24.R
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({ library(emphasis); library(ape) })

cat("emphasis version:", as.character(packageVersion("emphasis")), "\n\n")

## ------------------------------------------------------------------
## Part A: analytic re-implementation of the schedule in R/de.R:696-737
## ------------------------------------------------------------------
sd_decay <- 0.85; sd_min_frac <- 0.01; patience <- 5L
n_decays_to_floor <- ceiling(log(sd_min_frac) / log(sd_decay))
cat(sprintf("A1  decays needed to reach floor: ceil(log(0.01)/log(0.85)) = %d  (0.85^%d = %.5f, 0.85^%d = %.5f)\n",
            n_decays_to_floor, n_decays_to_floor - 1, sd_decay^(n_decays_to_floor - 1),
            n_decays_to_floor, sd_decay^n_decays_to_floor))

## Simulate the loop for a given improving/non-improving sequence.
## Order inside iteration k (de.R): eval -> update plateau_count (k>1 only) ->
## stop check (plateau_count >= patience && sd at floor) -> resample -> decay if !improving.
## At k = 1 `improving` is FALSE (de.R:697) so a decay happens after iteration 1.
schedule <- function(improving_seq, max_iter) {
  sd <- 1; floor_ <- sd_min_frac; plateau <- 0L; stop_k <- NA_integer_
  sd_at_k <- numeric(max_iter)
  for (k in seq_len(max_iter)) {
    imp <- if (k > 1L) improving_seq[k] else FALSE
    if (k > 1L) plateau <- if (imp) 0L else plateau + 1L
    sd_at_k[k] <- sd                       # SD used to draw the population evaluated at k
    if (plateau >= patience && sd <= floor_ * (1 + 1e-8)) { stop_k <- k; break }
    if (!imp) sd <- max(sd * sd_decay, floor_)
  }
  list(stop_k = stop_k, sd_at_k = sd_at_k, sd_final = sd)
}
best_case <- schedule(rep(FALSE, 100), 100)          # never improving: fastest possible stop
cat(sprintf("A2  earliest possible plateau stop (never improving): iteration %d\n", best_case$stop_k))
s20 <- schedule(rep(FALSE, 20), 20)
cat(sprintf("A3  max_iter = 20, never improving: stop_k = %s, SD fraction at iteration 20 = %.4f (= 0.85^19)\n",
            as.character(s20$stop_k), s20$sd_at_k[20]))
s50 <- schedule(rep(FALSE, 50), 50)
cat(sprintf("A4  max_iter = 50 (estimate_rates default), never improving: stop_k = %d\n", s50$stop_k))
## how many improving iterations can max_iter = 50 tolerate and still stop by plateau?
tol_imp <- NA
for (m in 0:49) {
  seqv <- rep(FALSE, 50); if (m > 0) seqv[2:(m + 1)] <- TRUE
  if (is.na(schedule(seqv, 50)$stop_k)) { tol_imp <- m - 1; break }
}
cat(sprintf("A5  max_iter = 50 stops by plateau only if at most %d of iterations 2..49 improved (>= tol)\n", tol_imp))
stopifnot(is.na(s20$stop_k), best_case$stop_k == 30L)

## ------------------------------------------------------------------
## Part B: empirical -- real CEM runs on a small CR tree
## ------------------------------------------------------------------
set.seed(24)
tr <- NULL
while (is.null(tr)) {
  t0 <- ape::rlineage(0.4, 0.1, Tmax = 8)
  t1 <- tryCatch(ape::drop.fossil(t0), error = function(e) NULL)
  if (!is.null(t1) && Ntip(t1) >= 15 && Ntip(t1) <= 30) tr <- t1
}
cat(sprintf("\nTree: %d tips, crown age %.2f\n", Ntip(tr), max(branching.times(tr))))
lb <- c(0, 0); ub <- c(1, 0.5)   # tight box: keeps augmentation cheap

run_cem <- function(max_iter, reps, num_particles = 20, num_trees = 1, max_time = 120) {
  out <- lapply(seq_len(reps), function(r) {
    tm <- system.time(
      fit <- estimate_rates(tr, method = "cem", model = "cr",
                            control = list(max_iter = max_iter, num_particles = num_particles,
                                           num_trees = num_trees, num_threads = 1L,
                                           lower_bound = lb, upper_bound = ub,
                                           max_time = max_time, verbose = FALSE)))[3]
    bl <- fit$details$best_loglik
    d  <- diff(bl)
    improving <- c(FALSE, is.finite(d) & d >= 1e-4)      # de.R:697-701 rule, k=1 never improving
    n_dec <- min(sum(!improving[seq_len(length(bl) - 1L)]), n_decays_to_floor)  # decays before last eval
    data.frame(rep = r, max_iter = max_iter, converged = fit$details$converged,
               k_ran = length(bl), n_improving = sum(improving),
               sd_frac_at_last = 0.85^n_dec, pars1 = fit$pars[1], pars2 = fit$pars[2],
               loglik = fit$loglik, secs = round(tm, 1))
  })
  do.call(rbind, out)
}

cat("\nB1  estimate_rates(method='cem') with the pipeline's max_iter = 20:\n")
b1 <- run_cem(20, 2); print(b1, row.names = FALSE)
cat("\nB2  estimate_rates(method='cem') with the estimate_rates default max_iter = 50:\n")
b2 <- run_cem(50, 2); print(b2, row.names = FALSE)

## B3: the pipeline's own CEM stage (its defaults: max_iter 20, 50 particles, 5 trees)
cat("\nB3  emphasis_pipeline(stages = c('bounds','cem')), pipeline max_iter default (20), 20 particles x 1 tree:\n")
b3 <- lapply(1:1, function(r) {
  tm <- system.time(
    pp <- emphasis_pipeline(tr, model = "cr", stages = c("bounds", "cem"),
                            control = list(num_threads = 1L, max_time = 120,
                                           cem = list(num_particles = 20, num_trees = 1)),
                            verbose = FALSE))[3]
  bl <- pp$fits$cem$details$best_loglik
  data.frame(rep = r, converged = pp$fits$cem$details$converged, k_ran = length(bl),
             n_improving = sum(c(FALSE, diff(bl) >= 1e-4)), secs = round(tm, 1))
})
b3 <- do.call(rbind, b3); print(b3, row.names = FALSE)


## B4: can the rule fire at all? max_iter = 100 with 1 and with 5 trees per particle.
## The plateau needs 5 consecutive iterations with best_loglik improvement < 1e-4; with a
## noisy best-of-population fhat that rarely happens, whatever max_iter is.
cat("\nB4  max_iter = 100, num_trees = 1 (3 reps) and num_trees = 5 (2 reps):\n")
b4a <- run_cem(100, 3, num_trees = 1, max_time = 100)
b4b <- run_cem(100, 2, num_trees = 5, max_time = 100)
b4  <- rbind(cbind(num_trees = 1, b4a), cbind(num_trees = 5, b4b)); print(b4, row.names = FALSE)
run_len <- function(bl) { imp <- c(FALSE, diff(bl) >= 1e-4); r <- rle(!imp); max(r$lengths[r$values]) }

## Summary assertions
cat("\nSUMMARY\n")
cat(sprintf("  max_iter=20 runs stopping on 'max_iter': %d/%d\n", sum(b1$converged == "max_iter"), nrow(b1)))
cat(sprintf("  pipeline CEM stage stopping on 'max_iter': %d/%d\n", sum(b3$converged == "max_iter"), nrow(b3)))
cat(sprintf("  max_iter=50 runs: converged = %s ; k_ran = %s ; n_improving = %s\n",
            paste(b2$converged, collapse = ","), paste(b2$k_ran, collapse = ","),
            paste(b2$n_improving, collapse = ",")))
cat(sprintf("  max_iter=100 runs (trees=%s): converged = %s ; k_ran = %s ; n_improving = %s\n",
            paste(b4$num_trees, collapse = ","), paste(b4$converged, collapse = ","),
            paste(b4$k_ran, collapse = ","), paste(b4$n_improving, collapse = ",")))

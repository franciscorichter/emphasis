## H11 — BDI: no finite-weight filtering; log_sum sign at lambda = 0.
## Self-contained. Run: Rscript dev/audit/checks/H11.R
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages(library(emphasis))
options(width = 120)

cat("== Part A: eval_logf on hand-built trees (linear link, dd model_bin = c(1,0,0)) ==\n")
# Tree data frame in the BDI/C++ convention: forward times, crown at 0 (not a node),
# observed speciation nodes at brts 1..k with n = 2..k+1, closing sentinel at tp.
mk_tree <- function(k, tp = k + 1) {
  data.frame(brts = c(seq_len(k), tp), n = c(2:(k + 1), k + 2), t_ext = 1e11,
             pd = 0, tip_start = 0, id = c(seq_len(k) - 1L, -1L), parent_id = -1L)
}
lam <- function(p8, N) pmax(0, p8[1] + p8[2] * N)
mu  <- function(p8, N) pmax(0, p8[5] + p8[6] * N)
hand_loglik <- function(p8, tr) {
  # sum log lambda over non-closing nodes - sum dt * n * (lambda + mu), rates at END node
  dt <- diff(c(0, tr$brts))
  N  <- tr$n
  sum(log(lam(p8, N[-nrow(tr)]))) - sum(dt * N * (lam(p8, N) + mu(p8, N)))
}
ev <- function(p8, tr) emphasis:::eval_logf(p8, list(tr), model = c(1L, 0L, 0L), link = 0L, rho = 1)$logf
p8 <- function(b0, bN, g0 = 0.1) c(b0, bN, 0, 0, g0, 0, 0, 0)

# A0: finite case, checks that the hand formula and the C++ agree (tree builder is right)
tr <- mk_tree(3); p <- p8(0.5, -0.05)
cat(sprintf("A0 finite 3-node tree: eval_logf = %.6f, hand = %.6f\n", ev(p, tr), hand_loglik(p, tr)))

# A1: one speciation node, lambda(2) = 0  (prod_ = 0, sum_ = 0 -> signbit(0)=FALSE -> +Inf)
tr <- mk_tree(1); p <- p8(0.2, -0.1)
cat(sprintf("A1 single node, lambda(N=2)=%g at the LAST +=  : eval_logf = %s   (hypothesis asserts -Inf)\n",
            lam(p, 2), format(ev(p, tr))))

# A2: two nodes, lambda(2)=0 then lambda(3)=0.1  (zero NOT last -> sum_ = -Inf -> -Inf)
tr <- mk_tree(2); p <- p8(-0.2, 0.1)
cat(sprintf("A2 lambda(2)=%g (first), lambda(3)=%g (last)    : eval_logf = %s\n",
            lam(p, 2), lam(p, 3), format(ev(p, tr))))

# A3: two nodes, lambda(2)=0.1 then lambda(3)=0     (zero last, sum_ = 0 -> +Inf)
tr <- mk_tree(2); p <- p8(0.3, -0.1)
cat(sprintf("A3 lambda(2)=%g (first), lambda(3)=%g (last)    : eval_logf = %s\n",
            lam(p, 2), lam(p, 3), format(ev(p, tr))))

# A4: 25 nodes, lambda = 10 on nodes 1..24 (prod_ crosses 1e21 -> sum_ > 0), lambda(26)=0 last
tr <- mk_tree(25); p <- p8(10 + 26 * (10 / 24) , -(10 / 24))   # lambda(N)= 10+26*s - s*N: lambda(26)=10? no: fix below
p <- c(26 * 0.5, -0.5, 0, 0, 0.1, 0, 0, 0)                       # lambda(N) = 13 - 0.5 N: lambda(2)=12 ... lambda(25)=0.5, lambda(26)=0
cat(sprintf("A4 25 nodes, lambda(2..25) in [%g,%g], lambda(26)=%g last: eval_logf = %s   (prod_ crossed upper threshold, sum_>0)\n",
            lam(p, 25), lam(p, 2), lam(p, 26), format(ev(p, tr))))

# A5: 30 nodes, lambda small (prod_ crosses 1e-19 -> sum_ < 0), lambda(31)=0 last.
# Exact binary fractions (1/128) so that beta0 + betaN*N is exactly 0 even under FMA contraction.
tr <- mk_tree(30); p <- c(31 / 128, -1 / 128, 0, 0, 0.1, 0, 0, 0)  # lambda(N) = (31 - N)/128
cat(sprintf("A5 30 nodes, lambda(2..30) in [%g,%g], lambda(31)=%g last: eval_logf = %s   (prod_ crossed lower threshold, sum_<0)\n",
            lam(p, 30), lam(p, 2), lam(p, 31), format(ev(p, tr))))
# A5b: same as A5 but with the 0.01-based parameters: FMA contraction of beta0 + betaN*N leaves a ~1e-17 residue
tr <- mk_tree(25); p <- c(26 * 0.01, -0.01, 0, 0, 0.1, 0, 0, 0)
cat(sprintf("A5b (0.26 - 0.01 N, N=26): R lambda = %g, eval_logf = %s   (finite => C++ eta was not exactly 0; FMA residue, not H11)\n",
            lam(p, 26), format(ev(p, tr))))

# A6: zero at a NON-last speciation node. Sequence of += : obs n=2 (.3), obs n=3 (.2), missing n=4 (.1),
# missing n=5 (0), [ext, ext skipped], obs n=4 (.1) last, closing sentinel excluded.
p  <- p8(0.5, -0.1)
tr <- data.frame(brts = 1:8, n = c(2, 3, 4, 5, 6, 5, 4, 5), t_ext = c(1e11, 1e11, 5, 6, 0, 0, 1e11, 1e11),
                 pd = 0, tip_start = c(0, 0, 3, 4, 3, 4, 0, 0), id = c(0L, 1L, 3L, 4L, 3L, 4L, 2L, -1L),
                 parent_id = c(-1L, -1L, 1L, 1L, 1L, 1L, -1L, -1L))
cat(sprintf("A6 zero at 4th += (n=5, lambda=%g), last += has lambda=%g: eval_logf = %s
",
            lam(p, 5), lam(p, 4), format(ev(p, tr))))

cat("\n== Part B: BDI dd augmentation at pars (1.5, -0.12, 0.4, 0); lambda(N) = 0 for N >= 12.5 ==\n")
suppressPackageStartupMessages(library(DDD))
set.seed(11)
# ddmodel 1: lambda = la0 - (la0-mu0) N/K ; slope -0.12 => K = 1.1/0.12
sim <- DDD::dd_sim(pars = c(1.5, 0.4, 1.1 / 0.12), age = 6, ddmodel = 1)
brts <- sort(as.numeric(sim$brts), decreasing = TRUE)
cat(sprintf("tree: %d tips, crown age %.3f\n", length(brts) + 1L, brts[1]))
pars <- c(1.5, -0.12, 0.4, 0)
summ <- function(pars, S = 200L, max_missing = 30L, reps = 5L, label = "") {
  out <- t(sapply(seq_len(reps), function(r) {
    e <- emphasis:::.augment_tree_bdi(brts, pars, model_bin = c(1L, 0L, 0L),
                                     sample_size = S, max_missing = max_missing, link = 0L, rho = 1)
    lw <- e$weights; m <- max(lw); w <- exp(lw - m); w_norm <- w / sum(w) * length(w)
    maxN <- max(sapply(e$trees, function(t) max(t$n)))
    c(n_trees = length(e$trees), n_pinf = sum(e$logf == Inf), n_ninf = sum(e$logf == -Inf),
      n_finite = sum(is.finite(e$logf)), fhat = e$fhat, nan_wnorm = sum(is.nan(w_norm)),
      maxN = maxN, ess_finite = emphasis:::.ess_from_lw(lw))
  }))
  cat(label, "\n"); print(round(out, 3)); invisible(out)
}
b1 <- summ(pars, label = "B1: pars (1.5,-0.12,0.4,0), sample_size 200, 5 replicates")
# B2: milder DD where augmented N cannot reach the zero (lambda(N)=0 at N=40)
b2 <- summ(c(0.8, -0.02, 0.2, 0), label = "B2: pars (0.8,-0.02,0.2,0): zero at N=40, unreachable with max_missing=30")

cat("\n== Part C: what does m_cpp do with the NaN weights? ==\n")
e <- emphasis:::.augment_tree_bdi(brts, pars, model_bin = c(1L, 0L, 0L), sample_size = 200L,
                                 max_missing = 30L, link = 0L, rho = 1)
lw <- e$weights; w_norm <- exp(lw - max(lw)); w_norm <- w_norm / sum(w_norm) * length(w_norm)
cat(sprintf("logf: %d +Inf, %d -Inf, %d finite; fhat = %s; NaN in w_norm: %d/%d\n",
            sum(e$logf == Inf), sum(e$logf == -Inf), sum(is.finite(e$logf)), format(e$fhat),
            sum(is.nan(w_norm)), length(w_norm)))
lb <- c(0.01, -1, 0.001, 0); ub <- c(5, 0, 2, 0)
run_m <- function(w) {
  es <- list(trees = e$trees, weights = w, rejected = 0L, rejected_overruns = 0L, rejected_lambda = 0L,
             rejected_zero_weights = 0L, time = 0, fhat = 0)
  m <- emphasis:::m_cpp(e_step = es, init_pars = c(1.2, -0.1, 0.3, 0), plugin = "rpd1", lower_bound = lb,
             upper_bound = ub, xtol_rel = 1e-4, num_threads = 1L, model = c(1L, 0L, 0L), link = 0L, rho = 1)
  cat(sprintf("   [nlopt code %d]\n", m$nlopt))
  as.numeric(m$estimates)
}
# correct handling: +Inf logf is a zero-density tree -> weight 0
lw_fix <- ifelse(is.finite(e$logf), lw, -Inf)
w_fix <- exp(lw_fix - max(lw_fix)); w_fix <- w_fix / sum(w_fix) * length(w_fix)
cat("m_cpp estimates with the NaN w_norm actually produced by .mcem_bdi: ", paste(round(run_m(w_norm), 4), collapse = " "), "\n")
cat("m_cpp estimates with +Inf trees given weight 0 (correct handling):   ", paste(round(run_m(w_fix), 4), collapse = " "), "\n")
cat("m_cpp estimates with weights 1 (unweighted, for reference):          ", paste(round(run_m(rep(1, length(lw))), 4), collapse = " "), "\n")
fhat_fix <- log(mean(exp(lw_fix - max(lw_fix)))) + max(lw_fix)
cat(sprintf("fhat as computed (%s) vs fhat with +Inf treated as zero-density (%.4f)\n", format(e$fhat), fhat_fix))

cat("\n== Part D: estimate_rates(model='dd', sampling='bdi') short run from the generating pars ==\n")
fit <- estimate_rates(brts, method = "mcem", model = "dd", init_pars = pars,
                      control = list(lower_bound = lb, upper_bound = ub, sample_size = 200L, max_iter = 4L,
                                     max_missing = 30L, num_threads = 1L, sampling = "bdi", verbose = TRUE))
cat("names(fit):", paste(names(fit), collapse = " "), "\n")
traj <- if (!is.null(fit$mcem)) fit$mcem else if (!is.null(fit$details$mcem)) fit$details$mcem else NULL
if (!is.null(traj)) print(traj[, intersect(c(paste0("par", 1:4), "fhat", "num_trees"), names(traj))])
cat("stop_reason:", format(fit$details$stop_reason), " iterations:", format(fit$details$iterations), "\n")
cat("final pars:", paste(round(fit$pars, 4), collapse = " "), " loglik:", format(fit$loglik), "\n")
if (!is.null(fit$details$final_IS)) cat(sprintf("final_IS: fhat=%s ESS=%s n +Inf=%d n -Inf=%d\n",
   format(fit$details$final_IS$fhat), format(fit$details$final_IS$ESS),
   sum(fit$details$final_IS$logf == Inf), sum(fit$details$final_IS$logf == -Inf)))

cat("\n== Part E: why the M-step does not move — -Inf * 0 = NaN in the nlopt objective (M_step.cpp:49) ==\n")
# Same E-step draw `e` as Part C. Objective Q(theta) = sum_i loglik(theta, z_i) * w[i].
theta_cand <- c(1.4, -0.11, 0, 0, 0.35, 0, 0, 0)
lf_cand <- emphasis:::eval_logf(theta_cand, e$trees, model = c(1L, 0L, 0L), link = 0L, rho = 1)$logf
cat(sprintf("at candidate theta: %d trees have loglik=-Inf; sum(loglik * w_fix) = %s; sum over finite-logf trees only = %.4f\n",
            sum(lf_cand == -Inf), format(sum(lf_cand * w_fix)), sum((lf_cand * w_fix)[is.finite(lf_cand)])))
keep <- is.finite(e$logf)
cat(sprintf("max N among finite-logf trees = %d, among -Inf trees = %d (zero of lambda at theta_gen: N = 12.5)\n",
            max(sapply(e$trees[keep], function(t) max(t$n))), max(sapply(e$trees[!keep], function(t) max(t$n)))))
w_keep <- exp(lw[keep] - max(lw[keep])); w_keep <- w_keep / sum(w_keep) * length(w_keep)
run_m2 <- function(trees, w, init) {
  es <- list(trees = trees, weights = w, rejected = 0L, rejected_overruns = 0L, rejected_lambda = 0L,
             rejected_zero_weights = 0L, time = 0, fhat = 0)
  m <- emphasis:::m_cpp(e_step = es, init_pars = init, plugin = "rpd1", lower_bound = lb,
                        upper_bound = ub, xtol_rel = 1e-4, num_threads = 1L, model = c(1L, 0L, 0L), link = 0L, rho = 1)
  sprintf("%s  [nlopt code %d]", paste(round(as.numeric(m$estimates), 4), collapse = " "), m$nlopt)
}
cat("init = theta_gen (1.5,-0.12,0.4,0):\n")
cat("  all 200 trees, w_fix (zero weight on the ±Inf trees):", run_m2(e$trees, w_fix, pars), "\n")
cat("  109 finite-logf trees only, renormalised weights:   ", run_m2(e$trees[keep], w_keep, pars), "\n")

cat("\n== Part F: same fit with max_missing = 3 (augmented N <= 12, lambda(12) = 0.06 > 0: no zeros reachable) ==\n")
fit2 <- estimate_rates(brts, method = "mcem", model = "dd", init_pars = pars,
                       control = list(lower_bound = lb, upper_bound = ub, sample_size = 200L, max_iter = 4L,
                                      max_missing = 3L, num_threads = 1L, sampling = "bdi", verbose = FALSE))
print(fit2$details$mcem[, c("par1", "par2", "par5", "fhat", "delta_max")])
cat("stop_reason:", fit2$details$stop_reason, " final pars:", paste(round(fit2$pars, 4), collapse = " "), "\n")
if (!is.null(fit2$details$final_IS)) cat(sprintf("final_IS: n +Inf=%d n -Inf=%d\n",
   sum(fit2$details$final_IS$logf == Inf), sum(fit2$details$final_IS$logf == -Inf)))

cat("\n== Part G: fit at mild pars (0.8,-0.02,0.2,0), max_missing = 30: M-step moves when no tree hits lambda = 0 ==\n")
fit3 <- estimate_rates(brts, method = "mcem", model = "dd", init_pars = c(0.8, -0.02, 0.2, 0),
                       control = list(lower_bound = lb, upper_bound = ub, sample_size = 200L, max_iter = 4L,
                                      max_missing = 30L, num_threads = 1L, sampling = "bdi", verbose = FALSE))
print(fit3$details$mcem[, c("par1", "par2", "par5", "fhat", "delta_max")])
cat("stop_reason:", fit3$details$stop_reason, " final pars:", paste(round(fit3$pars, 4), collapse = " "), "\n")
if (!is.null(fit3$details$final_IS)) cat(sprintf("final_IS: n +Inf=%d n -Inf=%d\n",
   sum(fit3$details$final_IS$logf == Inf), sum(fit3$details$final_IS$logf == -Inf)))

cat("\n== Part H: 12 iterations from mild init (0.8,-0.02,0.2,0), max_missing = 30 — does an ordinary fit drift into the regime? ==\n")
t0 <- proc.time()[3]
fit4 <- estimate_rates(brts, method = "mcem", model = "dd", init_pars = c(0.8, -0.02, 0.2, 0),
                       control = list(lower_bound = lb, upper_bound = ub, sample_size = 200L, max_iter = 12L,
                                      max_missing = 30L, num_threads = 1L, sampling = "bdi", verbose = FALSE))
tr4 <- fit4$details$mcem
tr4$zero_at_N <- -tr4$par1 / tr4$par2
print(round(tr4[, c("par1", "par2", "par5", "fhat", "delta_max", "zero_at_N")], 4))
cat(sprintf("stop_reason: %s  iterations: %d  reported loglik: %s  (%.0fs)\n", fit4$details$stop_reason,
            nrow(tr4), format(fit4$loglik), proc.time()[3] - t0))
if (!is.null(fit4$details$final_IS)) cat(sprintf("final_IS: fhat=%s n +Inf=%d n -Inf=%d\n",
   format(fit4$details$final_IS$fhat), sum(fit4$details$final_IS$logf == Inf), sum(fit4$details$final_IS$logf == -Inf)))

cat("\n== Part I: at the Part-H final pars, does the M-step move iff no tree has logf = -Inf? (8 fresh E-steps) ==\n")
th <- as.numeric(fit4$details$pars)   # 8-vector
for (r in 1:8) {
  e5 <- emphasis:::.augment_tree_bdi(brts, th, model_bin = c(1L, 0L, 0L), sample_size = 200L, max_missing = 30L, link = 0L, rho = 1)
  lw5 <- e5$weights; w5 <- exp(lw5 - max(lw5)); w5 <- w5 / sum(w5) * length(w5)      # exactly bdi.R:796-799
  es <- list(trees = e5$trees, weights = w5, rejected = 0L, rejected_overruns = 0L, rejected_lambda = 0L,
             rejected_zero_weights = 0L, time = 0, fhat = 0)
  m <- emphasis:::m_cpp(e_step = es, init_pars = th, plugin = "rpd1", lower_bound = c(lb, rep(0, 4))[c(1,2,3,4,5,6,7,8)] * 0 + c(0.01, -1, 0, 0, 0.001, 0, 0, 0),
                        upper_bound = c(5, 0, 0, 0, 2, 0, 0, 0), xtol_rel = 1e-3, num_threads = 1L,
                        model = c(1L, 0L, 0L), link = 0L, rho = 1)
  est <- as.numeric(m$estimates)
  cat(sprintf("  rep %d: n(-Inf)=%3d n(+Inf)=%d  fhat=%9s  max|delta|=%.2e  nlopt=%d\n", r,
              sum(e5$logf == -Inf), sum(e5$logf == Inf), format(round(e5$fhat, 3)), max(abs(est - th)), m$nlopt))
}

args <- commandArgs(TRUE); lib <- args[1]; tag <- basename(lib)
.libPaths(c(lib, .libPaths())); suppressMessages(library(emphasis))
ev <- function(p8, tr, model = c(1L,0L,0L), link = 0L, rho = 1)
  emphasis:::eval_logf(p8, list(tr), model = model, link = link, rho = rho)
mk_tree <- function(k, tp = k + 1) data.frame(brts = c(seq_len(k), tp), n = c(2:(k + 1), k + 2),
  t_ext = 1e11, pd = 0, tip_start = 0, id = c(seq_len(k) - 1L, -1L), parent_id = -1L)
# augmented tree with missing + extinction nodes (A6 shape, k obs nodes + 2 missing + 2 ext)
mk_aug <- function() data.frame(brts = 1:8, n = c(2,3,4,5,6,5,4,5),
  t_ext = c(1e11,1e11,5,6,0,0,1e11,1e11), pd = c(0,0,1.5,2,0,0,3,3.5),
  tip_start = c(0,0,3,4,3,4,0,0), id = c(0L,1L,3L,4L,3L,4L,2L,-1L), parent_id = c(-1L,-1L,1L,1L,1L,1L,-1L,-1L))
out <- list()

## 1. Random finite battery across links x model_bin x trees -> compare bitwise across builds
set.seed(1)
mb_list <- list(cr = c(0L,0L,0L), dd = c(1L,0L,0L), pd = c(0L,1L,0L), ep = c(0L,0L,1L), all = c(1L,1L,1L))
trees <- list(t1 = mk_tree(1), t3 = mk_tree(3), t10 = mk_tree(10), aug = mk_aug(), t40 = mk_tree(40))
res <- list(); i <- 0
for (link in 0:2) for (mb in names(mb_list)) for (tn in names(trees)) for (r in 1:3) {
  p <- c(runif(1, 0.5, 3), rnorm(1, 0, 0.02), rnorm(1, 0, 0.05), rnorm(1, 0, 0.05),
         runif(1, 0.05, 0.5), rnorm(1, 0, 0.01), rnorm(1, 0, 0.02), rnorm(1, 0, 0.02))
  if (link == 0L) p[c(2,6)] <- -abs(p[c(2,6)]) * 0.1  # keep linear rates positive on N <= 42
  rho <- if (r == 3) 0.7 else 1
  e <- ev(p, trees[[tn]], mb_list[[mb]], link, rho)
  i <- i + 1; res[[i]] <- data.frame(link = link, mb = mb, tree = tn, r = r, rho = rho, logf = e$logf, logg = e$logg)
}
out$battery <- do.call(rbind, res)

## 2. Zero / negative / degenerate rates
Z <- list()
z <- function(name, val) Z[[name]] <<- val
p8 <- function(b0, bN, g0 = 0.1) c(b0, bN, 0, 0, g0, 0, 0, 0)
z("lin_zero_first",  ev(p8(-0.2, 0.1), mk_tree(2))$logf)
z("lin_zero_last",   ev(p8(0.3, -0.1), mk_tree(2))$logf)
z("lin_zero_only",   ev(p8(0.2, -0.1), mk_tree(1))$logf)
z("lin_zero_middle_of_10", ev(c(0.7, -0.1, 0,0, 0.1,0,0,0), mk_tree(10))$logf)   # lambda(7)=0 at node 6 of 10
z("lin_neg_all",     ev(p8(0.1, -0.1), mk_tree(2))$logf)
z("lin_beta0_zero_cr", ev(p8(0, 0), mk_tree(3), model = c(0L,0L,0L))$logf)          # cr with lambda = 0
z("lin_zero_upper_cross", ev(c(104, -4, 0,0, 0.1,0,0,0), mk_tree(25))$logf)
z("lin_zero_lower_cross", ev(c(31/128, -1/128, 0,0, 0.1,0,0,0), mk_tree(30))$logf)
z("lin_zero_aug_nonlast", ev(p8(0.5, -0.1), mk_aug())$logf)
z("lin_zero_aug_last",    ev(p8(0.4, -0.1), mk_aug())$logf)                      # lambda(4)=0 on the last obs node; lambda(5) < 0 too
z("lin_mu_zero_ext",  ev(c(0.5, -0.02, 0,0, 0, 0,0,0), mk_aug())$logf)          # mu = 0 on extinction nodes (clamped path)
z("lin_zero_closing_only", ev(c(0.7, -0.1, 0,0,0.1,0,0,0), mk_tree(5))$logf)   # lambda(7)=0 only on closing sentinel n=7 -> must stay finite
z("lin_rho_zero_rate", ev(p8(0.3, -0.1), mk_tree(2), rho = 0.5)$logf)
z("lin_zero_ep_model", ev(p8(0.3, -0.1), mk_tree(2), model = c(1L,0L,1L))$logf)
z("lin_zero_pd_model", ev(c(0.3, 0, -0.1, 0, 0.1,0,0,0), mk_aug(), model = c(0L,1L,0L))$logf) # beta_M negative: lambda hits 0 on pd/n>=3
# exponential link
z("exp_underflow_last", ev(c(1358.6, -690.8, 0,0, 0.1,0,0,0), mk_tree(2), link = 1L)$logf) # lambda(2)=1e-10, lambda(3)=exp(-713.8)~1e-310 (denormal, not zero)
z("exp_true_zero_last", ev(c(2000, -1000, 0,0, 0.1,0,0,0), mk_tree(2), link = 1L)$logf)   # lambda(2)=1, lambda(3)=exp(-1000)=0
z("exp_true_zero_first", ev(c(-3000, 1000, 0,0, 0.1,0,0,0), mk_tree(2), link = 1L)$logf)  # lambda(2)=0, lambda(3)=1
z("exp_prod_underflow_to_zero", ev(c(1400 + 2*725.3 - 2*(-23.03) , 0, 0,0,0.1,0,0,0)[1:8] * 0 + c(1381.57, -702.27, 0,0,0.1,0,0,0), mk_tree(2), link = 1L)$logf) # lambda(2)=1e-10, lambda(3)=1e-315 -> prod 1e-325 = 0
z("exp_overflow", ev(c(800, 0, 0,0, 0.1,0,0,0), mk_tree(2), link = 1L)$logf)
z("exp_big_tree_ok", ev(c(2, -0.05, 0,0, -1, 0,0,0), mk_tree(60), link = 1L)$logf)
# gaussian link
z("gau_intercept_zero", ev(c(0, 0.1, 0,0, 0.1,0,0,0), mk_tree(2), link = 2L)$logf)
z("gau_intercept_neg",  ev(c(-0.5, 0.1, 0,0, 0.1,0,0,0), mk_tree(2), link = 2L)$logf)
z("gau_ok",             ev(c(0.5, 0.1, 0,0, 0.1,0,0,0), mk_tree(5), link = 2L)$logf)
# NaN / Inf parameters
z("nan_beta0", ev(c(NaN, -0.1, 0,0, 0.1,0,0,0), mk_tree(2))$logf)
z("na_beta0",  ev(c(NA_real_, -0.1, 0,0, 0.1,0,0,0), mk_tree(2))$logf)
z("inf_beta0_lin", ev(c(Inf, -0.1, 0,0, 0.1,0,0,0), mk_tree(2))$logf)
z("neginf_beta0_lin", ev(c(-Inf, -0.1, 0,0, 0.1,0,0,0), mk_tree(2))$logf)
# tiny positive rates (finite log expected)
z("lin_tiny_rate_1e-300", ev(c(1e-300 + 2*1e-301, -1e-301, 0,0, 0.1,0,0,0), mk_tree(1))$logf)
z("lin_rate_denormal", ev(c(5e-324*3, -5e-324, 0,0,0.1,0,0,0), mk_tree(1))$logf)
# multi-tree
z("multi", emphasis:::eval_logf(p8(0.3, -0.1), list(mk_tree(1), mk_tree(2), mk_tree(3)), model = c(1L,0L,0L), link = 0L, rho = 1)$logf)
# empty-ish: single closing node (no speciation events)
z("closing_only", ev(p8(0.3, -0.1), data.frame(brts = 1, n = 2, t_ext = 1e11, pd = 0, tip_start = 0, id = -1L, parent_id = -1L))$logf)
out$zero <- Z
saveRDS(out, sprintf("/Users/pancho/Code/emphasis/dev/audit/review/1.1_adv_%s.rds", tag))
cat("done", tag, "\n")

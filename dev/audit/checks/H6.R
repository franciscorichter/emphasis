## H6 — thinning proposal density vs. `Model::sampling_prob` for M/D models.
##
## Claim under test: for models with beta_M/gamma_M or beta_D/gamma_D != 0 the
## quantity `logg` returned by the C++ (`sampling_prob`, model.hpp:255-328) is
## not the log density of the augmentation the thinning sampler
## (augment_tree.cpp:128-191 + model.hpp:217-252) actually drew, because
##   (a) extinction_time() draws lifetimes with extinction_rate() (no D) on the
##       proxy node whose pd is still 0 during augmentation (model.hpp:221),
##   (b) nh_rate() uses extinction_rate_ep / speciation_rate_ep on a copy of the
##       proxy node with pd recomputed at t (model.hpp:243-248),
##   (c) sampling_prob credits lambda/mu evaluated on the FINISHED node (final
##       pd, focal_tip_start) and a compensator that is not int nh(t) dt.
##
## Method: rebuild in R, from the returned augmented trees,
##   q_cpp   = an R transcription of sampling_prob        (must equal logg to 1e-8)
##   q_true  = the density of the thinning process as coded (nh(t) with the
##             state the sampler had at time t, truncated-exponential lifetime
##             with the mu the sampler drew with, uniform parent choice)
## For CR and DD (N-only) q_true == q_cpp up to the parent-choice term (H5),
## which validates the transcription. For M/D models the two differ; the
## difference is decomposed per node and its effect on fhat is measured.
## Everything is rho = 1 (the validation-study setting).

.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages(library(emphasis))
options(width = 120)

T_EXT_TIP <- 10e10; T_EXT_UNS <- 5e10
is_ext <- function(te) te == 0
is_tip <- function(te) te == T_EXT_TIP
is_uns <- function(te) te == T_EXT_UNS
is_mis <- function(te) !(is_ext(te) | is_tip(te) | is_uns(te))
lmexp  <- function(x) { m <- max(x); m + log(mean(exp(x - m))) }

## ---- rates, transcribed from model.hpp:143-212 ----------------------------
link_fn <- function(link, eta) if (link == 1L) exp(eta) else pmax(0, eta)
m_cov <- function(nd) if (nd$n > 0) nd$pd / nd$n else 0
e_s   <- function(nd) {
  if (is_ext(nd$t_ext)) return(nd$brts - nd$tip_start)
  if (nd$parent_id >= 0) return(nd$brts - nd$focal_tip_start)
  m_cov(nd)
}
lam_noD <- function(p, nd, link) link_fn(link, p[1] + p[2]*nd$n + p[3]*m_cov(nd))
mu_noD  <- function(p, nd, link) link_fn(link, p[5] + p[6]*nd$n + p[7]*m_cov(nd))
lam_ep  <- function(p, nd, link) { M <- m_cov(nd); D <- e_s(nd) - M
                                   link_fn(link, p[1] + p[2]*nd$n + p[3]*M + p[4]*D) }
mu_ep   <- function(p, nd, link) { M <- m_cov(nd); D <- e_s(nd) - M
                                   link_fn(link, p[5] + p[6]*nd$n + p[7]*M + p[8]*D) }
lam_of <- function(p, nd, link, useD) if (useD) lam_ep(p, nd, link) else lam_noD(p, nd, link)
mu_of  <- function(p, nd, link, useD) if (useD) mu_ep(p, nd, link)  else mu_noD(p, nd, link)
exp_integral <- function(a, b, t1, t2) if (abs(b) < 1e-12) exp(a)*(t2 - t1) else exp(a)*(exp(b*t2) - exp(b*t1))/b
row_node <- function(df, i) as.list(df[i, ])

## ---- q_cpp: R transcription of Model::sampling_prob (model.hpp:255-328) ----
q_cpp_R <- function(df, p, model_bin, link, rho = 1) {
  useD <- model_bin[3] == 1L; ep_exp <- useD && link == 1L
  TT <- df$brts[nrow(df)]; inte <- 0; logg <- 0; prev <- 0
  tips <- df$n[1]; Ne <- 0
  sbE <- if (ep_exp) df$n[1] else 0
  parent_terms <- numeric(0)
  for (i in seq_len(nrow(df))) {
    nd <- row_node(df, i)
    lambda <- lam_of(p, nd, link, useD); mu <- max(mu_of(p, nd, link, useD), 1e-10)
    dt <- nd$brts - prev
    int_segment <- dt - rho*(1/mu)*(exp(-mu*(TT - nd$brts)) - exp(-mu*(TT - prev)))
    if (ep_exp && dt > 0) {
      Nseg <- nd$n; Mseg <- if (Nseg > 0) nd$pd/Nseg else 0
      A_lam <- p[1] + p[2]*Nseg + (p[3] - p[4])*Mseg
      lam_integral <- sbE * exp_integral(A_lam, p[4], prev, nd$brts)
      inte <- inte + (lam_integral/dt) * int_segment
    } else {
      inte <- inte + nd$n * lambda * int_segment
    }
    tips <- tips + is_tip(nd$t_ext); Ne <- Ne - is_ext(nd$t_ext)
    if (is_mis(nd$t_ext)) {
      lifespan <- nd$t_ext - nd$brts
      parent_terms <- c(parent_terms, -log(2*tips + Ne))
      logg <- logg + log(nd$n*mu*lambda) - mu*lifespan - log(2*tips + Ne); Ne <- Ne + 1
    } else if (is_uns(nd$t_ext)) {
      parent_terms <- c(parent_terms, -log(2*tips + Ne))
      logg <- logg + log(nd$n*lambda*(1 - rho)) - mu*(TT - nd$brts) - log(2*tips + Ne); Ne <- Ne + 1
    }
    if (ep_exp) {
      if (is_mis(nd$t_ext) || is_uns(nd$t_ext)) sbE <- sbE + exp(-p[4]*nd$brts)
      else if (is_ext(nd$t_ext)) sbE <- sbE - exp(-p[4]*nd$tip_start)
      else if (i != nrow(df)) sbE <- sbE + exp(-p[4]*nd$brts)
    }
    prev <- nd$brts
  }
  list(logg = logg - inte, inte = inte, parent_terms = parent_terms)
}

## ---- q_true: density of the thinning process as coded ----------------------
## State at time t = final tree minus every augmented species with t_spec > t
## (insertions happen in time order, augment_tree.cpp:136-189). Fields as the
## sampler saw them: pd = 0 stored on every node (make_node / create_tree),
## tip_start = 0 for observed nodes and t_spec for augmented ones,
## focal_tip_start = 0 everywhere (compute_pendant_pd runs only at the end).
q_true_R <- function(df, p, model_bin, link, rho = 1) {
  stopifnot(rho == 1)
  useD <- model_bin[3] == 1L
  TT <- df$brts[nrow(df)]
  ext <- is_ext(df$t_ext); mis <- is_mis(df$t_ext); tip <- is_tip(df$t_ext)
  spec_start <- ifelse(ext, df$tip_start, ifelse(mis, df$brts, -Inf))   # when the row entered the tree
  ts_sampler <- ifelse(ext, NA, ifelse(mis, df$brts, 0))                 # tip_start used during augmentation
  sp <- !ext
  n_alive <- function(t) 2 + sum(sp & df$brts < t) - sum(ext & df$brts < t)
  pd_s    <- function(t) { k <- sp & df$brts <= t & df$t_ext > t; sum(t - ts_sampler[k]) }
  K_alive <- function(t) sum(sp & df$brts < t & df$t_ext > t)             # candidate parents (nodes only)
  proxy_i <- function(t) { st <- which(spec_start < t); cand <- st[df$brts[st] >= t]
                           if (length(cand)) min(cand) else max(st) }
  proxy_node <- function(t, pd) { i <- proxy_i(t); nd <- row_node(df, i)
                                  nd$n <- n_alive(t); nd$pd <- pd; nd$focal_tip_start <- 0
                                  if (mis[i]) nd$tip_start <- df$brts[i]; nd }
  nh <- function(t) {                                   # model.hpp:238-252
    cur <- proxy_node(t, pd_s(t))
    lambda <- lam_of(p, cur, link, useD); mu <- max(mu_of(p, cur, link, useD), 1e-10)
    lambda * cur$n * (1 - exp(-mu*(TT - t)))
  }
  ## compensator: integrate nh over pieces where the state is constant
  bp <- sort(unique(c(0, df$brts))); bp <- bp[bp <= TT]
  inte <- 0
  for (j in seq_len(length(bp) - 1)) if (bp[j+1] - bp[j] > 1e-12)
    inte <- inte + integrate(Vectorize(nh), bp[j], bp[j+1], rel.tol = 1e-10, subdivisions = 500L)$value
  ## event terms
  ev <- data.frame()
  logq <- -inte; logq_cpp_parent <- -inte
  for (i in which(mis)) {
    ts <- df$brts[i]; l <- df$t_ext[i] - ts; r <- TT - ts
    cur   <- proxy_node(ts, pd_s(ts))
    lam_nh <- lam_of(p, cur, link, useD); mu_nh <- max(mu_of(p, cur, link, useD), 1e-10)
    cur0  <- cur; cur0$pd <- 0                          # extinction_time(): raw node, pd still 0, no-D rate
    mu_d  <- mu_noD(p, cur0, link); if (mu_d <= 0) mu_d <- 1e-10
    K <- K_alive(ts)
    nh_ts  <- lam_nh * cur$n * (1 - exp(-mu_nh * r))
    ldtr   <- log(mu_d) - mu_d*l - log(1 - exp(-mu_d*r))
    term   <- log(nh_ts) + ldtr - (if (K > 0) log(K) else 0)
    logq <- logq + term
    nd <- row_node(df, i)
    ev <- rbind(ev, data.frame(t_spec = ts, life = l, n = cur$n, K = K,
                               lam_nh = lam_nh, lam_credit = lam_of(p, nd, link, useD),
                               mu_draw = mu_d, mu_nh = mu_nh,
                               mu_credit = max(mu_of(p, nd, link, useD), 1e-10),
                               pd_sampler = pd_s(ts), pd_final = nd$pd,
                               tip_start_final = nd$tip_start, parent_id = nd$parent_id))
  }
  list(logq = logq, inte = inte, ev = ev)
}

## ---- driver -----------------------------------------------------------------
brts <- c(5, 3.6, 2.9, 2.2, 1.5, 0.9, 0.4)        # 8-tip tree, crown age 5
draw <- function(p8, model_bin, link, N) {
  emphasis:::augment_trees(brts = brts, pars = p8, sample_size = as.integer(N),
                           maxN = 200L*N, max_missing = 10000L, max_lambda = 500,
                           num_threads = 1L, model = as.integer(model_bin),
                           link = as.integer(link), rho = 1)
}
analyse <- function(label, p8, model_bin, link, N = 150) {
  E <- draw(p8, model_bin, link, N)
  n_mis <- sapply(E$trees, function(d) sum(is_mis(d$t_ext)))
  d_cpp <- d_true <- d_true_cppparent <- d_inte <- numeric(N)
  ev_all <- data.frame()
  for (i in seq_len(N)) {
    d <- E$trees[[i]]
    qc <- q_cpp_R(d, p8, model_bin, link); qt <- q_true_R(d, p8, model_bin, link)
    d_cpp[i]  <- qc$logg - E$logg[i]                       # transcription check
    d_true[i] <- E$logg[i] - qt$logq                       # H6 (+H5 parent term)
    ## replace the sampler's -log K by the C++ -log(2 tips + Ne) to isolate H6 from H5
    parent_true <- if (nrow(qt$ev)) sum(ifelse(qt$ev$K > 0, -log(qt$ev$K), 0)) else 0
    d_true_cppparent[i] <- E$logg[i] - (qt$logq - parent_true + sum(qc$parent_terms))
    d_inte[i] <- qc$inte - qt$inte
    if (nrow(qt$ev)) ev_all <- rbind(ev_all, cbind(tree = i, qt$ev))
  }
  cat(sprintf("\n=== %s  pars8 = (%s)  model_bin = (%s) link = %d  N = %d, mean #missing = %.2f\n",
              label, paste(signif(p8, 3), collapse = ","), paste(model_bin, collapse = ","), link, N, mean(n_mis)))
  cat(sprintf("  transcription check  max|q_cpp_R - logg_cpp|      = %.2e\n", max(abs(d_cpp))))
  cat(sprintf("  logg_cpp - logq_true (parent term H5 included)   : mean %+.4f  sd %.4f  range [%+.4f, %+.4f]\n",
              mean(d_true), sd(d_true), min(d_true), max(d_true)))
  cat(sprintf("  logg_cpp - logq_true (C++ parent term, H6 only)  : mean %+.4f  sd %.4f  max|.| %.2e\n",
              mean(d_true_cppparent), sd(d_true_cppparent), max(abs(d_true_cppparent))))
  cat(sprintf("  compensator: inte_cpp - int nh dt                : mean %+.4f  sd %.4f  max|.| %.2e\n",
              mean(d_inte), sd(d_inte), max(abs(d_inte))))
  if (nrow(ev_all)) {
    cat(sprintf("  per augmented node (n = %d): max|lam_nh - lam_credit| = %.3g, max|mu_draw - mu_nh| = %.3g, max|mu_draw - mu_credit| = %.3g, max|pd_s - pd_final| = %.3g\n",
                nrow(ev_all), max(abs(ev_all$lam_nh - ev_all$lam_credit)), max(abs(ev_all$mu_draw - ev_all$mu_nh)),
                max(abs(ev_all$mu_draw - ev_all$mu_credit)), max(abs(ev_all$pd_sampler - ev_all$pd_final))))
    cat(sprintf("  augmented nodes with parent_id == -1 (tip_start reset to 0 by compute_pendant_pd): %d of %d\n",
                sum(ev_all$parent_id == -1), nrow(ev_all)))
  }
  ## fhat with the C++ logg vs with the sampler's true density (same trees)
  lw_cpp  <- E$logf - E$logg
  lw_fix  <- E$logf - (E$logg - d_true_cppparent)          # H6 corrected, H5 term kept
  lw_full <- E$logf - (E$logg - d_true)                    # H6 + H5 corrected
  ess <- function(lw) { w <- exp(lw - max(lw)); sum(w)^2/sum(w^2) }
  cat(sprintf("  fhat: C++ %.4f (ESS %.0f) | H6-corrected %.4f (ESS %.0f) | H5+H6-corrected %.4f (ESS %.0f)\n",
              lmexp(lw_cpp), ess(lw_cpp), lmexp(lw_fix), ess(lw_fix), lmexp(lw_full), ess(lw_full)))
  invisible(list(E = E, ev = ev_all, d_true = d_true, d6 = d_true_cppparent, lw_cpp = lw_cpp, lw_fix = lw_fix))
}

set.seed(1)
## 1. validation of q_true on models where the hypothesis predicts agreement
r_cr <- analyse("CR exp-link",     c(log(0.6), 0, 0, 0, log(0.25), 0, 0, 0), c(0L,0L,0L), 1L)
r_dd <- analyse("DD linear N-only", c(0.9, -0.05, 0, 0, 0.15, 0, 0, 0),        c(1L,0L,0L), 0L)
## 2. M-only (beta_M, gamma_M), exponential link
r_m  <- analyse("M-only exp",       c(log(0.5), 0, 0.15, 0, log(0.25), 0, 0.5, 0), c(0L,1L,0L), 1L)
## 3. D-only ("d"), exponential link — the ep_exp compensator branch
r_d  <- analyse("D-only (d) exp",   c(log(0.5), 0, 0, 0.2, log(0.25), 0, 0, 0.4), c(0L,0L,1L), 1L)
## 4. D-only, linear link — the non-exp branch with speciation_rate_ep
r_dl <- analyse("D-only (d) linear", c(0.6, 0, 0, 0.1, 0.25, 0, 0, 0.15),        c(0L,0L,1L), 0L)
## 5. mu-side only: gamma_M != 0, beta_M = 0 — isolates (a) vs (c) on the lifetime density
r_gm <- analyse("M-only exp, gamma_M only", c(log(0.5), 0, 0, 0, log(0.25), 0, 0.8, 0), c(0L,1L,0L), 1L)

## ---- (a) empirically: which mu did the lifetimes come from? -----------------
## gamma_M = 0.8 with pd = 0 makes mu_draw = 0.25 for every node, while the
## credited mu_credit = 0.25*exp(0.8*M) varies with M. Compare the observed
## mean lifetime with its expectation under each candidate rate.
ev <- r_gm$ev
etrunc <- function(mu, r) 1/mu - r*exp(-mu*r)/(1 - exp(-mu*r))
cat(sprintf("\n(a) lifetimes (n = %d): observed mean %.4f | E under mu_draw (pd=0) %.4f | E under mu_credit %.4f | mu_credit range [%.3f, %.3f]\n",
            nrow(ev), mean(ev$life), mean(etrunc(ev$mu_draw, 5 - ev$t_spec)), mean(etrunc(ev$mu_credit, 5 - ev$t_spec)),
            min(ev$mu_credit), max(ev$mu_credit)))
se <- sd(ev$life)/sqrt(nrow(ev)); cat(sprintf("    s.e. of observed mean = %.4f\n", se))

## ---- theta-dependence of the gap and of fhat --------------------------------
cat("\ntheta-dependence (M-only exp, N = 400 per point): gap = fhat_cpp - fhat_H6corrected on the same trees\n")
grid <- expand.grid(beta_M = c(0, 0.15, 0.3), gamma_M = c(0, 0.5, 1.0))
res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(k) {
  p8 <- c(log(0.5), 0, grid$beta_M[k], 0, log(0.25), 0, grid$gamma_M[k], 0)
  r <- analyse(sprintf("grid bM=%.2f gM=%.2f", grid$beta_M[k], grid$gamma_M[k]), p8, c(0L,1L,0L), 1L, N = 400)
  gap <- lmexp(r$lw_cpp) - lmexp(r$lw_fix)
  bs <- replicate(200, { j <- sample(length(r$lw_cpp), replace = TRUE); lmexp(r$lw_cpp[j]) - lmexp(r$lw_fix[j]) })
  data.frame(beta_M = grid$beta_M[k], gamma_M = grid$gamma_M[k], fhat_cpp = lmexp(r$lw_cpp), fhat_fix = lmexp(r$lw_fix),
             gap = gap, gap_bs_se = sd(bs), mean_node_gap = mean(r$d6), max_node_gap = max(abs(r$d6)))
}))
print(res, digits = 4)

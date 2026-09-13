## H8 replicate: independent check on a different tree (12 tips), different
## theta (beta_N != 0), both samplers (linear- and exp-link augmentation), and
## fully observed trees (0 missing lineages).
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressMessages(library(emphasis))
T_TIP <- 10e10; T_EXT <- 0
is_ext <- function(t) t == T_EXT
rate_fun <- function(link) if (link == 0) function(b0, eta) pmax(0, b0 + eta) else if (link == 1) function(b0, eta) exp(b0 + eta) else function(b0, eta) b0 * exp(-0.5 * (eta - 1)^2)
node_rates <- function(p, node, link) {
  f <- rate_fun(link)
  M <- if (node$n > 0) node$pd / node$n else 0
  E <- if (is_ext(node$t_ext)) node$brts - node$tip_start else if (node$parent_id >= 0) node$brts - node$focal_tip_start else M
  D <- E - M
  c(lambda = f(p[1], p[2] * node$n + p[3] * M + p[4] * D), mu = f(p[5], p[6] * node$n + p[7] * M + p[8] * D))
}
## exact integral of max(0, a + b t) via fine trapezoid + kink handling (independent of H8.R's int_relu)
int_relu_num <- function(a, b, t1, t2) {
  sapply(a, function(ai) {
    ks <- -ai / b; pts <- sort(unique(c(t1, t2, if (is.finite(ks) && ks > t1 && ks < t2) ks)))
    s <- 0
    for (j in seq_len(length(pts) - 1)) { u1 <- pts[j]; u2 <- pts[j+1]; v1 <- max(0, ai + b*u1); v2 <- max(0, ai + b*u2); s <- s + 0.5*(v1+v2)*(u2-u1) }
    s })
}
logf_R <- function(p, tr, link, compensator = c("code", "exact")) {
  compensator <- match.arg(compensator); f <- rate_fun(link); nn <- nrow(tr)
  ev <- 0; inte <- 0; prev <- 0; alive_ts <- c(0, 0)
  for (i in seq_len(nn)) {
    node <- tr[i, ]; dt <- node$brts - prev; r <- node_rates(p, node, link)
    if (dt > 0) {
      if (compensator == "code") inte <- inte + dt * node$n * (r[["lambda"]] + r[["mu"]]) else {
        M <- if (node$n > 0) node$pd / node$n else 0; ts <- alive_ts
        if (link == 0) {
          inte <- inte + sum(int_relu_num(p[1] + p[2]*node$n + p[3]*M - p[4]*(ts + M), p[4], prev, node$brts)) +
                         sum(int_relu_num(p[5] + p[6]*node$n + p[7]*M - p[8]*(ts + M), p[8], prev, node$brts))
        } else {
          tot <- function(t) sapply(t, function(tt) { Dv <- (tt - ts) - M
            sum(f(p[1], p[2]*node$n + p[3]*M + p[4]*Dv)) + sum(f(p[5], p[6]*node$n + p[7]*M + p[8]*Dv)) })
          inte <- inte + integrate(tot, prev, node$brts, rel.tol = 1e-10, subdivisions = 2000L)$value
        }
      }
    }
    if (is_ext(node$t_ext)) { ev <- ev + log(max(r[["mu"]], 1e-300)); k <- which(abs(alive_ts - node$tip_start) < 1e-9)[1]; if (is.na(k)) stop("unknown lineage"); alive_ts <- alive_ts[-k] }
    else if (i != nn) { ev <- ev + log(r[["lambda"]]); alive_ts <- c(alive_ts, node$brts) }
    prev <- node$brts
  }
  c(loglik = ev - inte, inte = inte)
}
cpp_logf <- function(p, tr, link, model = c(0L,0L,1L)) emphasis:::eval_logf(p, list(tr), model = as.integer(model), link = as.integer(link), rho = 1)$logf

brts <- c(15, 12.1, 10.3, 9.0, 7.7, 6.2, 5.1, 4.4, 3.0, 2.2, 1.1)   # 12 tips
p_gen <- c(0.4, -0.01, 0, 0.05, 0.1, 0, 0, 0.02)                        # nd-style, beta_N != 0
for (samp in c(0L, 1L)) {
  aug <- emphasis:::augment_trees(brts, p_gen, sample_size = 30, maxN = 4000, max_missing = 300, max_lambda = 1e6,
                                  num_threads = 1, model = c(0L,0L,1L), link = samp, rho = 1)
  trees <- aug$trees; nmiss <- sapply(trees, function(t) sum(is_ext(t$t_ext)))
  cat(sprintf("\n=== sampler link=%d: %d trees, missing per tree %s ===\n", samp, length(trees), paste(range(nmiss), collapse="-")))
  d_exp  <- sapply(trees, function(t) cpp_logf(p_gen, t, 1) - logf_R(p_gen, t, 1, "exact")[["loglik"]])
  cat(sprintf(" exp-link anchor  max|C++ - exact| = %.2e\n", max(abs(d_exp))))
  d_linc <- sapply(trees, function(t) cpp_logf(p_gen, t, 0) - logf_R(p_gen, t, 0, "code")[["loglik"]])
  d_line <- sapply(trees, function(t) cpp_logf(p_gen, t, 0) - logf_R(p_gen, t, 0, "exact")[["loglik"]])
  cat(sprintf(" linear: max|C++ - code| = %.2e ; C++ - exact: mean %+.3f sd %.3f range [%+.3f, %+.3f]; cor(nmiss) %.2f\n",
              max(abs(d_linc)), mean(d_line), sd(d_line), min(d_line), max(d_line), suppressWarnings(cor(d_line, nmiss))))
  if (any(nmiss == 0)) cat(sprintf(" fully observed trees (0 missing): Delta = %s\n", paste(sprintf("%+.4f", d_line[nmiss==0]), collapse=" ")))
}
## fully observed 12-tip tree built by hand (no augmentation at all)
aug0 <- emphasis:::augment_trees(brts, c(0.4,0,0,0,0.0,0,0,0), sample_size = 1, maxN = 4000, max_missing = 0, max_lambda = 1e6,
                                 num_threads = 1, model = c(0L,0L,1L), link = 0L, rho = 1)
t0 <- aug0$trees[[1]]; stopifnot(sum(is_ext(t0$t_ext)) == 0)
for (bD in c(0.02, 0.05, 0.1)) { p <- p_gen; p[4] <- bD
  cat(sprintf(" observed-only tree, beta_D=%.2f: C++ linear = %.5f, exact = %.5f, Delta = %+.4f  (exp link Delta = %+.1e)\n", bD,
      cpp_logf(p, t0, 0), logf_R(p, t0, 0, "exact")[["loglik"]], cpp_logf(p, t0, 0) - logf_R(p, t0, 0, "exact")[["loglik"]],
      cpp_logf(p, t0, 1) - logf_R(p, t0, 1, "exact")[["loglik"]]))
}
## dd sanity on this tree
cat(sprintf(" dd (model 1,0,0) linear on observed tree: Delta = %.2e\n",
    cpp_logf(c(0.5,-0.02,0,0,0.1,0,0,0), t0, 0, model=c(1L,0L,0L)) - logf_R(c(0.5,-0.02,0,0,0.1,0,0,0), t0, 0, "exact")[["loglik"]]))

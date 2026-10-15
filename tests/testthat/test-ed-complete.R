# On a complete simulated tree the raw-ED likelihood is the one of the definition.
#
# Fair-proportion ED of a lineage s alive at t: the sum, over the stretches of
# branch on the path from the crown to s, of each stretch's length divided by
# the number of lineages alive at t below it.  Between events only pendant
# edges grow, so every ED_s grows at unit rate; at a split or an extinction the
# counts of the stretches above it change.  The reference below computes ED
# from that definition alone, straight off the simulator's L-table (no C++),
# and the complete-data log-likelihood of the exponential ed model from it:
# log rate of the event's lineage just before each event, minus the integral of
# the summed rates with ED growing at unit rate on each segment.

.ed_forest <- function(L, T) {
  b <- T - L[, 1]; b[1:2] <- 0
  list(birth = b, death = ifelse(L[, 4] == -1, Inf, T - L[, 4]),
       parent = c(NA, NA, match(L[-(1:2), 2], L[, 3])), n = nrow(L))
}

.ed_def <- function(f, t, al) {
  ex <- f$birth <= t
  hp <- which(ex & !is.na(f$parent))
  kids <- split(hp, factor(f$parent[hp], levels = seq_len(f$n)))
  A <- integer(f$n)                       # alive lineages in each subtree, at t
  for (i in order(f$birth, decreasing = TRUE)) if (ex[i]) A[i] <- as.integer(al[i]) + sum(A[kids[[i]]])
  stretch <- function(l, x, y, self) {    # lineage l from x to y, cut at its children's births
    cs <- kids[[l]]; cb <- f$birth[cs]
    cuts <- sort(unique(c(x, cb[cb > x & cb < y], y))); s <- 0
    for (j in seq_len(length(cuts) - 1)) {
      nd <- self + sum(A[cs[cb >= cuts[j + 1]]])
      if (nd > 0) s <- s + (cuts[j + 1] - cuts[j]) / nd
    }
    s
  }
  out <- rep(NA_real_, f$n)
  for (s in which(al)) {
    e <- stretch(s, f$birth[s], t, 1L); l <- s
    while (!is.na(f$parent[l])) { p <- f$parent[l]; e <- e + stretch(p, f$birth[p], f$birth[l], as.integer(al[p])); l <- p }
    out[s] <- e
  }
  out
}

.ref_ed_logf <- function(L, T, b0, bE, g0, gE) {
  f <- .ed_forest(L, T)
  ev <- rbind(data.frame(t = f$birth[-(1:2)], type = rep(1L, f$n - 2), row = 3:f$n),
              data.frame(t = f$death[is.finite(f$death)], type = rep(2L, sum(is.finite(f$death))),
                         row = which(is.finite(f$death))))
  ev <- ev[order(ev$t), ]
  I <- function(a, b, dt) if (abs(b) < 1e-12) exp(a) * dt else exp(a) * expm1(b * dt) / b
  times <- c(0, ev$t, T); ll <- 0
  for (k in seq_len(nrow(ev) + 1)) {
    t0 <- times[k]; t1 <- times[k + 1]
    al <- f$birth <= t0 & f$death >= t1
    e0 <- .ed_def(f, t0, al)
    ll <- ll - sum(I(b0 + bE * e0[al], bE, t1 - t0)) - sum(I(g0 + gE * e0[al], gE, t1 - t0))
    if (k > nrow(ev)) break
    r <- ev$row[k]
    if (ev$type[k] == 1L) ll <- ll + b0 + bE * (e0[f$parent[r]] + t1 - t0)
    else                  ll <- ll + g0 + gE * (e0[r] + t1 - t0)
  }
  ll
}

# The complete frame (as in test-d-bookkeeping.R): a birth node per lineage, an
# extinction node per death, the present, all replayed from the L-table.
.ed_complete_frame <- function(L, T) {
  f <- .ed_forest(L, T)
  ev <- rbind(data.frame(t = f$birth[-(1:2)], type = rep(1L, f$n - 2), row = 3:f$n),
              data.frame(t = f$death[is.finite(f$death)], type = rep(2L, sum(is.finite(f$death))),
                         row = which(is.finite(f$death))))
  ev <- ev[order(ev$t), ]; prow <- f$parent
  ids <- c(-2L, -3L, seq_len(f$n - 2L) - 1L)
  ts <- c(0, 0, rep(NA_real_, f$n - 2L)); alive <- c(TRUE, TRUE, rep(FALSE, f$n - 2L))
  out <- vector("list", nrow(ev) + 1L)
  for (k in seq_len(nrow(ev))) {
    t <- ev$t[k]; r <- ev$row[k]; a <- which(alive); n <- length(a); pd <- sum(t - ts[a])
    if (ev$type[k] == 1L) {
      p <- prow[r]
      out[[k]] <- data.frame(brts = t, n = n, t_ext = if (is.finite(f$death[r])) f$death[r] else 10e10,
                             pd = pd, tip_start = t, focal_tip_start = ts[p], clade = 0L,
                             id = ids[r], parent_id = ids[p])
      ts[p] <- t; ts[r] <- t; alive[r] <- TRUE
    } else {
      out[[k]] <- data.frame(brts = t, n = n, t_ext = 0, pd = pd, tip_start = ts[r],
                             focal_tip_start = ts[r], clade = 0L, id = ids[r], parent_id = -1L)
      alive[r] <- FALSE
    }
  }
  a <- which(alive)
  out[[nrow(ev) + 1L]] <- data.frame(brts = T, n = length(a), t_ext = 10e10, pd = sum(T - ts[a]),
                                     tip_start = T, focal_tip_start = -1, clade = 0L,
                                     id = -1L, parent_id = -1L)
  do.call(rbind, out)
}

test_that("on a complete simulated ed tree eval_logf is the log-likelihood of fair-proportion ED", {
  skip_on_cran()
  T <- 7
  mb <- emphasis:::.resolve_model("ed")
  gens <- list(c(log(0.8), -0.25, log(0.08), 0), c(log(0.25), 0.17, log(0.03), 0))
  worst <- 0; n_trees <- 0L; n_ext <- 0L
  set.seed(44)
  for (g in gens) for (rep in 1:4) {
    s <- simulate_tree(pars = g, max_t = T, model = "ed", link = "exponential", rho = 1,
                       max_lin = 300L, max_tries = 20L, num_threads = 1L)
    if (!identical(s$status, "done")) next
    df <- .ed_complete_frame(s$L, T)
    for (q in list(g, g + c(0.2, -0.1, 0.3, 0.05), g + c(-0.3, 0.08, -0.2, -0.04))) {
      got <- eval_logf(emphasis:::.expand_pars(q, mb), list(df), model = mb, link = 1L, rho = 1)$logf
      ref <- .ref_ed_logf(s$L, T, q[1], q[2], q[3], q[4])
      worst <- max(worst, abs(got - ref))
    }
    n_trees <- n_trees + 1L; n_ext <- n_ext + sum(df$t_ext == 0)
  }
  expect_gt(n_trees, 5L)
  expect_gt(n_ext, 0L)
  expect_lt(worst, 1e-8)
})

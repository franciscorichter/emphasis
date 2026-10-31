# The mark-sampling proposal for latent marks (dev/latent-marks/design.md, phase 3); internal.
#
# The proposal is the exact posterior of a nearby model in which a hidden lineage keeps the mark of
# the observed branch it hangs from, and so do its descendants: marks change only at observed
# splits.  Under that model, with a mark s on an observed branch and constant rates
# lambda_s = lambda_0 c_s, mu_s = mu_0 c_s along it,
#   - hidden subtrees attach along the branch as a Poisson process of intensity 2 lambda_s q_s(t),
#     q_s(t) the probability that a lineage born at t leaves no descendant at T; the 2 counts the
#     two orientations of a split (either daughter may be the one that survives), which the
#     L-table records as one, the new row being the hidden one;
#   - each hidden subtree is a birth-death process conditioned to die before T: a lineage at u has
#     birth rate lambda_s q_s(u) and death rate mu_s / q_s(u);
#   - the branch's factor in the reconstructed likelihood is exp(-int (mu_s + lambda_s (1 - 2 q_s)))
#     (with sigma = 0 the product over branches is the textbook p1 likelihood).
# The marks are drawn first, jointly over the observed tree, by forward filtering and backward
# sampling (FFBS) over the grid; then the hidden lineages given the marks.  The marginal density of
# the hidden lineages is the ratio of two pruning passes over the observed tree,
#   g(h) = p(obs, h) / p(obs),
# each summing over the marks, so an importance weight needs no mixture:  log w = log L_marks(h) -
# log g(h), with L_marks the latent-rate model's complete-tree likelihood (.eval_logf_marks).
# With sigma = 0 the proposal is the exact posterior of cr.

# ---- closed forms for cr on one branch (forward time t from the crown; s = T - t)
.pm_q <- function(t, lam, mu, T) {          # P(a lineage born at t has no descendant at T)
  s <- T - t; d <- lam - mu
  if (abs(d) < 1e-10 * lam) return(lam * s / (1 + lam * s))
  mu * (1 - exp(-d * s)) / (lam - mu * exp(-d * s))
}
# int_a^b lambda (1 - q(t)) dt
.pm_int_l1q <- function(a, b, lam, mu, T) {
  sa <- T - a; sb <- T - b; d <- lam - mu
  if (abs(d) < 1e-10 * lam) return(log((1 + lam * sa) / (1 + lam * sb)))
  if (d > 0) return(d * (sa - sb) + log((lam - mu * exp(-d * sa)) / (lam - mu * exp(-d * sb))))
  log((lam * exp(d * sa) - mu) / (lam * exp(d * sb) - mu))
}
# int_a^b mu / q(t) dt  (diverges as b -> T)
.pm_int_muq <- function(a, b, lam, mu, T) {
  sa <- T - a; sb <- T - b; d <- lam - mu
  if (abs(d) < 1e-10 * lam) return(lam * (sa - sb) + log(sa / sb))   # mu = lam: mu / q = lam + 1 / s
  if (d > 0) return(mu * (sa - sb) + d * (sa - sb) + log(-expm1(-d * sa)) - log(-expm1(-d * sb)))
  mu * (sa - sb) + log(expm1(d * sa) / expm1(d * sb))
}

# ---- the observed tree as branches (one per phylo edge), forward times
.pm_branches <- function(phy) {
  T <- max(ape::node.depth.edgelength(phy))
  dep <- ape::node.depth.edgelength(phy)          # distance from the root
  E <- phy$edge; nt <- ape::Ntip(phy)
  br <- data.frame(e = seq_len(nrow(E)), from = E[, 1], to = E[, 2], a = dep[E[, 1]], b = dep[E[, 2]])
  br$tip <- br$to <= nt
  br$parent <- match(br$from, br$to)              # NA for the two crown branches
  kids <- split(br$e, br$from)
  br$c1 <- br$c2 <- NA_integer_
  for (i in which(!br$tip)) { k <- kids[[as.character(br$to[i])]]; br$c1[i] <- k[1]; br$c2[i] <- k[2] }
  # postorder: children before parents
  ord <- integer(0); visit <- function(i) { if (!br$tip[i]) { visit(br$c1[i]); visit(br$c2[i]) }; ord <<- c(ord, i) }
  crown <- which(is.na(br$parent)); for (i in crown) visit(i)
  list(br = br, T = T, post = ord, crown = crown)
}

.lse <- function(x) { m <- max(x); if (!is.finite(m)) return(m); m + log(sum(exp(x - m))) }
# log of (P %*% exp(v)) row-wise, stable
.lPv <- function(lP, v) { m <- max(v); if (!is.finite(m)) return(rep(-Inf, nrow(lP))); log(exp(lP) %*% exp(v - m))[, 1] + m }

# ---- one pruning pass: per branch and state, the log factor `fac[e, s]`; returns the partials and
# the root value.  Internal branches add log lambda_s for their end split and the two children.
.pm_prune <- function(B, fac, lam_s, lP, lpi) {
  S <- ncol(fac); L <- matrix(NA_real_, nrow(B$br), S)
  for (i in B$post) {
    v <- fac[i, ]
    if (!B$br$tip[i]) v <- v + log(lam_s) + .lPv(lP, L[B$br$c1[i], ]) + .lPv(lP, L[B$br$c2[i], ])
    L[i, ] <- v
  }
  root <- lpi + .lPv(lP, L[B$crown[1], ]) + .lPv(lP, L[B$crown[2], ])
  list(L = L, root = .lse(root), rootv = root)
}

# reconstructed factor of each branch per state: -int (mu_s + lambda_s (1 - 2 q_s))
#   = (lambda_s - mu_s) (b - a) - 2 int lambda_s (1 - q_s)
.pm_fac_obs <- function(B, lam_s, mu_s) {
  S <- length(lam_s); br <- B$br
  f <- matrix(0, nrow(br), S)
  for (s in seq_len(S)) f[, s] <- (lam_s[s] - mu_s[s]) * (br$b - br$a) - 2 * .pm_int_l1q(br$a, br$b, lam_s[s], mu_s[s], B$T)
  f
}

# ---- FFBS: a joint draw of the branch marks from the approximate posterior
.pm_ffbs <- function(B, pr, lP, lpi) {
  st <- integer(nrow(B$br))
  r <- sample.int(length(lpi), 1, prob = exp(pr$rootv - max(pr$rootv)))
  draw_child <- function(i, sp) {
    lw <- lP[sp, ] + pr$L[i, ]; st[i] <<- sample.int(length(lw), 1, prob = exp(lw - max(lw)))
    if (!B$br$tip[i]) { draw_child(B$br$c1[i], st[i]); draw_child(B$br$c2[i], st[i]) }
  }
  for (i in B$crown) draw_child(i, r)
  st
}

# ---- hidden subtrees on one branch with rates (lam, mu): rows (birth, death, parent_row) where
# parent_row 0 is the observed branch itself and k > 0 the k-th hidden row
.pm_hidden_branch <- function(a, b, lam, mu, T) {
  rows <- list(); n <- 0L
  # attachments: thinning of a rate-2 lam Poisson process by q(t)
  t <- a
  repeat {
    t <- t + stats::rexp(1, 2 * lam); if (t >= b) break
    if (stats::runif(1) < .pm_q(t, lam, mu, T)) {
      n <- n + 1L; rows[[n]] <- c(t, NA, 0)
      stack <- n
      while (length(stack)) {
        k <- stack[1]; stack <- stack[-1]; u <- rows[[k]][1]
        repeat {
          E <- stats::rexp(1)
          H <- function(x) lam * (x - u) - .pm_int_l1q(u, x, lam, mu, T) + .pm_int_muq(u, x, lam, mu, T) - E
          hi <- T - 1e-12 * max(1, T)
          tau <- if (H(hi) <= 0) hi else stats::uniroot(H, c(u, hi), tol = 1e-12)$root
          qq <- .pm_q(tau, lam, mu, T)
          pb <- lam * qq / (lam * qq + mu / qq)
          if (stats::runif(1) < pb) {          # birth: a daughter, the lineage continues
            n <- n + 1L; rows[[n]] <- c(tau, NA, k); stack <- c(stack, n); u <- tau
          } else { rows[[k]][2] <- tau; break }
        }
      }
    }
  }
  if (!n) return(matrix(numeric(0), 0, 3))
  do.call(rbind, rows)
}

# ---- complete factor of each branch per state, given its hidden subtrees: the cr complete-data
# log-likelihood of the branch and everything hanging from it, rates (lam_s, mu_s), for the
# orientation the L-table records.  The proposal density of a draw is 2^k times the ratio of the
# pruned values (k attachments on observed branches); the target carries the same 2^k (the
# likelihood of the latent-rate model does not depend on the orientation), so it cancels in the
# weight and is left out of both.
.pm_fac_aug <- function(B, hid, lam_s, mu_s) {
  S <- length(lam_s); br <- B$br; f <- matrix(0, nrow(br), S)
  for (i in seq_len(nrow(br))) {
    h <- hid[[i]]; len <- br$b[i] - br$a[i]; nb <- 0; nd <- 0; hl <- 0
    if (nrow(h)) { nb <- nrow(h); nd <- nrow(h); hl <- sum(h[, 2] - h[, 1]) }
    f[i, ] <- -(lam_s + mu_s) * (len + hl) + nb * log(lam_s) + nd * log(mu_s)
  }
  f
}

# ---- the observed tree plus hidden lineages as an L-table (DDD layout: birth age, parent label,
# label, death age or -1) and its complete-tree frame for eval_logf
.pm_ltable <- function(B, hid) {
  br <- B$br; T <- B$T
  lab <- integer(nrow(br)); rows <- list(); nl <- 0L
  # the crown: two lineages, labels -1 and 2 (DDD's crown convention), then each internal node
  # continues its first child on the parent's label and gives the second a new one
  next_lab <- 3L
  newrow <- function(birth_t, parent, label, death_t) { nl <<- nl + 1L; rows[[nl]] <<- c(T - birth_t, parent, label, if (is.na(death_t)) -1 else T - death_t) }
  assign_lab <- function(i, l) {
    lab[i] <<- l
    if (!br$tip[i]) {
      assign_lab(br$c1[i], l)
      l2 <- next_lab; next_lab <<- next_lab + 1L
      newrow(br$b[i], l, l2, NA); assign_lab(br$c2[i], l2)
    }
  }
  newrow(0, 0, -1, NA); newrow(0, -1, 2, NA)
  assign_lab(B$crown[1], -1L); assign_lab(B$crown[2], 2L)
  for (i in seq_len(nrow(br))) {
    h <- hid[[i]]; if (!nrow(h)) next
    hl <- integer(nrow(h))
    for (k in seq_len(nrow(h))) {
      hl[k] <- next_lab; next_lab <- next_lab + 1L
      par <- if (h[k, 3] == 0) lab[i] else hl[h[k, 3]]
      newrow(h[k, 1], par, hl[k], h[k, 2])
    }
  }
  L <- do.call(rbind, rows); L[1, 1] <- T; L[2, 1] <- T
  L
}

# ---- n draws at parameters (lambda_0, mu_0) with the marks of a preset; returns the trees (frames
# for eval_logf), log g and the per-draw counts
.pm_draws <- function(phy, lambda0, mu0, marks, n) {
  B <- .pm_branches(phy)
  lam_s <- lambda0 * exp(marks$log_clam); mu_s <- mu0 * exp(marks$log_cmu)
  lP <- if (marks$S > 1L) log(marks$P) else matrix(0, 1, 1); lpi <- log(marks$pi)
  po <- .pm_prune(B, .pm_fac_obs(B, lam_s, mu_s), lam_s, lP, lpi)
  out <- vector("list", n); logg <- numeric(n); nh <- integer(n)
  for (j in seq_len(n)) {
    st <- .pm_ffbs(B, po, lP, lpi)
    hid <- lapply(seq_len(nrow(B$br)), function(i) .pm_hidden_branch(B$br$a[i], B$br$b[i], lam_s[st[i]], mu_s[st[i]], B$T))
    pa <- .pm_prune(B, .pm_fac_aug(B, hid, lam_s, mu_s), lam_s, lP, lpi)
    logg[j] <- pa$root - po$root
    L <- .pm_ltable(B, hid)
    out[[j]] <- .ltable_complete_frame(L, B$T)
    nh[j] <- sum(vapply(hid, nrow, 1L))
  }
  list(trees = out, logg = logg, n_hidden = nh, log_pobs = po$root)
}

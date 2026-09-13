suppressWarnings(suppressMessages(devtools::load_all(".", quiet = TRUE)))
ns <- asNamespace("emphasis")

# Pendant PD from the edge list, independent of the code under test.
# Alive just before t: an edge that started strictly before t and has not
# ended before t.  That is the lineage set the node's own `n` counts.
true_P <- function(phy, t) {
  h  <- ape::node.depth.edgelength(phy)
  st <- h[phy$edge[, 1]]; en <- h[phy$edge[, 2]]
  a  <- st < t & en >= t
  sum(t - st[a])
}
true_N <- function(phy, t) {
  h  <- ape::node.depth.edgelength(phy)
  st <- h[phy$edge[, 1]]; en <- h[phy$edge[, 2]]
  sum(st < t & en >= t)
}
true_ts <- function(phy, t) {
  h  <- ape::node.depth.edgelength(phy)
  st <- h[phy$edge[, 1]]; en <- h[phy$edge[, 2]]
  st[st < t & en >= t]
}

set.seed(4)
phy  <- ape::rphylo(12, 0.6, 0)
brts <- ns$.extract_brts(phy)
pts  <- ns$.pts(brts)
cat("n brts:", length(brts), " n pts:", length(pts), "\n")

a <- augment_trees(as.numeric(brts), c(0.3, 0, 0, 0, 0, 0, 0, 0),
                   1L, 500L, 100L, 1e6, 1L,
                   model = c(0L, 0L, 0L), link = 0L, rho = 1,
                   parent_tip_start = pts)
df <- a$trees[[1]]
cat("nodes:", nrow(df), " (observed:", length(brts), ")\n")
print(df)

M_pkg <- df$pd / df$n
E_pkg <- df$brts - df$focal_tip_start
ev    <- seq_len(nrow(df) - 1L)     # last node is the terminal marker

M_true <- vapply(df$brts[ev], function(t) true_P(phy, t) / true_N(phy, t), 0)
N_true <- vapply(df$brts[ev], function(t) true_N(phy, t), 0)
cat("\nmax |n_pkg - N_true|   :", max(abs(df$n[ev] - N_true)), "\n")
cat("max |M_pkg - M_true|   :", max(abs(M_pkg[ev] - M_true)), "\n")
cat("M_pkg / t              :", paste(round(M_pkg[ev] / df$brts[ev], 3), collapse = " "), "\n")

D_true <- vapply(seq_along(ev), function(k) {
  t <- df$brts[k]
  (t - pts[k]) - M_true[k]
}, 0)
D_pkg <- E_pkg[ev] - M_pkg[ev]
cat("max |D_pkg - D_true|   :", max(abs(D_pkg - D_true)), "\n")

# sum_s D(s, t) = 0 over the alive lineages, using the package's own M
sd0 <- vapply(ev, function(k) {
  t <- df$brts[k]
  sum((t - true_ts(phy, t)) - M_pkg[k])
}, 0)
cat("max |sum_s D(s,t)|     :", max(abs(sd0)), "\n")

cat("\n--- no-topology fallback (bare branching times) ---\n")
a0 <- augment_trees(as.numeric(brts), c(0.3, 0, 0, 0, 0, 0, 0, 0),
                    1L, 500L, 100L, 1e6, 1L,
                    model = c(0L, 0L, 0L), link = 0L, rho = 1)
d0 <- a0$trees[[1]]
cat("M0 / t :", paste(round((d0$pd / d0$n) / d0$brts, 3), collapse = " "), "\n")
cat("tip_start:", paste(round(d0$tip_start, 3), collapse = " "), "\n")
cat("focal    :", paste(round(d0$focal_tip_start, 3), collapse = " "), "\n")
cat("D0       :", paste(round((d0$brts - d0$focal_tip_start) - d0$pd / d0$n, 6), collapse = " "), "\n")

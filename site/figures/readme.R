# The README's figures, drawn from the package itself.
#   Rscript site/figures/readme.R          (from the repository root; writes man/figures/readme-*.png)
# Colours are the docs site's (site/style.css): ink, accent, and one hue for what is hidden.

devtools::load_all(quiet = TRUE)
library(ape)

INK    <- "#16201d"
INK2   <- "#3f4a45"
ACCENT <- "#1d6d5a"
HIDDEN <- "#6c63b5"
RULE   <- "#dde3de"
OUT    <- "man/figures"

png_out <- function(name, w, h) {
  ragg::agg_png(file.path(OUT, name), width = w, height = h, units = "in", res = 200, background = "white")
}

# Edges with no descendant alive at the present: extinct (or, in an augmented tree, hidden).
hidden_edges <- function(tr, t_present) {
  depth <- node.depth.edgelength(tr)
  n <- Ntip(tr)
  alive_tip <- depth[seq_len(n)] > t_present - 1e-6
  alive_node <- logical(n + tr$Nnode); alive_node[seq_len(n)] <- alive_tip
  for (e in postorder(tr)) {   # children before parents
    if (alive_node[tr$edge[e, 2]]) alive_node[tr$edge[e, 1]] <- TRUE
  }
  !alive_node[tr$edge[, 2]]
}

# hidden_col = NA draws the tree with its extinct lineages invisible: the reconstructed tree
# in exactly the layout of the full history, so the two differ only by what is missing.
draw_tree <- function(tr, t_present, title, sub, hidden_col) {
  tr <- ladderize(tr, right = FALSE)
  h <- hidden_edges(tr, t_present)
  plot(tr, show.tip.label = FALSE, edge.color = ifelse(h, if (is.na(hidden_col)) "#FFFFFF00" else hidden_col, INK),
       edge.width = 1.5, x.lim = c(0, t_present), no.margin = FALSE, root.edge = FALSE)
  mtext(title, side = 3, line = 1.4, adj = 0, col = INK, font = 2, cex = 0.95)
  coloured_sub(sub)
  axis(1, at = c(0, t_present), labels = c("crown", "present"), col = RULE, col.axis = INK2,
       cex.axis = 0.9, tck = -0.02, mgp = c(3, 0.4, 0), padj = 0)
}

# A subtitle in pieces, each piece in its own colour: list(c(text, colour), ...).
# The colour of the words replaces a legend.
coloured_sub <- function(pieces, cex = 0.8) {
  x <- par("usr")[1]
  for (pc in pieces) {
    mtext(pc[1], side = 3, line = 0.25, at = x, adj = 0, col = pc[2], cex = cex, font = if (pc[2] == INK2) 1 else 2)
    # mtext's cex is absolute; strwidth's is relative to par("cex")
    x <- x + strwidth(pc[1], cex = cex / par("cex"), font = if (pc[2] == INK2) 1 else 2)
  }
}

# ---- Figure 1: what a reconstructed tree hides, and one draw of it ---------------------
truth <- c(log(0.7), -0.012, log(0.3), 0)       # dd: lambda = 0.7 exp(-0.012 N), mu = 0.3
# the first seed whose clade has 35-60 survivors (the compiled sampler's draws differ between
# platforms, so the seed is searched for rather than fixed)
for (seed in 1:500) {
  set.seed(seed)
  sim <- simulate_tree(pars = truth, max_t = 10, model = "dd")
  if (identical(sim$status, "done") && Ntip(sim$tes) >= 35 && Ntip(sim$tes) <= 60) break
}
obs <- sim$tes
aug <- NULL
for (k in 1:50) {
  aug <- suppressWarnings(simulate_tree(obs, pars = truth, model = "dd", n_trees = 1, method = "bdi"))$tas
  if (!is.null(aug)) break
}
stopifnot(!is.null(aug))
cat("figure 1: seed", seed, "survivors", Ntip(obs), "\n")

n_obs <- Ntip(obs); n_ext <- Ntip(sim$tas) - n_obs; n_aug <- Ntip(aug) - n_obs
png_out("readme-hidden.png", 10, 3.8)
par(mfrow = c(1, 3), mar = c(2.2, 1.2, 3.4, 1.6), family = "sans")
draw_tree(sim$tas, 10, "The full history",
          list(c(sprintf("%d survivors", n_obs), INK), c(" and ", INK2), c(sprintf("%d extinct lineages", n_ext), HIDDEN)), HIDDEN)
draw_tree(sim$tas, 10, "The data: the reconstructed tree",
          list(c(sprintf("the same %d survivors", n_obs), INK), c(", nothing else", INK2)), NA)
draw_tree(aug, 10, "One draw of the hidden lineages",
          list(c(sprintf("%d extinct lineages", n_aug), HIDDEN), c(sprintf(" proposed (the history had %d)", n_ext), INK2)), HIDDEN)
dev.off()

# ---- Figure 2: the covariate path on a tree with an age effect -------------------------
# speciation lambda = 0.5 exp(-0.4 D): lineages that have waited longer than the clade's
# mean speciate less. covariate_path() orders the covariates by when they enter.
for (seed in 1:500) {
  set.seed(seed)
  s <- simulate_tree(pars = c(log(0.5), -0.4, log(0.05), 0), max_t = 12, model = "d")
  if (identical(s$status, "done") && Ntip(s$tes) >= 120 && Ntip(s$tes) <= 250) break
}
cp <- covariate_path(s$tes, covariates = c("N", "D", "ED", "K"), criterion = "BIC")
P <- cp$speciation$path
cat("figure 2: seed", seed, "tips", Ntip(s$tes), "entry", cp$speciation$entry_order,
    "active", cp$speciation$active, "D at choice", round(cp$speciation$coef[["D"]], 3), "\n")

png_out("readme-path.png", 7, 3.8)
par(mar = c(3.2, 4.2, 3.6, 3), mgp = c(2, 0.6, 0), family = "sans")
covs <- c("N", "D", "ED", "K")
cols <- c(N = "#9aa5a0", D = ACCENT, ED = "#3f4a45", K = "#6f7b76")
yl <- range(c(as.matrix(P[, covs]), -0.4), na.rm = TRUE) * 1.08
plot(NA, xlim = rev(range(P$gamma)), ylim = yl, axes = FALSE,
     xlab = "dgLARS tuning parameter (decreasing: more covariates enter)", ylab = "speciation coefficient",
     col.lab = INK2, cex.lab = 0.82)
abline(h = -0.4, col = ACCENT, lty = 3, lwd = 1.5)
g_choice <- P$gamma[which.min(P$BIC)]
abline(v = g_choice, col = INK2, lty = 2)
for (v in covs) lines(P$gamma, P[[v]], col = cols[[v]], lwd = if (v == "D") 2.6 else 1.6)
last <- nrow(P)
# labels at the left, where each line is still at zero except D; spread so they do not collide
# N, ED and K: each labelled at the end of the path, spread vertically so the labels do not collide
others <- setdiff(covs, "D")
ends <- vapply(others, function(v) P[[v]][last], 0)
ord <- order(ends); ylab_pos <- ends[ord]
for (i in seq_along(ylab_pos)[-1]) ylab_pos[i] <- max(ylab_pos[i], ylab_pos[i - 1] + 0.05)
for (i in seq_along(ord)) mtext(others[ord][i], side = 4, at = ylab_pos[i], las = 1, line = 0.3, col = cols[[others[ord][i]]], font = 2, cex = 0.8)
text(max(P$gamma), 0.02, "N, ED and K: zero until after the BIC choice", adj = c(0, 0), col = INK2, cex = 0.75)
mid <- which.min(abs(P$gamma - mean(range(P$gamma))))
text(P$gamma[mid], P$D[mid], "D", pos = 1, col = ACCENT, font = 2, cex = 0.9, offset = 0.5)
text(g_choice, yl[1], "chosen by BIC", pos = 2, col = INK2, cex = 0.75, offset = 0.3)
d_at <- P$D[which.min(P$BIC)]
points(g_choice, d_at, pch = 21, bg = ACCENT, col = "white", cex = 1.3)
text(g_choice, d_at, sprintf("D = %.2f", d_at), adj = c(1.15, 1.6), col = ACCENT, font = 2, cex = 0.8)
text(max(P$gamma), -0.4, "true effect of D", adj = c(0, -0.5), col = ACCENT, cex = 0.75)
axis(1, col = RULE, col.axis = INK2, cex.axis = 0.75, tck = -0.02)
axis(2, col = RULE, col.axis = INK2, cex.axis = 0.75, tck = -0.02, las = 1)
mtext(sprintf("%s enters first, and BIC keeps %s", cp$speciation$entry_order[1], paste(cp$speciation$active, collapse = " + ")), side = 3, line = 1.6, adj = 0, col = INK, font = 2, cex = 0.95)
mtext(sprintf("a %d-tip tree simulated with speciation 0.5 exp(-0.4 D); the path over N, D, ED and K", Ntip(s$tes)),
      side = 3, line = 0.5, adj = 0, col = INK2, cex = 0.75)
dev.off()

# capability-map.R -- the capability map of the overview page (section 10).
#
# What the experiment record says emphasis can estimate, by quantity and tree size. Every cell
# is typed by hand from the experiment it names; nothing is computed. Run from a directory
# that holds the experiments repository's sim/figstyle.R and sim/palette*.R:
#
#   Rscript capability-map.R <out.png>
#
# Categories (ordered): reliable (the experiment's rule is met), partly (met for one sign, one
# direction or a subset of fits), not estimable (trees of that size do not carry the quantity),
# not measured (no current measurement; the label names the experiment that will make it).
# A ring marks a cell whose measurement E45 is repeating under the conditioned likelihood.
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
source("sim/figstyle.R")
args <- commandArgs(TRUE)
out <- if (length(args)) args[1] else "capability-map.png"

sizes <- c(25, 50, 100, 200, 400, 800, 1600)
LEV <- c("reliable", "partly", "not estimable", "not measured")
FILL <- c("reliable" = DIVERGING[7], "partly" = DIVERGING[5],
          "not estimable" = DIVERGING[2], "not measured" = PANEL)
TXT <- c("reliable" = "#FFFFFF", "partly" = INK, "not estimable" = INK, "not measured" = GREY)

cell <- function(row, size, cat, lab, ring = FALSE)
  data.frame(row = row, size = size, cat = cat, lab = lab, ring = ring)
nm <- function(row, sz, lab) do.call(rbind, lapply(sz, function(s) cell(row, s, "not measured", lab)))
R <- "reliable"; P <- "partly"; NI <- "not estimable"

# ---- panel A: one row per question a user asks
A_rows <- c("Speciation rate (cr)",
            "Diversity dependence of speciation",
            "ED effect on speciation (centred)",
            "K effect on speciation",
            "D effect on speciation",
            "Extinction rate, constant rates, turnover 0.10",
            "Extinction slopes on covariates",
            "Choice among the ladder's models")
r <- A_rows
A <- rbind(
  # E1: MCEM at the exact maximum under cr and dd, 20-100 tips (unconditioned, linear link);
  # E42: the same on real trees of 50-130 tips.
  cell(r[1], 25, R, "E1", TRUE), cell(r[1], 50, R, "E1", TRUE), cell(r[1], 100, R, "E1", TRUE),
  nm(r[1], c(200, 400, 800), "E45"), nm(r[1], 1600, ""),
  cell(r[2], 25, R, "E1", TRUE), cell(r[2], 50, R, "E1\nE42", TRUE),
  cell(r[2], 100, R, "E1\nE42", TRUE), nm(r[2], c(200, 400), "E45"), nm(r[2], c(800, 1600), ""),
  # E24: edc recovered (0.96, 1.02) and selected (0.95, 0.85), no false positive, 100 tips;
  # E11: ned selected on 0.05 and 0.00 of null trees near 100 tips.
  nm(r[3], 25, ""), nm(r[3], 50, "E46"), cell(r[3], 100, R, "E24\nE11"),
  nm(r[3], c(200, 400), "E46"), nm(r[3], c(800, 1600), ""),
  # E24: k meets every half at -0.35, sampler-limited at +0.35; E41: nk never meets the stopping
  # rule on real trees of 50-65 tips.
  nm(r[4], 25, ""), cell(r[4], 50, P, "E41"), cell(r[4], 100, P, "E24"),
  nm(r[4], c(200, 400), "E46"), nm(r[4], c(800, 1600), ""),
  # D: only complete trees are current (E44); every reconstructed-tree D result is superseded.
  nm(r[5], 25, ""), nm(r[5], c(50, 100, 200, 400), "E46"), nm(r[5], c(800, 1600), ""),
  # E41: 0 of 16 real trees exclude zero; E32: bd_ML at zero on 37 of 40; E38: detected on
  # at most 33% of trees up to 1600 tips.
  nm(r[6], 25, "E45"), cell(r[6], 50, NI, "E41"),
  cell(r[6], 100, NI, "E38", TRUE), cell(r[6], 200, NI, "E32"),
  cell(r[6], 400, NI, "E38", TRUE), cell(r[6], 800, NI, "E38", TRUE),
  cell(r[6], 1600, NI, "E38"),
  # E41: slopes run to the box edge; E24: fitted slopes far from their generating zero;
  # E32/E35: slope on the box limit at 0.10, not estimable by rule at 0.25.
  nm(r[7], 25, ""), cell(r[7], 50, NI, "E41"), cell(r[7], 100, NI, "E24"),
  cell(r[7], 200, NI, "E32\nE35"), nm(r[7], 400, "E46"), nm(r[7], c(800, 1600), ""),
  # E24 and E11 at 100 tips: edc and k selected, no invented ED effect; ned, nd and the
  # ladder's own selection rates are superseded (E47).
  nm(r[8], 25, ""), nm(r[8], 50, "E47"), cell(r[8], 100, P, "E24\nE11"),
  nm(r[8], 200, "E47"), nm(r[8], c(400, 800, 1600), ""))

# ---- panel B: E38's exact profile likelihood. reliable = within a factor of 2 on 80% of
# trees; partly = zero excluded on 80%, not bounded within 2; not estimable = zero not
# excluded on 80%. The cr cells up to 800 tips at turnover up to 0.40 are being remeasured (E45).
B_rows <- c("constant rates, turnover 0.10", "constant rates, turnover 0.25", "constant rates, turnover 0.40",
            "constant rates, turnover 0.50", "constant rates, turnover 0.60",
            "diversity dependence, turnover 0.10-0.40")
b <- function(row, cats) {
  s <- c(100, 200, 400, 800, 1600)
  ring <- row %in% B_rows[1:3] & s <= 800
  rbind(nm(row, c(25, 50), ""), data.frame(row = row, size = s, cat = cats, lab = "", ring = ring))
}
B <- rbind(
  b(B_rows[1], c(NI, NI, NI, NI, NI)),
  b(B_rows[2], c(NI, NI, NI, NI, P)),
  b(B_rows[3], c(NI, NI, P, P, R)),
  b(B_rows[4], c(P, P, P, R, R)),
  b(B_rows[5], c(P, P, R, R, R)),
  nm(B_rows[6], c(25, 50), ""), cell(B_rows[6], 100, R, ""), cell(B_rows[6], 200, R, ""),
  nm(B_rows[6], c(400, 800, 1600), ""))

draw <- function(d, rows, title, subtitle, legend = FALSE) {
  d$x <- match(d$size, sizes)
  d$y <- length(rows) + 1 - match(d$row, rows)
  d$cat <- factor(d$cat, levels = LEV)
  ggplot(d, aes(x, y)) +
    geom_tile(aes(fill = cat), colour = NA, width = 0.94, height = 0.88) +
    geom_text(aes(label = lab, colour = cat), size = 2.45, lineheight = 0.85) +
    geom_point(data = d[d$ring, ], aes(x = x + 0.32, y = y + 0.22, shape = "measured again, conditioned on survival (E45)"),
               size = 1.8, stroke = 0.6, colour = INK, fill = "#FFFFFF") +
    scale_shape_manual(values = 21, name = NULL) +
    scale_fill_manual(values = FILL, drop = FALSE, name = NULL) +
    scale_colour_manual(values = TXT, guide = "none", drop = FALSE) +
    scale_x_continuous(breaks = seq_along(sizes), labels = sizes, position = "top",
                       expand = c(0, 0)) +
    scale_y_continuous(breaks = seq_along(rows), labels = rev(rows), expand = c(0, 0)) +
    coord_cartesian(xlim = c(0.5, length(sizes) + 0.5), ylim = c(0.5, length(rows) + 0.5)) +
    labs(x = "tips in the reconstructed tree", y = NULL, title = title, subtitle = subtitle) +
    theme_emphasis(base_size = 8.5) +
    theme(panel.grid = element_blank(), panel.grid.major = element_blank(), axis.text.y = element_text(colour = INK, size = 7.6),
          axis.text.x = element_text(colour = GREY, size = 7.4),
          axis.title.x = element_text(colour = GREY, size = 7.2, hjust = 0.5),
          plot.title = element_text(size = 8.8, face = "bold"),
          plot.subtitle = element_text(size = 7.2),
          legend.position = if (legend) "top" else "none",
          legend.key.size = unit(9, "pt"), legend.text = element_text(size = 7.4),
          legend.justification = "left", legend.box = "vertical", legend.box.just = "left", legend.location = "plot", legend.box.spacing = unit(2, "pt"),
          legend.margin = margin(0, 0, 2, 0), legend.spacing.x = unit(3, "pt"), legend.spacing.y = unit(0, "pt")) +
    guides(fill = guide_legend(order = 1, nrow = 1, override.aes = list(colour = NA)),
           shape = guide_legend(order = 2))
}

pA <- draw(A, A_rows, "A. Each quantity a user asks about, by tree size",
           "Label: the experiment behind the cell; on a grey cell, the one that will measure it.",
           legend = TRUE)
pB <- draw(B, B_rows, "B. The extinction rate by turnover, constant rates and diversity dependence (exact likelihood, E38)",
           "Same colours. Reliable: within a factor of 2 on 80% of trees; partly: zero excluded on 80%, not bounded.")

p <- (pA / pB) + plot_layout(heights = c(8, 6)) +
  plot_annotation(
    title = "emphasis estimates speciation and its diversity dependence at 25-100 tips;\nunder constant rates at turnover 0.10, extinction is not estimable up to 1600 tips",
    subtitle = "Each cell is what the experiment record says about one quantity at one tree size.",
    theme = theme_emphasis(base_size = 8.5) +
      theme(plot.title = element_text(size = 10, face = "bold", lineheight = 1.05),
            plot.subtitle = element_text(size = 7.4),
            plot.caption = element_text(size = 7, colour = GREY)))
ggsave(out, p, width = 6.4, height = 6.9, dpi = 170, bg = "#FFFFFF")
cat("wrote", out, "\n")

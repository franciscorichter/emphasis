# mfc_table.R -- the ESS comparison as a LaTeX table, from dev/mfc_ess.csv.
args <- commandArgs(trailingOnly = TRUE)
IN  <- if (length(args) >= 1) args[1] else "dev/mfc_ess.csv"
d <- read.csv(IN, stringsAsFactors = FALSE)
med <- function(x) stats::median(x, na.rm = TRUE)
cols <- c(ess_thin_ed = "thinning, ED-aware", ess_thin_mf = "thinning, mean-field",
          ess_bdi = "conditional, mean-field ED", ess_bdi_fold = "conditional, ED folded into the intercept")
a <- do.call(rbind, lapply(split(d, d$n_target), function(s) data.frame(
  n = s$n_target[1], tips = med(s$n_tips), cells = nrow(s),
  thin_ed = med(s$ess_thin_ed), thin_mf = med(s$ess_thin_mf),
  bdi = med(s$ess_bdi), bdi_fold = med(s$ess_bdi_fold),
  sec_thin = med(s$sec_thin_ed), sec_bdi = med(s$sec_bdi))))
cat("\\begin{center}\\small\n\\begin{tabular}{rrrrrrrr}\n\\toprule\n")
cat("& & \\multicolumn{2}{c}{thinning} & \\multicolumn{2}{c}{conditional} & \\multicolumn{2}{c}{seconds} \\\\\n")
cat("\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\\cmidrule(lr){7-8}\n")
cat("target & tips & ED-aware & mean-field & mean-field $\\overline{\\ED}$ & folded & thinning & cond. \\\\\n\\midrule\n")
for (i in seq_len(nrow(a))) with(a[i, ], cat(sprintf(
  "%d & %.0f & %.1f & %.1f & \\textbf{%.1f} & %.1f & %.1f & %.1f \\\\\n",
  n, tips, thin_ed, thin_mf, bdi, bdi_fold, sec_thin, sec_bdi)))
cat("\\bottomrule\n\\end{tabular}\n\\end{center}\n")
cat("\n% gain of the mean-field conditional over the better thinning arm:\n")
for (i in seq_len(nrow(a))) with(a[i, ], cat(sprintf("%%   %d tips: %.1fx\n", n, bdi / max(thin_ed, thin_mf))))

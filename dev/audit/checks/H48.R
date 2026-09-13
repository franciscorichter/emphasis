## H48: augmented lineages with parent_id == -1 get tip_start = 0 (speciation node)
## while their extinction node keeps tip_start = t_spec; P(t) over-counts by t_spec.
## Self-contained; ~10 s.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages(library(emphasis))

T_EXT_TIP <- 1e11; T_EXT_UNS <- 5e10; T_EXT_EXT <- 0.0
is_ext  <- function(df) df$t_ext == T_EXT_EXT
is_miss <- function(df) !(df$t_ext %in% c(T_EXT_TIP, T_EXT_UNS, T_EXT_EXT))

# crown age 4, 7 tips; first non-crown branching at forward time 4-3.2 = 0.8
brts <- c(4, 3.2, 2.5, 1.8, 1.1, 0.5)
pars8 <- c(0.2, 0, 0, 0, 0.5, 0, 0, 0)     # cr: lambda = 0.2, mu = 0.5 (linear link)

raw <- emphasis:::augment_trees(brts, pars8, sample_size = 400L, maxN = 200000L,
                                max_missing = 1L, max_lambda = 500, num_threads = 1L,
                                model = c(0L, 0L, 0L), link = 0L, rho = 1.0)
trees <- raw$trees
cat("draws returned:", length(trees), "\n")

# hand computation of pendant PD at time tm with a chosen per-lineage start
hand_P <- function(df, tm, ts) {
  alive <- !is_ext(df) & df$brts <= tm & df$t_ext > tm
  sum(tm - ts[alive])
}

# classify draws
n_miss <- sapply(trees, function(d) sum(is_miss(d)))
pid    <- sapply(trees, function(d) { m <- which(is_miss(d)); if (length(m)) d$parent_id[m[1]] else NA })
tspec  <- sapply(trees, function(d) { m <- which(is_miss(d)); if (length(m)) d$brts[m[1]] else NA })
cat("draws with 1 missing lineage:", sum(n_miss == 1), "\n")
cat("  of which parent_id == -1:", sum(n_miss == 1 & pid == -1), "\n")
cat("  t_spec range for parent_id == -1:", range(tspec[n_miss == 1 & pid == -1]), " (first obs branching at 0.8)\n")
cat("  t_spec range for parent_id >= 0 :", range(tspec[n_miss == 1 & pid >= 0]), "\n")

check_draw <- function(df) {
  m  <- which(is_miss(df)); e <- which(is_ext(df))
  ts_code <- df$tip_start                                   # what C++ stored
  ts_true <- ifelse(is_miss(df), df$brts, df$tip_start)     # augmented lineage born at t_spec
  P_code  <- sapply(df$brts, function(tm) hand_P(df, tm, ts_code))
  P_true  <- sapply(df$brts, function(tm) hand_P(df, tm, ts_true))
  c(t_spec          = df$brts[m],
    ts_spec_node    = df$tip_start[m],
    ts_ext_node     = df$tip_start[e],
    E_ext_node      = df$brts[e] - df$tip_start[e],          # e_s for extinction node
    max_abs_pd_vs_code_formula = max(abs(df$pd - P_code)),   # reproduces C++ pd?
    max_abs_pd_vs_true         = max(abs(df$pd - P_true)),   # over-count
    n_nodes_alive_over        = sum(abs(df$pd - P_true) > 1e-12))
}

sel_m1 <- which(n_miss == 1 & pid == -1)
sel_p  <- which(n_miss == 1 & pid >= 0)
res_m1 <- t(sapply(trees[sel_m1], check_draw))
res_p  <- t(sapply(trees[sel_p],  check_draw))

cat("\n== parent_id == -1 draws (", nrow(res_m1), ") ==\n")
cat("tip_start on speciation node == 0 for all:", all(res_m1[, "ts_spec_node"] == 0), "\n")
cat("tip_start on extinction node == t_spec for all:", all(abs(res_m1[, "ts_ext_node"] - res_m1[, "t_spec"]) < 1e-12), "\n")
cat("pd reproduced by the code formula (ts=0) max|diff|:", max(res_m1[, "max_abs_pd_vs_code_formula"]), "\n")
cat("pd vs true pendant (ts=t_spec): max|diff| - t_spec, over draws:",
    range(res_m1[, "max_abs_pd_vs_true"] - res_m1[, "t_spec"]), "\n")
cat("nodes over-counted per draw (min/median/max):",
    min(res_m1[, "n_nodes_alive_over"]), median(res_m1[, "n_nodes_alive_over"]), max(res_m1[, "n_nodes_alive_over"]), "\n")
print(head(round(res_m1, 4), 5))

cat("\n== parent_id >= 0 draws (", nrow(res_p), ") ==\n")
cat("tip_start on speciation node == t_spec for all:", all(abs(res_p[, "ts_spec_node"] - res_p[, "t_spec"]) < 1e-12), "\n")
cat("pd vs true pendant max|diff|:", max(res_p[, "max_abs_pd_vs_true"]), "\n")

## Effect on logf / logg: score a parent_id == -1 draw as returned vs a corrected copy
## (tip_start = t_spec on the speciation node, pd recomputed with that start).
fix_draw <- function(df) {
  m <- which(is_miss(df))
  df$tip_start[m] <- df$brts[m]
  df$pd <- sapply(df$brts, function(tm) hand_P(df, tm, df$tip_start))
  df
}
score <- function(df, pars, model, link) {
  r <- emphasis:::eval_logf(pars, list(df), model = as.integer(model), link = as.integer(link), rho = 1.0)
  c(logf = r$logf, logg = r$logg)
}
cases <- list(
  cr_linear   = list(pars = c(0.2, 0, 0, 0, 0.5, 0, 0, 0),          model = c(0,0,0), link = 0),
  dd_linear   = list(pars = c(0.5, -0.02, 0, 0, 0.3, 0, 0, 0),      model = c(1,0,0), link = 0),
  cr_exp      = list(pars = c(log(0.2), 0, 0, 0, log(0.5), 0, 0, 0), model = c(0,0,0), link = 1),
  m_linear    = list(pars = c(0.3, 0, 0.1, 0, 0.3, 0, 0.05, 0),     model = c(0,1,0), link = 0),
  nd_linear   = list(pars = c(0.3, -0.01, 0, 0.1, 0.3, 0, 0, 0.05), model = c(1,0,1), link = 0),
  nd_exp      = list(pars = c(log(0.3), -0.01, 0, 0.2, log(0.3), 0, 0, 0.1), model = c(1,0,1), link = 1),
  nd_gauss    = list(pars = c(0.3, 0.1, 0, 0.3, 0.3, 0.1, 0, 0.2),  model = c(1,0,1), link = 2)
)
cat("\n== eval_logf on parent_id == -1 draws: as returned minus corrected (max |diff| over", length(sel_m1), "draws) ==\n")
for (nm in names(cases)) {
  cs <- cases[[nm]]
  d <- t(sapply(trees[sel_m1], function(df) {
    a <- score(df, cs$pars, cs$model, cs$link); b <- score(fix_draw(df), cs$pars, cs$model, cs$link)
    a - b
  }))
  cat(sprintf("%-10s max|dlogf| = %.4g   max|dlogg| = %.4g\n", nm, max(abs(d[, "logf"])), max(abs(d[, "logg"]))))
}

## Prevalence with a realistic max_missing (how often does parent_id == -1 occur?)
raw2 <- emphasis:::augment_trees(brts, pars8, sample_size = 200L, maxN = 200000L,
                                 max_missing = 1000L, max_lambda = 500, num_threads = 1L,
                                 model = c(0L, 0L, 0L), link = 0L, rho = 1.0)
all_miss <- do.call(rbind, lapply(raw2$trees, function(d) d[is_miss(d), c("brts", "parent_id", "tip_start")]))
cat("\n== max_missing = 1000:", nrow(all_miss), "augmented lineages in", length(raw2$trees), "draws;",
    sprintf("%.1f%%", 100 * mean(all_miss$parent_id == -1)), "have parent_id == -1;",
    "all of those have tip_start == 0:", all(all_miss$tip_start[all_miss$parent_id == -1] == 0),
    "; all others tip_start == brts:", all(abs(all_miss$tip_start - all_miss$brts)[all_miss$parent_id >= 0] < 1e-12), "\n")
cat("t_spec of parent_id == -1 lineages: all < 0.8 (first obs branching)?",
    all(all_miss$brts[all_miss$parent_id == -1] < 0.8), "; max =", max(all_miss$brts[all_miss$parent_id == -1]), "\n")
cat("parent_id == -1 lineages born before 0.8 / all lineages born before 0.8:",
    sum(all_miss$parent_id == -1), "/", sum(all_miss$brts < 0.8), "\n")

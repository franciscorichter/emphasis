# =============================================================================
# MCEM simulation study — parameter recovery
#
# For each scenario we:
#   1. Simulate N_SIM phylogenetic trees forward in time.
#   2. Estimate parameters with estimate_rates(method = "mcem").
#   3. Compare estimates to the known generative parameters.
#
# The general diversification model is:
#   lambda(t) = max(0, beta_0  + beta_N*N  + beta_P*P)
#   mu(t)     = max(0, gamma_0 + gamma_N*N + gamma_P*P)
#   where P = pendant PD (sum of pendant edge lengths of alive lineages)
#
# Parameters: c(beta_0, [beta_N], [beta_P], gamma_0, [gamma_N], [gamma_P])
#
# Three scenarios
# ┌──────┬────────────────────────────────────────────────────┬──────────────┐
# │  CR  │ constant rate — 2 params (beta_0, gamma_0)         │ self-consistent│
# │  DD  │ diversity-dependent — 4 params                     │ self-consistent│
# │  PD  │ PD-dependent — 4 params                            │ self-consistent│
# └──────┴────────────────────────────────────────────────────┴──────────────┘
#
# Output (OUTDIR):
#   01_recovery_boxplot.pdf     — estimates vs true value per parameter
#   02_scatter_true_vs_est.pdf  — scatter estimated vs true
#   03_bias_rmse.pdf            — bias and RMSE bar charts
#   04_mcem_traces.pdf          — log-likelihood convergence traces
#   05_par_traces.pdf           — parameter traces over MCEM iterations
#   06_tree_sizes.pdf           — tree size distribution per scenario
#   07_error_vs_ntips.pdf       — absolute error vs number of tips
#   summary_table.csv           — bias / RMSE / SD summary
#   results.rds                 — raw list of all fitted results
# =============================================================================


# ── Dependencies ─────────────────────────────────────────────────────────────

pkgs <- c("ggplot2", "patchwork")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

library(emphasis)
library(ggplot2)
library(patchwork)


# ── Global settings (edit here) ───────────────────────────────────────────────

set.seed(42)

N_SIM   <- 30      # replicates per scenario (increase for a full study)
MAX_T   <- 8       # crown age for forward simulation
OUTDIR  <- file.path("inst", "scripts", "sim_study_output")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# MCEM tuning — deliberately fast for a sweep; tighten for final results
MCEM_CTRL <- list(
  sample_size = 100,
  tol         = 0.5,   # SE(fhat) threshold for convergence
  burnin      = 5,
  max_missing = 1e4,
  max_lambda  = 500,
  num_threads = 1L
)


# ── Scenario definitions ──────────────────────────────────────────────────────
#
# Parameters for simulate_tree() and estimate_rates() use the same compact
# ordering: c(beta_0, [beta_N/beta_P], gamma_0, [gamma_N/gamma_P])
# Length = 2 + 2*sum(model_bin).

scenarios <- list(

  CR = list(
    label       = "CR (constant rate)",
    model       = "cr",
    sim_pars    = c(0.5, 0.1),
    true_pars   = c(beta_0 = 0.50, gamma_0 = 0.10),
    lower_bound = c(0.10, 0.01),
    upper_bound = c(1.50, 0.50)
  ),

  DD = list(
    label       = "DD (diversity-dependent)",
    model       = "dd",
    sim_pars    = c(0.8, -0.010, 0.15, 0.0),
    true_pars   = c(beta_0 = 0.80, beta_N = -0.010, gamma_0 = 0.15, gamma_N = 0.000),
    lower_bound = c(0.20, -0.100, 0.01, -0.050),
    upper_bound = c(2.00,  0.010, 0.50,  0.050)
  ),

  PD = list(
    label       = "PD (PD-dependent)",
    model       = "pd",
    sim_pars    = c(0.8, -0.015, 0.15, 0.0),
    true_pars   = c(beta_0 = 0.80, beta_P = -0.015, gamma_0 = 0.15, gamma_P = 0.000),
    lower_bound = c(0.20, -0.100, 0.01, -0.050),
    upper_bound = c(2.00,  0.050, 0.50,  0.050)
  )
)


# ── Simulation + estimation ───────────────────────────────────────────────────

run_scenario <- function(sc, n_sim, max_t, ctrl) {
  cat(sprintf("\n── %s  (%d replicates) ──\n", sc$label, n_sim))
  pb <- utils::txtProgressBar(min = 0, max = n_sim, style = 3)
  rows <- vector("list", n_sim)

  for (k in seq_len(n_sim)) {
    # --- forward simulation ---------------------------------------------------
    tree <- tryCatch(
      simulate_tree(pars = sc$sim_pars, max_t = max_t,
                    model = sc$model, max_tries = 300, useDDD = TRUE),
      error = function(e) NULL
    )

    if (is.null(tree) || tree$status != "done" ||
        is.null(tree$brts) || length(tree$brts) < 2L) {
      utils::setTxtProgressBar(pb, k)
      next
    }

    n_tips <- length(tree$brts) + 1L   # splits + 1 = extant tips

    # --- MCEM estimation -------------------------------------------------------
    fit <- tryCatch(
      estimate_rates(tree        = tree,
                     method      = "mcem",
                     model       = sc$model,
                     lower_bound = sc$lower_bound,
                     upper_bound = sc$upper_bound,
                     control     = ctrl),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      utils::setTxtProgressBar(pb, k)
      next
    }

    # fit$details is the .mcem_dynamic_fresh() return:
    #   $mcem       — data.frame(par1, par2, par3, par4, fhat, sample_size, time)
    #   $pars       — final parameter vector
    #   $iterations — number of iterations
    #   $se         — final SE(fhat)
    mcem_df <- fit$details$mcem

    rows[[k]] <- list(
      rep        = k,
      n_tips     = n_tips,
      est        = fit$pars,           # named numeric(4)
      loglik     = fit$loglik,
      iterations = fit$details$iterations,
      se         = fit$details$se,
      mcem_df    = mcem_df
    )

    utils::setTxtProgressBar(pb, k)
  }
  close(pb)
  Filter(Negate(is.null), rows)
}

all_results <- lapply(scenarios, run_scenario,
                      n_sim = N_SIM, max_t = MAX_T, ctrl = MCEM_CTRL)

cat("\n\nSimulation complete.\n")
n_ok <- sapply(all_results, length)
for (nm in names(n_ok))
  cat(sprintf("  %-4s  %d / %d succeeded\n", nm, n_ok[nm], N_SIM))

saveRDS(list(results = all_results, scenarios = scenarios,
             N_SIM = N_SIM, MAX_T = MAX_T, ctrl = MCEM_CTRL),
        file.path(OUTDIR, "results.rds"))


# ── Build tidy data frames ────────────────────────────────────────────────────

# Collect all parameter names across scenarios
ALL_PAR_NAMES <- unique(unlist(lapply(scenarios, function(sc) names(sc$true_pars))))

# (a) estimates — long: one row per replicate × parameter
est_long <- do.call(rbind, lapply(names(all_results), function(sc_name) {
  sc <- scenarios[[sc_name]]
  pnames <- names(sc$true_pars)
  do.call(rbind, lapply(all_results[[sc_name]], function(r) {
    data.frame(
      scenario  = sc_name,
      label     = sc$label,
      rep       = r$rep,
      n_tips    = r$n_tips,
      parameter = pnames,
      estimate  = as.numeric(r$est[pnames]),
      true_val  = as.numeric(sc$true_pars[pnames]),
      stringsAsFactors = FALSE
    )
  }))
}))
est_long$parameter <- factor(est_long$parameter, levels = ALL_PAR_NAMES)
est_long$error     <- est_long$estimate - est_long$true_val

# (b) MCEM iteration traces — long: one row per replicate × iteration × variable
build_trace <- function(sc_name, par_long = FALSE) {
  sc <- scenarios[[sc_name]]
  do.call(rbind, lapply(all_results[[sc_name]], function(r) {
    df <- r$mcem_df
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df$scenario  <- sc_name
    df$rep       <- r$rep
    df$iteration <- seq_len(nrow(df))
    df
  }))
}

traces_all <- do.call(rbind, lapply(names(scenarios), build_trace))
traces_all <- traces_all[is.finite(traces_all$fhat), , drop = FALSE]

# parameter traces — melt par columns to long, using scenario-specific names
par_trace_long <- do.call(rbind, lapply(names(scenarios), function(sc_name) {
  sc <- scenarios[[sc_name]]
  pnames <- names(sc$true_pars)
  d_sc <- traces_all[traces_all$scenario == sc_name, , drop = FALSE]
  if (nrow(d_sc) == 0) return(NULL)
  do.call(rbind, lapply(seq_along(pnames), function(j) {
    col <- paste0("par", j)
    if (!(col %in% names(d_sc))) return(NULL)
    d <- d_sc[, c("scenario", "rep", "iteration", col), drop = FALSE]
    d$parameter <- pnames[j]
    d$value     <- d[[col]]
    d[[col]]    <- NULL
    d
  }))
}))
par_trace_long$parameter <- factor(par_trace_long$parameter, levels = ALL_PAR_NAMES)

# add true value for reference lines
par_trace_long$true_val <- mapply(
  function(sc, par) scenarios[[sc]]$true_pars[par],
  par_trace_long$scenario, as.character(par_trace_long$parameter)
)

# (c) per-replicate summary
rep_summary <- do.call(rbind, lapply(names(all_results), function(sc_name) {
  do.call(rbind, lapply(all_results[[sc_name]], function(r) {
    data.frame(scenario   = sc_name,
               rep        = r$rep,
               n_tips     = r$n_tips,
               loglik     = r$loglik,
               iterations = r$iterations,
               se         = r$se,
               stringsAsFactors = FALSE)
  }))
}))

# (d) bias / RMSE summary table
bias_rmse <- do.call(rbind, lapply(split(est_long,
    list(est_long$scenario, est_long$parameter)), function(d) {
  data.frame(
    scenario  = d$scenario[1],
    parameter = as.character(d$parameter[1]),
    true_val  = d$true_val[1],
    n         = nrow(d),
    mean_est  = round(mean(d$estimate, na.rm = TRUE), 5),
    sd_est    = round(sd(d$estimate,   na.rm = TRUE), 5),
    bias      = round(mean(d$error,    na.rm = TRUE), 5),
    rmse      = round(sqrt(mean(d$error^2, na.rm = TRUE)), 5),
    stringsAsFactors = FALSE
  )
}))
bias_rmse$parameter <- factor(bias_rmse$parameter, levels = ALL_PAR_NAMES)
bias_rmse <- bias_rmse[order(bias_rmse$scenario, bias_rmse$parameter), ]
row.names(bias_rmse) <- NULL

cat("\n── Summary table ─────────────────────────────────────────────────────────\n")
print(bias_rmse, row.names = FALSE)
write.csv(bias_rmse, file.path(OUTDIR, "summary_table.csv"), row.names = FALSE)


# ── Shared helpers ────────────────────────────────────────────────────────────

# True-value data frame for geom_point overlays
true_df <- do.call(rbind, lapply(names(scenarios), function(sc_name) {
  sc <- scenarios[[sc_name]]
  pnames <- names(sc$true_pars)
  data.frame(scenario  = sc_name,
             parameter = factor(pnames, levels = ALL_PAR_NAMES),
             true_val  = as.numeric(sc$true_pars),
             stringsAsFactors = FALSE)
}))

scenario_colours <- setNames(
  RColorBrewer::brewer.pal(max(3, length(scenarios)), "Set2")[seq_along(scenarios)],
  names(scenarios)
)
# Fallback if RColorBrewer not available
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  scenario_colours <- setNames(c("#66C2A5","#FC8D62","#8DA0CB"), names(scenarios))
}

th <- theme_minimal(base_size = 13) +
  theme(strip.text      = element_text(face = "bold"),
        legend.position = "bottom")

subtitle_n <- sprintf("%d replicates per scenario, crown age = %d", N_SIM, MAX_T)


# ── Plot 1: Recovery boxplots ─────────────────────────────────────────────────

p1 <- ggplot(est_long, aes(x = scenario, y = estimate, fill = scenario)) +
  geom_boxplot(alpha = 0.75, outlier.size = 1.2, outlier.alpha = 0.5) +
  geom_point(data = true_df,
             aes(x = scenario, y = true_val),
             shape = 23, size = 3.5, fill = "red", colour = "black",
             inherit.aes = FALSE) +
  facet_wrap(~ parameter, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = scenario_colours, guide = "none") +
  labs(title    = "MCEM parameter recovery — boxplots",
       subtitle = subtitle_n,
       x = "Scenario", y = "Estimated value",
       caption  = "Red diamond = true generative value") +
  th

ggsave(file.path(OUTDIR, "01_recovery_boxplot.pdf"), p1, width = 10, height = 9)
cat("Saved 01_recovery_boxplot.pdf\n")


# ── Plot 2: Scatter estimated vs true ─────────────────────────────────────────

p2 <- ggplot(est_long, aes(x = true_val, y = estimate, colour = scenario)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey50", linewidth = 0.5) +
  geom_point(alpha = 0.65, size = 2) +
  facet_wrap(~ parameter, scales = "free", ncol = 2) +
  scale_colour_manual(values = scenario_colours) +
  labs(title    = "Estimated vs true parameter values",
       subtitle = "Diagonal = perfect recovery",
       x = "True value", y = "Estimated value", colour = "Scenario") +
  th

ggsave(file.path(OUTDIR, "02_scatter_true_vs_est.pdf"), p2, width = 10, height = 9)
cat("Saved 02_scatter_true_vs_est.pdf\n")


# ── Plot 3: Bias and RMSE bar charts ──────────────────────────────────────────

p3a <- ggplot(bias_rmse, aes(x = parameter, y = bias, fill = scenario)) +
  geom_col(position = position_dodge(width = 0.75), colour = "white", width = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = scenario_colours) +
  labs(title = "Bias (mean error)", x = "Parameter",
       y = "mean(est − true)", fill = NULL) +
  th + theme(legend.position = "right")

p3b <- ggplot(bias_rmse, aes(x = parameter, y = rmse, fill = scenario)) +
  geom_col(position = position_dodge(width = 0.75), colour = "white", width = 0.7) +
  scale_fill_manual(values = scenario_colours) +
  labs(title = "RMSE", x = "Parameter", y = "RMSE", fill = NULL) +
  th + theme(legend.position = "right")

ggsave(file.path(OUTDIR, "03_bias_rmse.pdf"),
       p3a / p3b, width = 10, height = 9)
cat("Saved 03_bias_rmse.pdf\n")


# ── Plot 4: MCEM log-likelihood convergence traces ────────────────────────────

# Up to 8 replicates per scenario for clarity
set.seed(1)
trace_sub <- do.call(rbind, lapply(names(scenarios), function(sc_name) {
  d    <- traces_all[traces_all$scenario == sc_name, ]
  reps <- unique(d$rep)
  sel  <- if (length(reps) > 8) sample(reps, 8) else reps
  d[d$rep %in% sel, ]
}))
trace_sub$rep_id <- paste(trace_sub$scenario, trace_sub$rep)

p4 <- ggplot(trace_sub,
             aes(x = iteration, y = fhat, group = rep_id, colour = scenario)) +
  geom_line(alpha = 0.55, linewidth = 0.5) +
  geom_point(data = trace_sub[!duplicated(paste(trace_sub$rep_id,
                                                 trace_sub$iteration == max(trace_sub$iteration))), ],
             aes(x = iteration, y = fhat), size = 0.8, alpha = 0.4) +
  facet_wrap(~ scenario, scales = "free_y", ncol = 1) +
  scale_colour_manual(values = scenario_colours, guide = "none") +
  labs(title    = "MCEM log-likelihood traces",
       subtitle = "Up to 8 replicates per scenario",
       x = "MCEM iteration",
       y = expression(hat(italic(f))[italic(t)])) +
  th + theme(legend.position = "none")

ggsave(file.path(OUTDIR, "04_mcem_traces.pdf"), p4, width = 8, height = 11)
cat("Saved 04_mcem_traces.pdf\n")


# ── Plot 5: Parameter convergence traces ──────────────────────────────────────

# One scenario at a time; use same 8-rep subset from above
par_trace_sub <- par_trace_long[
  paste(par_trace_long$scenario, par_trace_long$rep) %in% trace_sub$rep_id, ]
par_trace_sub$rep_id <- paste(par_trace_sub$scenario, par_trace_sub$rep)

p5_list <- lapply(names(scenarios), function(sc_name) {
  d <- par_trace_sub[par_trace_sub$scenario == sc_name, ]
  if (nrow(d) == 0) return(NULL)
  ggplot(d, aes(x = iteration, y = value, group = rep_id)) +
    geom_line(alpha = 0.5, linewidth = 0.4,
              colour = scenario_colours[[sc_name]]) +
    geom_hline(aes(yintercept = true_val), linetype = "dashed",
               colour = "red", linewidth = 0.6) +
    facet_wrap(~ parameter, scales = "free_y", ncol = 2) +
    labs(title = scenarios[[sc_name]]$label,
         x = "MCEM iteration", y = "Parameter value") +
    th + theme(strip.text = element_text(size = 9))
})
p5_list <- Filter(Negate(is.null), p5_list)

p5 <- patchwork::wrap_plots(p5_list, ncol = 1) +
  patchwork::plot_annotation(
    title   = "Parameter traces over MCEM iterations",
    caption = "Dashed red line = true generative value"
  )

ggsave(file.path(OUTDIR, "05_par_traces.pdf"),
       p5, width = 10, height = 5 * length(p5_list))
cat("Saved 05_par_traces.pdf\n")


# ── Plot 6: Tree size distribution ────────────────────────────────────────────

ntips_df <- unique(rep_summary[, c("scenario", "rep", "n_tips")])

p6 <- ggplot(ntips_df, aes(x = n_tips, fill = scenario)) +
  geom_histogram(binwidth = 3, colour = "white", alpha = 0.85) +
  facet_wrap(~ scenario, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = scenario_colours, guide = "none") +
  labs(title    = "Distribution of tree size across replicates",
       subtitle = sprintf("Crown age = %d, %d replicates", MAX_T, N_SIM),
       x = "Extant tips", y = "Count") +
  th + theme(legend.position = "none")

ggsave(file.path(OUTDIR, "06_tree_sizes.pdf"), p6, width = 7, height = 8)
cat("Saved 06_tree_sizes.pdf\n")


# ── Plot 7: Absolute error vs tree size ───────────────────────────────────────

p7 <- ggplot(est_long, aes(x = n_tips, y = abs(error), colour = scenario)) +
  geom_point(alpha = 0.35, size = 1.5) +
  geom_smooth(method = "loess", formula = y ~ x,
              se = FALSE, linewidth = 0.9) +
  facet_wrap(~ parameter, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = scenario_colours) +
  labs(title    = "Estimation error vs tree size",
       subtitle = "Loess smoother; smaller trees → larger errors expected",
       x = "Extant tips", y = "|Estimated − true|", colour = "Scenario") +
  th

ggsave(file.path(OUTDIR, "07_error_vs_ntips.pdf"), p7, width = 10, height = 9)
cat("Saved 07_error_vs_ntips.pdf\n")


# ── MCEM diagnostics plot ────────────────────────────────────────────────────
# Iterations to convergence and final SE(fhat)

p8a <- ggplot(rep_summary, aes(x = scenario, y = iterations, fill = scenario)) +
  geom_boxplot(alpha = 0.75, outlier.size = 1.2) +
  scale_fill_manual(values = scenario_colours, guide = "none") +
  labs(title = "MCEM iterations to convergence", x = NULL, y = "Iterations") +
  th

p8b <- ggplot(rep_summary, aes(x = scenario, y = se, fill = scenario)) +
  geom_boxplot(alpha = 0.75, outlier.size = 1.2) +
  geom_hline(yintercept = MCEM_CTRL$tol,
             linetype = "dashed", colour = "grey30", linewidth = 0.5) +
  annotate("text", x = 0.5, y = MCEM_CTRL$tol * 1.05,
           label = paste("tol =", MCEM_CTRL$tol),
           hjust = 0, size = 3.5, colour = "grey30") +
  scale_fill_manual(values = scenario_colours, guide = "none") +
  labs(title = "Final SE(fhat) at convergence", x = NULL, y = "SE(fhat)") +
  th

ggsave(file.path(OUTDIR, "08_mcem_diagnostics.pdf"),
       p8a / p8b, width = 8, height = 8)
cat("Saved 08_mcem_diagnostics.pdf\n")


# ── Print final summary ───────────────────────────────────────────────────────

cat("\n\n═══════════════════════════════════════════════════════════════\n")
cat("  Simulation study complete\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat(sprintf("  Scenarios:  %s\n", paste(names(scenarios), collapse = ", ")))
cat(sprintf("  Replicates: %d per scenario (%d requested)\n",
            sum(n_ok), N_SIM * length(scenarios)))
cat(sprintf("  Crown age:  %d\n", MAX_T))
cat(sprintf("  Output:     %s\n\n", normalizePath(OUTDIR)))

cat("  Bias / RMSE summary:\n\n")
print(bias_rmse[, c("scenario","parameter","true_val","mean_est","bias","rmse")],
      row.names = FALSE)

cat("\n  Output files:\n")
for (f in sort(list.files(OUTDIR))) cat(sprintf("    • %s\n", f))

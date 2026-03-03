#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(emphasis)
  library(tibble)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(purrr)
  library(tidyr)
  library(patchwork)
})

set.seed(20250301)

n_iter <- 30
max_t_values <- c(3, 6, 9)
max_lin <- 5e4
max_tries <- 50
inner_reps_fast <- 8
inner_reps_slow <- 8

out_dir <- file.path("results", "benchmarks")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

combo_defs <- tibble::tribble(
  ~combo_id, ~combo_label, ~model, ~max_t, ~pars,
  "CR-weak",  "CR (low rates): mu=0.12, lambda0=0.40",          "pd", 4, list(c(0.12, 0.40, 0.00, 0.00)),
  "CR-strong","CR (high rates): mu=0.08, lambda0=0.60",         "pd", 5, list(c(0.08, 0.60, 0.00, 0.00)),
  "DD-light", "DD-lite: betaN=-0.015, betaP=0",                 "pd", 5, list(c(0.10, 0.55, -0.015, 0.00)),
  "DD-strong","DD-strong: betaN=-0.035, betaP=0",               "pd", 6, list(c(0.08, 0.60, -0.035, 0.00)),
  "PD-mixed", "PD mixed: betaN=-0.025, betaP=0.015",            "pd", 6, list(c(0.08, 0.55, -0.025, 0.015)),
  "PD-heavy", "PD heavy: betaN=-0.040, betaP=0.030",            "pd", 7, list(c(0.06, 0.65, -0.040, 0.030))
)

combo_defs$combo_label <- factor(combo_defs$combo_label, levels = combo_defs$combo_label)

scenarios <- combo_defs |>
  tidyr::crossing(fast = c(TRUE, FALSE))

extract_tree_metrics <- function(res, fast_backend) {
  if (inherits(res, "try-error") || is.null(res)) {
    return(c(tree_size = 0, tip_count = 0, pd_total = 0))
  }

  tas <- res$tas
  tes <- res$tes
  L <- res$L

  if (!is.null(tas)) {
    tree_size <- nrow(tas$edge)
    tip_count <- if (!is.null(tes)) length(tes$tip.label) else length(tas$tip.label)
    pd_total <- sum(tas$edge.length)
  } else if (!is.null(L)) {
    tree_size <- nrow(L)
    tip_count <- sum(L[, 4] == -1)
    pd_total <- sum(L[, 3] - L[, 2])
  } else {
    tree_size <- 0
    tip_count <- 0
    pd_total <- 0
  }

  c(tree_size = tree_size, tip_count = tip_count, pd_total = pd_total)
}

run_single <- function(pars, max_t, model, fast) {
  reps <- if (fast) inner_reps_fast else inner_reps_slow
  run_records <- vector("list", reps)
  elapsed_each <- numeric(reps)

  for (j in seq_len(reps)) {
    start_time <- Sys.time()
    res <- try(suppressWarnings(
      simulate_tree(
        pars = pars,
        max_t = max_t,
        model = model,
        fast = fast,
        max_lin = max_lin,
        max_tries = max_tries,
        useDDD = TRUE
      )
    ))
    elapsed_each[j] <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    metrics <- extract_tree_metrics(res, fast)
    run_records[[j]] <- list(
      status = if (inherits(res, "try-error")) "error" else res$status,
      tree_size = metrics[["tree_size"]],
      tip_count = metrics[["tip_count"]],
      pd_total = metrics[["pd_total"]]
    )
  }

  tibble(
    inner_iter = seq_len(reps),
    status = vapply(run_records, `[[`, character(1), "status"),
    tree_size = vapply(run_records, `[[`, numeric(1), "tree_size"),
    tip_count = vapply(run_records, `[[`, numeric(1), "tip_count"),
    pd_total = vapply(run_records, `[[`, numeric(1), "pd_total"),
    elapsed = elapsed_each,
    elapsed_total = sum(elapsed_each)
  )
}

results <- pmap_dfr(
  scenarios,
  function(combo_id, combo_label, model, max_t, pars, fast) {
    message(sprintf("Running %s (%s) / fast=%s / max_t=%s", combo_id, model, fast, max_t))
    pars_vec <- unlist(pars)
    map_dfr(seq_len(n_iter), function(i) {
      run_single(pars_vec, max_t, model, fast) |>
        mutate(
          replicate = i,
          combo_id = combo_id,
          combo_label = combo_label,
          model = model,
          fast = fast,
          max_t = max_t,
          pars = list(pars_vec)
        )
    })
  }
)

csv_path <- file.path(out_dir, "simulate_tree_speed.csv")
write_csv(results, csv_path)
message("Saved raw timings to ", csv_path)

plot_data <- results |>
  mutate(
    backend = ifelse(fast, "C++ (fast)", "R (slow)"),
    elapsed_adj = pmax(elapsed, 1e-4),
    tree_size = tidyr::replace_na(tree_size, 0),
    tip_count = tidyr::replace_na(tip_count, 0),
    pd_total = tidyr::replace_na(pd_total, 0)
  )

done_only <- plot_data |> filter(status == "done")

status_summary <- plot_data |>
  count(combo_label, backend, status) |>
  group_by(combo_label, backend) |>
  mutate(prop = n / sum(n)) |>
  ungroup()

extinct_summary <- status_summary |>
  filter(status == "extinct") |>
  select(combo_label, backend, extinct_prop = prop) |>
  tidyr::complete(combo_label, backend, fill = list(extinct_prop = 0))

summary_table <- done_only |>
  group_by(combo_label, backend) |>
  summarise(
    median_elapsed = median(elapsed, na.rm = TRUE),
    mean_elapsed = mean(elapsed, na.rm = TRUE),
    median_tree_size = median(tree_size, na.rm = TRUE),
    mean_tip_count = mean(tip_count, na.rm = TRUE),
    mean_pd_total = mean(pd_total, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(
    status_summary |>
      group_by(combo_label, backend) |>
      summarise(done_rate = prop[status == "done"], .groups = "drop"),
    by = c("combo_label", "backend")
  )

summary_path <- file.path(out_dir, "simulate_tree_speed_summary.csv")
readr::write_csv(summary_table, summary_path)
message("Saved summary metrics to ", summary_path)

palette_backend <- c("C++ (fast)" = "#1b9e77", "R (slow)" = "#d95f02")

box_plot <- ggplot(plot_data, aes(x = backend, y = elapsed_adj, fill = backend)) +
  geom_boxplot(alpha = 0.85, width = 0.6, outlier.size = 0.6, color = "#222") +
  geom_point(aes(color = backend), position = position_jitter(width = 0.08),
             size = 0.5, alpha = 0.4, show.legend = FALSE) +
  facet_wrap(~combo_label, ncol = 3) +
  scale_y_log10(labels = scales::label_number(accuracy = 0.001)) +
  scale_fill_manual(values = palette_backend) +
  scale_color_manual(values = palette_backend) +
  labs(title = "Elapsed time across scenarios",
       subtitle = "Log-scale boxplots with jittered runs",
       x = NULL, y = "Elapsed time (s, log10)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "none")

size_plot <- ggplot(plot_data, aes(x = tree_size, y = elapsed_adj, color = backend)) +
  geom_point(alpha = 0.65, size = 0.8) +
  geom_smooth(se = FALSE, method = "loess", span = 0.75, linewidth = 0.7) +
  facet_wrap(~combo_label, ncol = 3, scales = "free_x") +
  scale_y_log10(labels = scales::label_number(accuracy = 0.001)) +
  scale_color_manual(values = palette_backend) +
  labs(title = "Tree size vs elapsed time",
       subtitle = "Points are successful runs; y-axis log scale",
       x = "Tree size (rows / edges)", y = "Elapsed time (s, log10)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

size_density_plot <- ggplot(plot_data, aes(x = tree_size, fill = backend)) +
  geom_density(alpha = 0.5, color = NA, adjust = 1.1) +
  facet_wrap(~combo_label, ncol = 3, scales = "free_x") +
  scale_fill_manual(values = palette_backend) +
  labs(title = "Tree size distributions",
       subtitle = "Overlapping densities by backend",
       x = "Tree size (rows / edges)", y = "Density") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

tip_density_plot <- ggplot(plot_data, aes(x = tip_count, fill = backend)) +
  geom_density(alpha = 0.5, color = NA, adjust = 1.1) +
  facet_wrap(~combo_label, ncol = 3, scales = "free_x") +
  scale_fill_manual(values = palette_backend) +
  geom_text(
    data = extinct_summary,
    aes(x = -Inf, y = Inf,
        label = paste0("extinct: ", scales::percent(extinct_prop, accuracy = 0.1)),
        color = backend),
    inherit.aes = FALSE,
    hjust = -0.05, vjust = 1.2, size = 3.1, show.legend = FALSE
  ) +
  scale_color_manual(values = palette_backend, guide = "none") +
  labs(title = "Species count distributions",
       subtitle = "Overlapping densities by backend",
       x = "Number of extant species", y = "Density") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

tip_pd_plot <- ggplot(done_only, aes(x = tip_count, y = pd_total, color = backend)) +
  geom_point(alpha = 0.65, size = 0.9) +
  geom_smooth(se = FALSE, method = "loess", linewidth = 0.8, span = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#999999") +
  facet_wrap(~combo_label, ncol = 3, scales = "free") +
  scale_color_manual(values = palette_backend) +
  labs(title = "Species vs phylogenetic diversity",
       subtitle = "Each point is one simulated tree (successful runs)",
       x = "Number of extant species", y = "Total phylogenetic diversity") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

combined_plot <- (tip_pd_plot / size_plot / size_density_plot / tip_density_plot) +
  patchwork::plot_layout(heights = c(1, 1, 0.7, 0.7)) +
  patchwork::plot_annotation(
    title = "simulate_tree() fast vs slow benchmarks",
    subtitle = sprintf("Each facet is a parameter combination (n_iter = %d, inner reps = %d)",
                      n_iter, inner_reps_fast)
  )

ggsave(file.path(out_dir, "simulate_tree_comparison_grid.png"), combined_plot,
       width = 12, height = 15, dpi = 220)

ggsave(file.path(out_dir, "simulate_tree_time_vs_size.png"), size_plot,
       width = 12, height = 10, dpi = 220)

ggsave(file.path(out_dir, "simulate_tree_size_density.png"), size_density_plot,
       width = 12, height = 8, dpi = 220)

ggsave(file.path(out_dir, "simulate_tree_tip_density.png"), tip_density_plot,
       width = 12, height = 8, dpi = 220)

message("Saved enhanced plots to ", out_dir)

combo_dir <- file.path(out_dir, "per_combo")
if (!dir.exists(combo_dir)) dir.create(combo_dir, recursive = TRUE)

purrr::walk(unique(plot_data$combo_id), function(cid) {
  combo_label <- plot_data$combo_label[plot_data$combo_id == cid][1]
  df_all <- plot_data |> filter(combo_id == cid)
  combo_extinct <- extinct_summary |> filter(combo_label == !!combo_label)
  scatter_combo <- ggplot(df_all, aes(x = tree_size, y = elapsed_adj, color = backend)) +
    geom_point(alpha = 0.7, size = 1.2) +
    geom_smooth(se = FALSE, method = "loess", span = 0.8, linewidth = 0.8) +
    scale_y_log10(labels = scales::label_number(accuracy = 0.001)) +
    scale_color_manual(values = palette_backend) +
    labs(title = as.character(combo_label),
         subtitle = "Elapsed time vs tree size (log y)",
         x = "Tree size (rows / edges)", y = "Elapsed time (s, log10)") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")

  density_combo <- ggplot(df_all, aes(x = tip_count, fill = backend)) +
    geom_density(alpha = 0.5, color = NA) +
    scale_fill_manual(values = palette_backend) +
    geom_text(
      data = combo_extinct,
      aes(x = -Inf, y = Inf,
          label = paste0("extinct: ", scales::percent(extinct_prop, accuracy = 0.1)),
          color = backend),
      inherit.aes = FALSE,
      hjust = -0.05, vjust = 1.2, size = 3.2, show.legend = FALSE
    ) +
    scale_color_manual(values = palette_backend, guide = "none") +
    labs(title = "Species count density",
         x = "Number of extant species", y = "Density") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")

  combo_panel <- (scatter_combo / density_combo) +
    patchwork::plot_layout(heights = c(1.4, 0.8))

  combo_path <- file.path(combo_dir, sprintf("simulate_tree_%s.png", cid))
  ggsave(combo_path, combo_panel, width = 7.5, height = 10, dpi = 220)
})

message("Saved per-scenario panels to ", combo_dir)

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

n_iter <- 12
max_t_values <- c(3, 6, 9)
max_lin <- 5e4
max_tries <- 50
inner_reps_fast <- 8
inner_reps_slow <- 1

out_dir <- file.path("results", "benchmarks")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

combo_defs <- tibble::tribble(
  ~combo_id, ~combo_label, ~model, ~max_t, ~pars,
  "CR-1", "CR: mu=0.15, lambda=0.35 (max_t=3)", "cr", 3, list(c(0.15, 0.35)),
  "CR-2", "CR: mu=0.10, lambda=0.55 (max_t=6)", "cr", 6, list(c(0.10, 0.55)),
  "CR-3", "CR: mu=0.05, lambda=0.75 (max_t=9)", "cr", 9, list(c(0.05, 0.75)),
  "PD-1", "PD: mu=0.10, lambda0=0.40, betaN=-0.02, betaP=0 (max_t=3)", "pd", 3, list(c(0.10, 0.40, -0.02, 0.00)),
  "PD-2", "PD: mu=0.08, lambda0=0.50, betaN=-0.03, betaP=0.02 (max_t=6)", "pd", 6, list(c(0.08, 0.50, -0.03, 0.02)),
  "PD-3", "PD: mu=0.06, lambda0=0.65, betaN=-0.04, betaP=0.03 (max_t=9)", "pd", 9, list(c(0.06, 0.65, -0.04, 0.03))
)

combo_defs$combo_label <- factor(combo_defs$combo_label, levels = combo_defs$combo_label)

scenarios <- combo_defs |>
  tidyr::crossing(fast = c(TRUE, FALSE)) |>
  dplyr::filter(!(combo_id == "CR-3" & fast == FALSE))

extract_tree_metrics <- function(res, fast_backend) {
  if (inherits(res, "try-error") || is.null(res)) {
    return(c(tree_size = NA_real_, tip_count = NA_real_))
  }

  if (fast_backend) {
    L <- res$L
    if (is.null(L)) {
      return(c(tree_size = NA_real_, tip_count = NA_real_))
    }
    tree_size <- nrow(L)
    tip_count <- sum(L[, 4] == -1)
  } else {
    tas <- res$tas
    if (is.null(tas)) {
      return(c(tree_size = NA_real_, tip_count = NA_real_))
    }
    tree_size <- nrow(tas$edge)
    tip_count <- length(tas$tip.label)
  }

  c(tree_size = tree_size, tip_count = tip_count)
}

run_single <- function(pars, max_t, model, fast) {
  reps <- if (fast) inner_reps_fast else inner_reps_slow
  run_records <- vector("list", reps)
  total_elapsed <- system.time({
    for (j in seq_len(reps)) {
      res <- try(suppressWarnings(
        simulate_tree(
          pars = pars,
          max_t = max_t,
          model = model,
          fast = fast,
          max_lin = max_lin,
          max_tries = max_tries,
          useDDD = FALSE
        )
      ))
      metrics <- extract_tree_metrics(res, fast)
      run_records[[j]] <- list(
        status = if (inherits(res, "try-error")) "error" else res$status,
        tree_size = metrics[["tree_size"]],
        tip_count = metrics[["tip_count"]]
      )
    }
  })[["elapsed"]]

  tibble(
    inner_iter = seq_len(reps),
    status = vapply(run_records, `[[`, character(1), "status"),
    tree_size = vapply(run_records, `[[`, numeric(1), "tree_size"),
    tip_count = vapply(run_records, `[[`, numeric(1), "tip_count"),
    elapsed = total_elapsed / reps,
    elapsed_total = total_elapsed,
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
  filter(status == "done") |>
  mutate(
    elapsed_adj = pmax(elapsed, 1e-4),
    backend = ifelse(fast, "C++ (fast)", "R (slow)")
  )

status_summary <- results |>
  count(combo_label, backend = ifelse(fast, "C++ (fast)", "R (slow)"), status) |>
  group_by(combo_label, backend) |>
  mutate(prop = n / sum(n)) |>
  ungroup()

box_plot <- ggplot(plot_data, aes(x = backend, y = elapsed, fill = backend)) +
  geom_boxplot(alpha = 0.8, width = 0.6, outlier.size = 0.6) +
  facet_grid(combo_label ~ ., scales = "free_y") +
  scale_fill_manual(values = c("C++ (fast)" = "#1b9e77", "R (slow)" = "#d95f02")) +
  labs(x = NULL, y = "Elapsed time (s)", title = "Elapsed time distribution per combo") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "none")

size_plot <- ggplot(plot_data, aes(x = tree_size, y = elapsed_adj, color = backend)) +
  geom_point(alpha = 0.6, size = 0.9) +
  facet_grid(combo_label ~ ., scales = "free_x") +
  scale_y_log10() +
  scale_color_manual(values = c("C++ (fast)" = "#1b9e77", "R (slow)" = "#d95f02")) +
  labs(x = "Tree size (rows in L-table / edges)", y = "Elapsed time (s, log10)",
       title = "Tree size vs elapsed time") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

status_plot <- ggplot(status_summary, aes(x = backend, y = prop, fill = status)) +
  geom_col(position = "stack", color = "white", width = 0.6) +
  facet_grid(combo_label ~ .) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = "% of outcomes", title = "Outcome proportions") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

combined_plot <- patchwork::wrap_plots(
  box_plot,
  size_plot + theme(legend.position = "none"),
  status_plot,
  ncol = 3
)

ggsave(file.path(out_dir, "simulate_tree_comparison_grid.png"), combined_plot,
       width = 12, height = 14, dpi = 200)

ggsave(file.path(out_dir, "simulate_tree_time_vs_size.png"), size_plot,
       width = 6, height = 14, dpi = 200)

message("Saved enhanced plots to ", out_dir)

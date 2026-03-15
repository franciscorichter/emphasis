#!/usr/bin/env Rscript
# ============================================================================
# Simulation study: parameter recovery + model selection
#
# 3 models (DD, PD, EP) x 10 trees each = 30 simulated trees.
# For each tree: estimate parameters under all 3 models, select by AIC.
# Report: confusion matrix, accuracy, parameter recovery boxplots.
#
# Usage:  Rscript simulation_study.R
# Output: console + simulation_study_results.RData
# ============================================================================

# ── Install dependencies ─────────────────────────────────────────
for (pkg in c("devtools", "ape", "mgcv", "ggplot2")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}
if (requireNamespace("emphasis", quietly = TRUE)) {
  try(remove.packages("emphasis"), silent = TRUE)
}
devtools::install_github(
  "franciscorichter/emphasis",
  force = TRUE, upgrade = "never"
)

library(emphasis)
library(ggplot2)

cat("emphasis version:",
    as.character(packageVersion("emphasis")), "\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n\n")

# ── Configuration ────────────────────────────────────────────────

n_trees_per_model <- 10
max_t             <- 10
link              <- "exponential"
n_cores           <- max(1L, parallel::detectCores() - 1L)

true_models <- list(
  DD = list(
    model = "dd",
    pars  = c(0.8, -0.02, 0.2, 0.0)
  ),
  PD = list(
    model = "pd",
    pars  = c(0.6, 0.005, 0.15, -0.003)
  ),
  EP = list(
    model = "ep",
    pars  = c(0.5, 0.05, 0.1, 0.01),
    link  = "linear"   # EP only supports linear link
  )
)

# GAM control for estimation
gam_ctrl <- list(n_grid = 100, sample_size = 50)

# ── 1. Simulate trees ───────────────────────────────────────────

cat("=== Phase 1: Simulating trees ===\n\n")
set.seed(2026)

all_trees    <- list()
tree_info    <- data.frame()

for (mod_name in names(true_models)) {
  cfg      <- true_models[[mod_name]]
  mod_link <- if (!is.null(cfg$link)) cfg$link else link

  cat(sprintf("  [%s] simulating %d trees (link=%s)...\n",
              mod_name, n_trees_per_model, mod_link))

  count <- 0L
  attempts <- 0L
  while (count < n_trees_per_model && attempts < 500) {
    attempts <- attempts + 1L
    sim <- simulate_tree(
      pars = cfg$pars, model = cfg$model,
      max_t = max_t, link = mod_link
    )
    if (sim$status == "done" &&
        length(sim$tes$tip.label) >= 5) {
      count <- count + 1L
      key <- sprintf("%s_%02d", mod_name, count)
      all_trees[[key]] <- sim
      tree_info <- rbind(tree_info, data.frame(
        key       = key,
        true_model = mod_name,
        n_tips    = length(sim$tes$tip.label),
        link      = mod_link,
        stringsAsFactors = FALSE
      ))
      cat(sprintf("    tree %d: %d tips\n",
                  count, length(sim$tes$tip.label)))
    }
  }
}

cat(sprintf("\nTotal trees: %d\n\n", length(all_trees)))

# ── 2. Estimate under all models ────────────────────────────────

cat("=== Phase 2: Parameter estimation ===\n\n")

fit_models <- c("dd", "pd", "ep")
fit_links  <- c(dd = "exponential",
                pd = "exponential",
                ep = "linear")

all_fits    <- list()
fit_results <- data.frame()

for (i in seq_len(nrow(tree_info))) {
  key  <- tree_info$key[i]
  tree <- all_trees[[key]]
  cat(sprintf(
    "  [%d/%d] %s (%d tips)\n",
    i, nrow(tree_info), key,
    tree_info$n_tips[i]
  ))

  for (mod in fit_models) {
    mod_link <- fit_links[[mod]]

    # Auto-detect bounds
    ab <- tryCatch(
      auto_bounds(
        tree, model = mod, link = mod_link,
        n_probe = 200,
        train_surv_gam = FALSE,
        num_threads = n_cores,
        verbose = FALSE
      ),
      error = function(e) NULL
    )

    if (is.null(ab) || ab$n_feasible < 5) {
      cat(sprintf("    %s: bounds failed\n", mod))
      fit_results <- rbind(fit_results, data.frame(
        key = key, true_model = tree_info$true_model[i],
        fit_model = toupper(mod), link = mod_link,
        loglik = NA, AIC = NA,
        stringsAsFactors = FALSE
      ))
      next
    }

    ctrl <- gam_ctrl
    ctrl$lower_bound <- ab$lower_bound
    ctrl$upper_bound <- ab$upper_bound

    fit <- tryCatch(
      estimate_rates(
        tree, method = "gam",
        model = mod, link = mod_link,
        control = ctrl
      ),
      error = function(e) NULL
    )

    if (!is.null(fit)) {
      cat(sprintf(
        "    %s: loglik=%.2f AIC=%.2f\n",
        mod, fit$loglik, fit$AIC
      ))
      fit_key <- paste(key, toupper(mod), sep = "_")
      all_fits[[fit_key]] <- fit
      fit_results <- rbind(fit_results, data.frame(
        key = key,
        true_model = tree_info$true_model[i],
        fit_model = toupper(mod),
        link = mod_link,
        loglik = fit$loglik,
        AIC = fit$AIC,
        stringsAsFactors = FALSE
      ))
    } else {
      cat(sprintf("    %s: estimation failed\n", mod))
      fit_results <- rbind(fit_results, data.frame(
        key = key, true_model = tree_info$true_model[i],
        fit_model = toupper(mod), link = mod_link,
        loglik = NA, AIC = NA,
        stringsAsFactors = FALSE
      ))
    }
  }
  cat("\n")
}

# ── 3. Model selection ──────────────────────────────────────────

cat("=== Phase 3: Model selection ===\n\n")

selected <- data.frame()
for (key in unique(fit_results$key)) {
  sub <- fit_results[fit_results$key == key, ]
  true_mod <- sub$true_model[1]
  valid <- sub[is.finite(sub$AIC), ]
  if (nrow(valid) > 0) {
    best <- valid$fit_model[which.min(valid$AIC)]
  } else {
    best <- NA
  }
  selected <- rbind(selected, data.frame(
    key = key, true_model = true_mod,
    selected = best,
    stringsAsFactors = FALSE
  ))
}

# Confusion matrix
ok <- !is.na(selected$selected)
model_levels <- c("DD", "PD", "EP")
conf <- table(
  True = factor(selected$true_model[ok], levels = model_levels),
  Selected = factor(selected$selected[ok], levels = model_levels)
)

cat("Confusion matrix (rows=true, cols=selected):\n")
print(conf)

accuracy <- sum(diag(conf)) / sum(conf)
cat(sprintf("\nOverall accuracy: %.0f%%\n", 100 * accuracy))

cat("\nPer-model recall:\n")
recall <- diag(conf) / pmax(rowSums(conf), 1)
print(round(recall, 2))

# ── 4. Parameter recovery ───────────────────────────────────────

cat("\n=== Phase 4: Parameter recovery ===\n\n")

recovery <- data.frame()
for (mod_name in names(true_models)) {
  cfg       <- true_models[[mod_name]]
  true_pars <- cfg$pars
  par_names <- names(
    estimate_rates_control("gam", length(true_pars))
  )[seq_along(true_pars)]
  # Get par names from first successful fit
  mod_key <- toupper(mod_name)
  for (key in tree_info$key[tree_info$true_model == mod_name]) {
    fit_key <- paste(key, mod_key, sep = "_")
    fit <- all_fits[[fit_key]]
    if (!is.null(fit)) {
      par_names <- names(fit$pars)
      break
    }
  }

  for (key in tree_info$key[tree_info$true_model == mod_name]) {
    fit_key <- paste(key, mod_key, sep = "_")
    fit <- all_fits[[fit_key]]
    if (!is.null(fit)) {
      for (j in seq_along(true_pars)) {
        recovery <- rbind(recovery, data.frame(
          model     = mod_name,
          parameter = par_names[j],
          true_val  = true_pars[j],
          estimate  = fit$pars[j],
          stringsAsFactors = FALSE
        ))
      }
    }
  }
}

# Summary table
if (nrow(recovery) > 0) {
  recovery$error <- recovery$estimate - recovery$true_val
  summ <- aggregate(
    cbind(estimate, error) ~ model + parameter + true_val,
    data = recovery, FUN = function(x) {
      c(mean = mean(x), sd = sd(x))
    }
  )
  cat("Parameter recovery summary:\n")
  print(summ, digits = 4)
}

# ── 5. Plots ────────────────────────────────────────────────────

cat("\n=== Phase 5: Generating plots ===\n\n")

# 5a. Parameter recovery boxplots
if (nrow(recovery) > 0) {
  ref <- unique(recovery[, c("model", "parameter", "true_val")])

  p1 <- ggplot(recovery, aes(x = parameter, y = estimate)) +
    geom_boxplot(fill = "steelblue", alpha = 0.4) +
    geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
    geom_point(
      data = ref,
      aes(x = parameter, y = true_val),
      colour = "red", shape = 4,
      size = 3, stroke = 1.5
    ) +
    facet_wrap(~model, scales = "free") +
    labs(
      title = "Parameter recovery",
      subtitle = "Red X = true value",
      y = "Estimate"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(
      angle = 45, hjust = 1
    ))

  ggsave(
    "simulation_study_recovery.pdf",
    p1, width = 10, height = 6
  )
  cat("  Saved simulation_study_recovery.pdf\n")
}

# 5b. Confusion matrix heatmap
conf_pct <- sweep(conf, 1, pmax(rowSums(conf), 1), "/")
conf_df <- as.data.frame(as.table(conf_pct))
names(conf_df) <- c("True", "Selected", "Recall")

p2 <- ggplot(
  conf_df,
  aes(x = Selected, y = True, fill = Recall)
) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * Recall)),
            size = 5) +
  scale_fill_gradient(
    low = "white", high = "steelblue",
    limits = c(0, 1)
  ) +
  labs(
    title = sprintf(
      "Model selection confusion (accuracy=%.0f%%)",
      100 * accuracy
    ),
    x = "Selected model",
    y = "True model"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(
  "simulation_study_confusion.pdf",
  p2, width = 5, height = 4
)
cat("  Saved simulation_study_confusion.pdf\n")

# ── Save ─────────────────────────────────────────────────────────

save(
  tree_info, fit_results, selected, conf,
  recovery, all_fits, true_models,
  file = "simulation_study_results.RData"
)
cat("\nSaved simulation_study_results.RData\n")

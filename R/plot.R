PDTT_plot <- function(tree) {
    ct <- max(tree$brts)
    times <- seq(0, ct, length.out = 1000)
    times <- sort(times, decreasing = TRUE)
    PD <- NULL
    for (i in seq_along(times)) {
        G <- GPD(times[i], tree[-1, ])

        PD <- rbind(PD,
                   data.frame(time = rep(times[i], nrow(G)),
                              P = colSums(G) / (nrow(G) - 1),
                              lineage = as.character(seq_len(nrow(G)))))
    }
    g1 <- ggplot2::ggplot(PD) +
          ggplot2::geom_line(ggplot2::aes_string(x = "time",
                                                 y = "P",
                                                 colour = "lineage",
                                                 alpha = 0.5)) +
          ggplot2::theme(legend.position = "none",
                         panel.grid.major = ggplot2::element_blank(),
                         panel.grid.minor = ggplot2::element_blank(),
                         panel.background = ggplot2::element_blank(),
                         axis.line = ggplot2::element_line(colour = "black"))
    return(g1)
}

learning.curves <- function(...) {
  dfs <- list(...)

  for(df in dfs) {
    if(!all(c('par1', 'par2', 'par3', 'par4') %in% names(df))) {
      stop("Each data frame should have columns: par1, par2, par3, and par4")
    }
  }

  combined_df <- dplyr::bind_rows(lapply(seq_along(dfs), function(i) {
    dfs[[i]] %>%
      dplyr::mutate(iteration = 1:nrow(dfs[[i]]), algorithm = paste("alg", i, sep = "")) %>%
      tidyr::gather(key = "parameter", value = "value", -iteration, -algorithm)
  }))

  p <- ggplot2::ggplot(combined_df, ggplot2::aes(x = iteration, y = value, color = algorithm)) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~parameter, scales = "free", ncol = 2) +
    ggplot2::labs(x = "Iteration number") +
    ggplot2::scale_color_discrete(name = "Algorithm") +
    ggplot2::theme_minimal()

  print(p)
}

lik_vis <- function(...) {
  vec_list <- list(...)

  df <- data.frame(iteration = seq_len(max(sapply(vec_list, length))))

  for (i in seq_along(vec_list)) {
    vec_name <- names(vec_list)[i]
    if (is.null(vec_name)) vec_name <- paste0("vec", i)

    df[[vec_name]] <- vec_list[[i]]
  }

  df_long <- tidyr::gather(df, key = "vector", value = "value", -iteration)

  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = iteration, y = value, color = vector)) +
    ggplot2::geom_line() +
    ggplot2::labs(title = "Plot for multiple vectors",
         x = "Iteration number",
         y = "Value")+
    ggplot2::theme_minimal()

  print(p)
}

process_and_plot <- function(file_path) {
  load(file_path)

  computation_time <- time1 - time0
  cat("Data loaded\n")
  computation_hours <- computation_time['elapsed'] %/% 3600
  computation_minutes <- (computation_time['elapsed'] %% 3600) %/% 60
  computation_seconds <- computation_time['elapsed'] %% 60

  print(paste("Simulation took ",computation_hours, "hours,", computation_minutes, "minutes and", round(computation_seconds), "seconds"))

  parameters <- c(betaN_interval = betaN_interval, betaP_interval =  betaP_interval,
                lambda_interval = lambda_interval, max_lin = max_lin, max_tries = max_tries,
                mu_interval = mu_interval, n_trees = n_trees, phylo_name = phylo_name)

  for (param_name in names(parameters)) {
    print(paste(param_name, ":", parameters[param_name]))
  }

  cat("Extracting trees information...\n")

  parameters <- results$param

  num_simulations <- length(results$trees)
  num_completed <- sum(!sapply(results$loglik_estimation, inherits, "try-error"))

  completed_indices <- which(!sapply(results$loglik_estimation, inherits, "try-error"))

  trees = results$trees[completed_indices]

  num_extinct_per_tree <- sapply(trees, count_extinct_species)

  cat("Number of simulations attempted:", num_simulations, "\n")
  cat("Number of simulations completed successfully:", num_completed, "\n")

  params_completed <- data.frame(mu = results$param$mu,
  lambda = results$param$lambda,
  betaN = results$param$betaN,
  betaP = results$param$betaP,
  loglik = unlist(results$loglik_estimation[completed_indices]),
  nspecies =  num_extinct_per_tree
  )

  cat("Creating plots...\n")

  plot1 = ggplot2::ggplot(params_completed, ggplot2::aes(x = mu, y = lambda, color = loglik)) +
    ggplot2::geom_point(alpha = 0.2) +
    ggplot2::scale_color_gradient(low = "blue", high = "red") +
    ggplot2::labs(title = "log likelihood estimation", x = "mu", y = "lambda") +
    ggplot2::theme_minimal()

  plot2 = ggplot2::ggplot(params_completed, ggplot2::aes(x = betaN, y = betaP, color = loglik)) +
    ggplot2::geom_point(alpha = 0.2) +
    ggplot2::scale_color_gradient(low = "blue", high = "red") +
    ggplot2::labs(title = "log likelihood estimation", x = "betaN", y = "betaP") +
    ggplot2::theme_minimal()

  plot3 = ggplot2::ggplot(params_completed, ggplot2::aes(x = mu, y = lambda, color = nspecies)) +
    ggplot2::geom_point(alpha = 0.2) +
    ggplot2::scale_color_gradient(low = "blue", high = "red") +
    ggplot2::labs(title = "Number of extinct species", x = "mu", y = "lambda") +
    ggplot2::theme_minimal()

  plot4 = ggplot2::ggplot(params_completed, ggplot2::aes(x = betaN, y = betaP, color = nspecies)) +
    ggplot2::geom_point(alpha = 0.2) +
    ggplot2::scale_color_gradient(low = "blue", high = "red") +
    ggplot2::labs(title = "Number of extinct species", x = "betaN", y = "betaP") +
    ggplot2::theme_minimal()

  combined_plots <- (plot1 | plot2) / (plot3 | plot4)

  print(combined_plots)

  plot5 = ggplot2::ggplot(params_completed, ggplot2::aes(x = nspecies, y =  loglik)) +
    ggplot2::geom_point(alpha = 0.2) +
    ggplot2::theme_minimal()

  print(plot5)
}

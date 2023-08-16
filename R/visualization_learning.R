learning.curves <- function(...) {
  dfs <- list(...)
  
  # Ensure all data.frames have the expected columns
  for(df in dfs) {
    if(!all(c('par1', 'par2', 'par3', 'par4') %in% names(df))) {
      stop("Each data frame should have columns: par1, par2, par3, and par4")
    }
  }
  
  # Convert data.frames into a long format and bind them together
  combined_df <- bind_rows(lapply(seq_along(dfs), function(i) {
    dfs[[i]] %>%
      mutate(iteration = 1:nrow(dfs[[i]]), algorithm = paste("alg", i, sep = "")) %>%
      gather(key = "parameter", value = "value", -iteration, -algorithm)
  }))
  
  # Plot
  p <- ggplot(combined_df, aes(x = iteration, y = value, color = algorithm)) +
    geom_line() +
    facet_wrap(~parameter, scales = "free", ncol = 2) +
    labs(x = "Iteration number") +
    scale_color_discrete(name = "Algorithm") +
    theme_minimal()
  
  print(p)
}

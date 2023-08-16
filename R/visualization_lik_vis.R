lik_vis <- function(...) {
  vec_list <- list(...)
  
  # Convert list of vectors to a data.frame
  df <- data.frame(iteration = seq_len(max(sapply(vec_list, length))))
  
  for (i in seq_along(vec_list)) {
    vec_name <- names(vec_list)[i]
    if (is.null(vec_name)) vec_name <- paste0("vec", i)
    
    df[[vec_name]] <- vec_list[[i]]
  }
  
  # Convert to long format for plotting
  df_long <- df %>% gather(key = "vector", value = "value", -iteration)
  
  p <- ggplot(df_long, aes(x = iteration, y = value, color = vector)) + 
    geom_line() + 
    labs(title = "Plot for multiple vectors",
         x = "Iteration number",
         y = "Value")+
    theme_minimal()
  
  print(p)
}

# Cr
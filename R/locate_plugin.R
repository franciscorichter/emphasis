#' function to locate a diversification model
#' @param PlugIn the name of the plugin
#' @export
#' @return path of the plugin.
locate_plugin <- function(PlugIn) {
  path <- getLoadedDLLs()[[PlugIn]][["path"]]
  if (!is.character(path)) {
    msg <- paste0(c("Can't locate plugin. Did you load library '",
                   PlugIn, "'?"))
    warning(msg)
    return(PlugIn)  # could be absolute path
  }
  return(path)
}

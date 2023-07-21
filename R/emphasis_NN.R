.onLoad <- function(libname, pkgname) {
  model_path <- system.file("data", "NNtest.pt", package = "emphasis")
  assign("cnn_ltt", torch::torch_load(model_path), envir = .GlobalEnv)
}
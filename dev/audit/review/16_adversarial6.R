.libPaths(c(Sys.getenv("EMPH_LIB"), .libPaths()))
suppressMessages(library(emphasis))
brts1 <- c(5,3.822147,3.059192,1.966150,1.604249,1.358484,0.622046,0.467531,
           0.393533,0.285495,0.276148,0.225482,0.201295,0.154060,0.057677,
           0.044974,0.041380,0.038051,0.013895)
cat("build:", Sys.getenv("EMPH_LIB"), "\n")
# maxN is meaningless for the BDI sampler, but the new guard is applied to all
r <- tryCatch(estimate_rates(brts1, model="cr", method="mcem", init_pars=c(0.6,0.1),
   control=list(sampling="bdi", num_trees=300L, maxN=200L, max_iter=1L,
                lower_bound=c(0,0), upper_bound=c(2,1), num_threads=1L)),
   error=function(e) paste("ERROR:", conditionMessage(e)))
cat("bdi + maxN<num_trees ->", if (is.character(r)) r else
    sprintf("ok, pars=%s", paste(round(r$pars,4), collapse=",")), "\n")

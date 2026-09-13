lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, .libPaths()))
suppressMessages(library(emphasis))
set.seed(5)
tp <- 5
lam <- c(runif(40, 1e-4, 5), 0.5,0.5,0.5,1e-8)
mu  <- c(runif(40, 1e-4, 5), 0.5-1e-13, 0.5-1e-14, 0.5-1e-11, 5e-9)
out <- numeric(0)
for (i in seq_along(lam)) {
  if (mu[i] >= lam[i]) next
  for (n in c(0L,1L,3L)) for (sg in list(c(0,2),c(1,4.9),c(2,3)))
    out <- c(out, emphasis:::.bdi_integral_cr(sg[1],sg[2],n,2L,lam[i],mu[i],tp),
                  emphasis:::.bdi_p_cr(sg[1],lam[i],mu[i],tp))
}
cat(sprintf("n=%d  md5=%s\n", length(out), digest::digest(out, algo="md5")))
saveRDS(out, file.path("/tmp", paste0("bit_", basename(lib), ".rds")))

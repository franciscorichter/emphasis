# bdi_invariance.R -- a bit-identity corpus for the conditional (BDI) sampler.
#
# Parent recording adds bookkeeping to .bdi_augment_one but must not move the
# random stream: every draw has to come out byte for byte what it was before.
# This captures the corpus, and, run again, compares against it.
#
#   Rscript dev/bdi_invariance.R capture   # write dev/bdi_invariance.rds
#   Rscript dev/bdi_invariance.R check     # compare against it
suppressMessages(devtools::load_all(".", quiet = TRUE))
mode <- commandArgs(trailingOnly = TRUE)[1]; if (is.na(mode)) mode <- "check"
STORE <- "dev/bdi_invariance.rds"

CELLS <- list(
  list(tag = "cr-25-rho1",   n = 25, model = "cr", rho = 1.0,  pars = c(0.5, 0.15)),
  list(tag = "cr-25-rho05",  n = 25, model = "cr", rho = 0.5,  pars = c(0.5, 0.15)),
  list(tag = "cr-60-rho1",   n = 60, model = "cr", rho = 1.0,  pars = c(0.5, 0.15)),
  list(tag = "dd-25-rho1",   n = 25, model = "dd", rho = 1.0,  pars = c(0.6, -0.01, 0.15, 0)),
  list(tag = "dd-60-rho1",   n = 60, model = "dd", rho = 1.0,  pars = c(0.6, -0.01, 0.15, 0)),
  list(tag = "dd-60-rho05",  n = 60, model = "dd", rho = 0.5,  pars = c(0.6, -0.01, 0.15, 0))
)

one <- function(cell) {
  set.seed(1000 + nchar(cell$tag))
  phy  <- ape::rcoal(cell$n)
  brts <- emphasis:::.extract_brts(phy)
  mb   <- emphasis:::.resolve_model(cell$model)
  p8   <- emphasis:::.expand_pars(cell$pars, mb)
  set.seed(7)
  a <- emphasis:::.augment_tree_bdi(brts, p8, model_bin = mb[1:3],
                                    sample_size = 25L, link = 0L, rho = cell$rho)
  list(tag = cell$tag,
       logg = a$logg, logf = a$logf, weights = a$weights, fhat = a$fhat,
       n_valid = a$n_valid, n_rejected = a$n_rejected,
       nrow = vapply(a$trees, nrow, 1L),
       brts = lapply(a$trees, function(t) t$brts),
       t_ext = lapply(a$trees, function(t) t$t_ext),
       n = lapply(a$trees, function(t) t$n))
}
now <- lapply(CELLS, one); names(now) <- vapply(CELLS, `[[`, "", "tag")

if (mode == "capture") {
  saveRDS(now, STORE)
  cat(sprintf("[invariance] captured %d cells, %d draws\n", length(now),
              sum(vapply(now, function(x) x$n_valid, 1L))))
  quit(save = "no")
}
if (!file.exists(STORE)) stop("no corpus at ", STORE, " -- run `capture` first")
was <- readRDS(STORE)
bad <- 0L
for (tag in names(was)) {
  d <- all.equal(was[[tag]], now[[tag]], tolerance = 0)
  if (isTRUE(d)) cat(sprintf("  ok   %s (%d draws)\n", tag, now[[tag]]$n_valid))
  else { bad <- bad + 1L; cat(sprintf("  MOVED %s\n", tag)); cat(paste0("        ", d, "\n")) }
}
cat(sprintf("[invariance] %d of %d cells bit-identical\n", length(was) - bad, length(was)))
if (bad) quit(save = "no", status = 1)

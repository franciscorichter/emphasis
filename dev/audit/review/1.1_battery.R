lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, .libPaths())); library(emphasis)
set.seed(7)

mk_tree <- function(k, tp = k + 1) {
  data.frame(brts = c(seq_len(k), tp), n = c(2:(k + 1), k + 2), t_ext = 1e11,
             pd = 0, tip_start = 0, id = c(seq_len(k) - 1L, -1L), parent_id = -1L)
}
# augmented tree with extinctions, nonzero pd / tip_start (exercises M and D slots)
aug_tree <- function(k, seed) {
  set.seed(seed)
  brts <- sort(runif(k, 0.05, 4)); tp <- 5
  n <- 2L + cumsum(sample(c(1L, -1L), k, TRUE, prob = c(.75, .25)))
  n <- pmax(n, 2L)
  text <- ifelse(runif(k) < 0.25, brts + runif(k, 0, 1), 1e11)
  data.frame(brts = c(brts, tp), n = c(n, max(n)),
             t_ext = c(text, 1e11), pd = c(runif(k, 0, 3), 0),
             tip_start = c(runif(k, 0, 1), 0),
             id = c(seq_len(k) - 1L, -1L),
             parent_id = c(rep(-1L, k), -1L))
}
ev <- function(p, tr, model, link, rho = 1)
  emphasis:::eval_logf(p, if (is.data.frame(tr)) list(tr) else tr,
                       model = as.integer(model), link = as.integer(link), rho = rho)

models <- list(cr = c(0L,0L,0L), dd = c(1L,0L,0L), m = c(1L,1L,0L), d = c(1L,1L,1L))
out <- list(); i <- 0
for (mn in names(models)) for (link in 0:2) for (tn in c("t3","t4","t8","a6","a12")) {
  tr <- switch(tn, t3 = mk_tree(3), t4 = mk_tree(4), t8 = mk_tree(8),
               a6 = aug_tree(6, 11), a12 = aug_tree(12, 23))
  pset <- list(
    base    = c(0.5, -0.02, 0.01, 0.01, 0.2, 0.005, 0.0, 0.0),
    zerob0  = c(0.0, 0.0, 0, 0, 0.2, 0, 0, 0),
    negb0   = c(-0.4, 0.0, 0, 0, 0.2, 0, 0, 0),
    mueqlam = c(0.3, 0, 0, 0, 0.3, 0, 0, 0),
    mugtlam = c(0.3, 0, 0, 0, 0.9, 0, 0, 0),
    tiny    = c(1e-12, 0, 0, 0, 1e-12, 0, 0, 0),
    huge    = c(500, 0, 0, 0, 1, 0, 0, 0),
    ddzero  = c(0.4, -0.1, 0, 0, 0.1, 0, 0, 0),
    expneg  = c(-800, 0, 0, 0, -800, 0, 0, 0),
    nanpar  = c(NaN, 0, 0, 0, 0.2, 0, 0, 0)
  )
  for (pn in names(pset)) {
    r <- try(ev(pset[[pn]], tr, models[[mn]], link), silent = TRUE)
    i <- i + 1
    out[[i]] <- data.frame(model = mn, link = link, tree = tn, pars = pn,
      logf = if (inherits(r, "try-error")) NA_real_ else r$logf[1],
      logg = if (inherits(r, "try-error")) NA_real_ else r$logg[1],
      err  = inherits(r, "try-error"))
  }
}
res <- do.call(rbind, out)
saveRDS(res, commandArgs(TRUE)[2])
cat("rows", nrow(res), "\n")

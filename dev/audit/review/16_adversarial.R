lib <- Sys.getenv("EMPH_LIB")
.libPaths(c(lib, .libPaths()))
suppressMessages(library(emphasis))
ns <- asNamespace("emphasis")
ex <- function(p, m) ns$.expand_pars(p, m)

brts1 <- c(5.000000, 3.822147, 3.059192, 1.966150, 1.604249, 1.358484,
           0.622046, 0.467531, 0.393533, 0.285495, 0.276148, 0.225482,
           0.201295, 0.154060, 0.057677, 0.044974, 0.041380, 0.038051,
           0.013895)
cr <- c(0L,0L,0L); dd <- c(1L,0L,0L); dmod <- c(0L,0L,1L)

run <- function(tag, pars, model, link, lb, ub, ss=20L, maxN=2000L, it=4L,
                max_missing=1e4, tol=1e-2, patience=3L, seed=7) {
  set.seed(seed)
  r <- tryCatch(ns$.mcem_dynamic_fresh(brts1, ex(pars, model), sample_size=ss,
        maxN=maxN, max_missing=max_missing, lower_bound=ex(lb,model),
        upper_bound=ex(ub,model), max_iter=it, xtol=1e-3, tol=tol,
        patience=patience, num_threads=1L, verbose=FALSE, model=model,
        link=link, max_time=60), error=function(e) paste("ERROR:", conditionMessage(e)))
  if (is.character(r)) { cat(sprintf("%-28s %s\n", tag, r)); return(invisible(NULL)) }
  nr <- if (is.null(r$mcem)) 0L else nrow(r$mcem)
  cat(sprintf("%-28s stop=%-15s iters=%s nrow=%s fin=%s n_fail=%s maxN=%s pars=%s\n",
      tag, r$stop_reason, r$iterations, nr, r$final_estep, r$n_failed, r$maxN,
      paste(sprintf("%.4g", ns$.contract_pars(r$pars, model)), collapse=",")))
  invisible(r)
}

cat("---- A. links x models ----\n")
run("cr/linear",      c(0.5,0.1), cr,   0L, c(0,0), c(2,1))
run("cr/exponential", log(c(0.5,0.1)), cr, 1L, c(-5,-5), c(2,2))
run("cr/gaussian",    c(0.5,0.1), cr,   2L, c(0,0), c(2,1))
run("dd/linear",      c(0.5,-0.01,0.1), dd, 0L, c(0,-1,0), c(2,0,1))
run("dd/exponential", c(-0.7,-0.01,-2.3), dd, 1L, c(-5,-1,-5), c(2,0,2))
run("d/linear",       c(0.5,0.01,0.1), dmod, 0L, c(0,-1,0), c(2,1,1))

cat("---- B. boundary rates ----\n")
run("mu == lambda",  c(0.5,0.5), cr, 0L, c(0,0), c(2,2))
run("mu > lambda",   c(0.3,0.9), cr, 0L, c(0,0), c(2,2))
run("mu ~ 0",        c(0.5,1e-8), cr, 0L, c(0,0), c(2,2))
run("lambda ~ 0",    c(1e-6,0.1), cr, 0L, c(0,0), c(2,2))
run("sample_size=1", c(0.5,0.1), cr, 0L, c(0,0), c(2,1), ss=1L)
run("maxN==ss",      c(0.5,0.1), cr, 0L, c(0,0), c(2,1), ss=20L, maxN=20L)
run("maxN<ss",       c(0.5,0.1), cr, 0L, c(0,0), c(2,1), ss=20L, maxN=5L, it=3L)
run("max_iter=1",    c(0.5,0.1), cr, 0L, c(0,0), c(2,1), it=1L)

cat("---- C. exponential link near log-rate 0 (metric floor) ----\n")
r <- run("exp near logl=0", c(0.02,-1.6), cr, 1L, c(-5,-5), c(2,2), it=12L, seed=11)
if (!is.null(r)) {
  m <- r$mcem
  print(round(m[, c("par1","par5","delta_max","abs_step")], 6))
}

cat("---- D. tiny / large trees ----\n")
set.seed(5); b_small <- sort(c(3, 1.2), decreasing=TRUE)
set.seed(5)
r <- tryCatch(ns$.mcem_dynamic_fresh(b_small, ex(c(0.5,0.1),cr), sample_size=10L,
      maxN=2000L, max_missing=1e4, lower_bound=ex(c(0,0),cr), upper_bound=ex(c(2,1),cr),
      max_iter=3L, xtol=1e-3, tol=1e-2, patience=3L, num_threads=1L, verbose=FALSE,
      model=cr, link=0L, max_time=60), error=function(e) conditionMessage(e))
cat("3-tip tree:", if (is.character(r)) r else sprintf("stop=%s iters=%d nrow=%d", r$stop_reason, r$iterations, nrow(r$mcem)), "\n")

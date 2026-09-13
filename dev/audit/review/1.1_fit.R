args <- commandArgs(TRUE); lib <- args[1]; tag <- basename(lib)
.libPaths(c(lib, .libPaths())); suppressMessages(library(emphasis)); suppressPackageStartupMessages(library(DDD))
set.seed(11)
sim <- DDD::dd_sim(pars = c(1.5, 0.4, 1.1 / 0.12), age = 6, ddmodel = 1)
brts <- sort(as.numeric(sim$brts), decreasing = TRUE)
pars <- c(1.5, -0.12, 0.4, 0); lb <- c(0.01, -1, 0.001, 0); ub <- c(5, 0, 2, 0)
rep_fit <- function(fit, lab) {
  fi <- fit$details$final_IS
  cat(sprintf("[%s] %s: stop=%s iters=%s pars=%s loglik=%s var=%s | final_IS: n=%d +Inf=%d -Inf=%d NaN=%d ESS=%s\n", tag, lab,
      format(fit$details$stop_reason), format(fit$details$iterations), paste(round(fit$pars, 4), collapse = " "),
      format(fit$loglik), format(fit$loglik_var),
      if (is.null(fi)) 0L else length(fi$logf), if (is.null(fi)) 0L else sum(fi$logf == Inf), if (is.null(fi)) 0L else sum(fi$logf == -Inf),
      if (is.null(fi)) 0L else sum(is.nan(fi$logf)), if (is.null(fi)) NA else format(fi$ESS)))
}
t0 <- proc.time()[3]
set.seed(1)
f1 <- try(estimate_rates(brts, method = "mcem", model = "dd", init_pars = pars,
        control = list(lower_bound = lb, upper_bound = ub, sample_size = 200L, max_iter = 4L,
                       max_missing = 30L, num_threads = 1L, sampling = "bdi", verbose = FALSE)))
if (!inherits(f1, "try-error")) rep_fit(f1, "dd/lin BDI from gen pars") else cat(tag, "BDI fit error:", f1, "\n")
set.seed(1)
f2 <- try(estimate_rates(brts, method = "mcem", model = "dd", init_pars = pars,
        control = list(lower_bound = lb, upper_bound = ub, sample_size = 200L, max_iter = 4L,
                       max_missing = 30L, num_threads = 1L, sampling = "dynamic_fresh", verbose = FALSE)))
if (!inherits(f2, "try-error")) rep_fit(f2, "dd/lin thinning from gen pars") else cat(tag, "thinning fit error:", f2, "\n")
# direct E-step draw: count of +Inf logf under each build
set.seed(3)
e <- emphasis:::.augment_tree_bdi(brts, pars, model_bin = c(1L,0L,0L), sample_size = 200L, max_missing = 30L, link = 0L, rho = 1)
cat(sprintf("[%s] .augment_tree_bdi: n_trees=%d +Inf=%d -Inf=%d NaN=%d fhat=%s names=%s\n", tag, length(e$trees),
    sum(e$logf == Inf), sum(e$logf == -Inf), sum(is.nan(e$logf)), format(e$fhat), paste(names(e), collapse = ",")))
cat(sprintf("[%s] elapsed %.0fs\n", tag, proc.time()[3] - t0))

# Refit a sample of the validation study's failing cells under the wave-1 rule
# (stop = "step") and the ascent rule (stop = "ascent"), and report, per cell:
# fraction converged, log-likelihood deficit against the exact MLE, iterations,
# and total E-step draws.
suppressWarnings(suppressMessages({
  devtools::load_all(".", quiet = TRUE)
  library(ape); library(DDD); library(parallel)
}))
VAL <- "/Users/pancho/Code/emphasis/dev/validation"
src <- new.env()
sys.source(file.path(VAL, "R", "00-common.R"), envir = src)
for (nm in ls(src)) assign(nm, get(nm, envir = src))

TR  <- readRDS(file.path(VAL, "data", "trees-main.rds"))$trees
REF <- readRDS(file.path(VAL, "data", "reference-main.rds"))$refs

CELLS <- rbind(
  data.frame(cell = "cr-n020-e09-l100", config = "C1", sampler = "bdi",           init = "far", model = "cr"),
  data.frame(cell = "cr-n050-e03-l100", config = "C1", sampler = "bdi",           init = "far", model = "cr"),
  data.frame(cell = "cr-n020-e09-l100", config = "C2", sampler = "dynamic_fresh", init = "far", model = "cr"),
  data.frame(cell = "cr-n050-e06-l100", config = "C2", sampler = "dynamic_fresh", init = "far", model = "cr"),
  data.frame(cell = "cr-n100-e06-l100", config = "C2", sampler = "dynamic_fresh", init = "far", model = "cr"),
  data.frame(cell = "cr-n050-e09-l100", config = "C4", sampler = "dynamic_fresh", init = "mle", model = "cr"),
  data.frame(cell = "cr-n100-e09-l100", config = "C4", sampler = "dynamic_fresh", init = "mle", model = "cr"),
  data.frame(cell = "cr-n050-e06-l100", config = "C4", sampler = "dynamic_fresh", init = "mle", model = "cr"),
  data.frame(cell = "dd-A",             config = "D2", sampler = "dynamic_fresh", init = "far", model = "dd"),
  data.frame(cell = "dd-A",             config = "D4", sampler = "dynamic_fresh", init = "mle", model = "dd"),
  stringsAsFactors = FALSE)

N_TREES <- as.integer(Sys.getenv("N_TREES", "4"))
# rule spec: "step" or "ascent@<ascent_tol>"
RULES   <- strsplit(Sys.getenv("RULES", "step,ascent@0.01"), ",")[[1]]
CELLSEL <- Sys.getenv("CELLSEL", "")
MAXIT   <- 400L
MAXTIME <- as.numeric(Sys.getenv("MAXTIME", "600"))

jobs <- list()
for (r in seq_len(nrow(CELLS))) {
  ids <- grep(paste0("^", CELLS$cell[r], "-t"), names(TR), value = TRUE)
  ids <- head(ids, N_TREES)
  for (id in ids) for (rule in RULES)
    jobs[[length(jobs) + 1L]] <- c(as.list(CELLS[r, ]), list(tree_id = id, rule = rule))
}
if (nzchar(CELLSEL)) {
  keep <- strsplit(CELLSEL, ",")[[1]]
  jobs <- Filter(function(j) paste(j$config, j$cell) %in% keep, jobs)
}
cat("jobs:", length(jobs), "\n")

one <- function(job) {
  tt <- TR[[job$tree_id]]; rf <- REF[[job$tree_id]]
  brts <- tt$brts
  model <- job$model
  box <- if (model == "cr") val_box_cr(tt$lambda_gen, 1) else
    val_box_dd(tt$lambda0_gen, tt$mu0_gen, tt$K_gen)
  init <- if (model == "cr") val_init_cr(job$init, tt$lambda_gen, rf$mle) else
    val_init_dd(job$init, tt$lambda0_gen, tt$mu0_gen, tt$K_gen, rf$mle)
  init <- pmin(pmax(init, box$lower), box$upper)
  ctrl <- list(lower_bound = box$lower, upper_bound = box$upper,
               sampling = job$sampler, sample_size = 200L, maxN = 4000L,
               max_iter = MAXIT, max_missing = 1e4, num_threads = 1L,
               max_time = MAXTIME, verbose = FALSE)
  spec <- strsplit(job$rule, "@", fixed = TRUE)[[1]]
  if (identical(spec[1], "ascent_nogrow")) {
    ctrl$stop <- "ascent"; ctrl$grow_sample <- FALSE
  } else {
    ctrl$stop <- spec[1]
  }
  if (length(spec) > 1L) ctrl$ascent_tol <- as.numeric(spec[2])
  set.seed(val_seed(paste(job$tree_id, job$config, job$rule)))
  t0 <- proc.time()[3]
  fit <- tryCatch(
    estimate_rates(brts, method = "mcem", model = model, init_pars = init,
                   control = ctrl, link = "linear", cond = NULL),
    error = function(e) NULL)
  el <- proc.time()[3] - t0
  if (is.null(fit)) return(data.frame(job, ok = FALSE, stringsAsFactors = FALSE))
  pars <- as.numeric(fit$pars)
  lx <- if (model == "dd") rf$lx else NULL
  ll_hat <- tryCatch(ll_exact(pars, brts, model = model, lx = lx),
                     error = function(e) NA_real_)
  tr <- fit$details$mcem
  draws <- if (!is.null(fit$total_draws) && is.finite(fit$total_draws))
    fit$total_draws else sum(tr$num_trees, na.rm = TRUE)
  data.frame(job, ok = TRUE,
             stop_reason = fit$stop_reason,
             converged = fit$stop_reason %in% c("converged", "converged_step"),
             iterations = fit$iterations,
             draws = draws,
             S_final = if (is.null(fit$sample_size_final)) NA_integer_ else
               fit$sample_size_final,
             delta_ell = ll_hat - rf$mle_loglik,
             secs = el,
             stringsAsFactors = FALSE)
}

nc <- as.integer(Sys.getenv("NCORES", "6"))
res <- mclapply(jobs, function(j) tryCatch(one(j), error = function(e)
  data.frame(as.data.frame(j), ok = FALSE, stringsAsFactors = FALSE)),
  mc.cores = nc, mc.preschedule = FALSE)
out <- do.call(rbind, lapply(res, function(x)
  if (is.data.frame(x)) x else NULL))
saveRDS(out, Sys.getenv("OUT", "dev/stopping-rule/measure-out.rds"))
cat("done:", nrow(out), "rows\n")

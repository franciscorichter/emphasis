.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths())); library(emphasis)
inc <- system.file("include", package = "RcppParallel")
cat("RcppParallel", as.character(packageVersion("RcppParallel")), "\n")
cat(grep("TBB_VERSION_MAJOR|TBB_VERSION_MINOR", readLines(file.path(inc, "tbb/tbb_stddef.h")), value = TRUE), sep = "\n")
cat("oneapi dir present:", dir.exists(file.path(inc, "oneapi")), "\n")
brts <- c(6, 4.5, 3.2, 2.5, 1.8, 1.1, 0.7, 0.3)
model <- c(0L,0L,0L)
ip <- emphasis:::.expand_pars(c(0.5, 0.1), model)
lb <- rep(0,8); ub <- c(2,0,0,0,1,0,0,0)
ncall <- 0L
cond_fail3 <- function(p) { ncall <<- ncall + 1L; if (ncall == 3L) stop("boom from R conditional"); 0 }
for (k in 1:5) {
  ncall <- 0L
  r <- tryCatch(emphasis:::em_cpp(brts, ip, 20L, 1000L, 1e4, 1e6, lb, ub, 1e-3, 2L, FALSE, model, 0L, 1.0, cond_fail3), error=function(e) paste("R error:", conditionMessage(e)))
  cat(sprintf("run %d: calls=%d -> %s\n", k, ncall, if (is.list(r)) "returned list" else r))
}
cond_warn <- function(p) { warning("warn in conditional"); 0 }
r <- tryCatch(emphasis:::em_cpp(brts, ip, 20L, 1000L, 1e4, 1e6, lb, ub, 1e-3, 2L, FALSE, model, 0L, 1.0, cond_warn), error=function(e) paste("R error:", conditionMessage(e)))
cat("warning-in-conditional:", if (is.list(r)) "ok" else r, "\n")
# Still healthy afterwards?
r <- emphasis:::em_cpp(brts, ip, 20L, 1000L, 1e4, 1e6, lb, ub, 1e-3, 2L, FALSE, model, 0L, 1.0, NULL)
cat("post-error plain fit ok, nlopt =", r$nlopt, "\n")

.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib-wave1", .libPaths())); library(emphasis)
brts_dd <- c(6, 4.848493, 4.401821, 3.108164, 3.073914, 2.835023, 1.838828, 0.50463)
p8 <- c(1.5,-0.12,0,0,0.4,0,0,0)
set.seed(11)
e <- emphasis:::.augment_tree_bdi(brts_dd, p8, model_bin=c(1L,0L,0L), sample_size=200L,
                                  max_missing=30L, link=0L, rho=1)
cat("n_valid",e$n_valid," n_nonfinite",e$n_nonfinite," n_rejected",e$n_rejected,
    " n_rej_mm",e$n_rejected_max_missing," attempts",e$n_attempts,
    " length(trees)",length(e$trees), " all finite logf:", all(is.finite(e$logf)),
    " sum(!is.finite(logf)) =", sum(!is.finite(e$logf)), "\n")

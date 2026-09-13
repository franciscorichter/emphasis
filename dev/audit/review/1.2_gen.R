.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
library(emphasis)
brts_cr <- c(4.087931,3.824527,3.380355,1.879518,1.837783,1.474373,1.182076,
             1.172046,1.067297,0.847818,0.516806,0.493594,0.415446,0.379772,
             0.37824,0.353089,0.318518,0.223403,0.222619,0.181101,0.172338,
             0.161752,0.076436,0.058077)
set.seed(1)
e <- emphasis:::.augment_tree_bdi(brts_cr, c(1,0,0,0,0.5,0,0,0),
        model_bin=c(0L,0L,0L), sample_size=40L, max_missing=1e4, link=0L, rho=1)
lw <- e$weights; w <- exp(lw-max(lw)); w <- w/sum(w)*length(w)
saveRDS(list(trees=e$trees, w=w, brts=brts_cr), "/tmp/1.2_trees.rds")
cat("n trees", length(e$trees), " sum_w", sum(w), "\n")

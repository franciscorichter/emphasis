.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressMessages(library(emphasis))
brts <- c(5.31, 4.2, 3.8, 3.1, 2.5, 2.2, 1.9, 1.4, 1.1, 0.8, 0.6, 0.4, 0.3, 0.2)
pm0 <- c(0.5,0,0,0, 0,0,0,0)
raw0 <- tryCatch(emphasis:::augment_trees(as.numeric(brts), pm0, sample_size = 50L, maxN = 5000L,
                     max_missing = 200L, max_lambda = 1e3, num_threads = 1L,
                     model = c(0L,0L,0L), link = 0L, rho = 1), error = function(e) e)
print(names(raw0))
print(raw0[setdiff(names(raw0), c("trees","logf","logg"))])
ne0 <- vapply(raw0$trees, function(tr) sum(tr$t_ext == 0), integer(1))
print(table(ne0))
lw <- raw0$logf - raw0$logg
print(summary(lw))
print(unlist(raw0)[1:3])

## H61 follow-up: where do the negative edges on the 2-tip BDI tas come from?
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({library(emphasis); library(ape)})
tr2 <- ape::read.tree(text = "(A:5,B:5);")
L <- DDD::phylo2L(tr2); print(L)
brts <- 5
a <- emphasis:::.augment_tree_bdi(brts, c(0.5,0,0,0,0.3,0,0,0), c(0L,0L,0L), sample_size = 1L)
df <- a$trees[[1]]
La <- emphasis:::.aug_to_Ltable(df, 5, brts, L); print(La)
tas <- DDD::L2phylo(La, dropextinct = FALSE); cat("min edge:", min(tas$edge.length), "\n")

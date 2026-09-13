.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths())); library(emphasis)
# float resolution of simulator times: duplicates / zero-length edges at max_t = 100
res <- data.frame(); L <- NULL
for (r in 1:40) {
  s <- tryCatch(simulate_tree(pars = c(0.1, 0.02), max_t = 100, model = "cr", max_lin = 20000, max_tries = 5), error = function(e) NULL)
  if (is.null(s) || is.null(s$tas)) next
  tas <- s$tas; L <- s$L
  bt <- ape::branching.times(tas)
  ages <- L[, 1]
  res <- rbind(res, data.frame(ntip = ape::Ntip(tas), nL = nrow(L),
     dup_bt = sum(duplicated(bt)),
     min_gap_bt = min(diff(sort(unique(bt)))),
     zero_edges = sum(tas$edge.length <= 0),
     min_edge = min(tas$edge.length),
     dup_birth_L = sum(duplicated(ages)),
     dup_birth_within_1e5 = sum(diff(sort(ages)) < 1e-5) ))
}
print(res); print(colSums(res[, -1]))
cat("float32 granularity at t~100:", 100 * 2^-23, "\n")
if (!is.null(L)) {
  x <- 100 - L[,1]   # forward birth dates
  # float32 representable? (round to 23-bit mantissa)
  f32 <- function(v) { e <- floor(log2(pmax(v, 1e-300))); round(v / 2^(e-23)) * 2^(e-23) }
  cat("forward dates exactly float32-representable (fraction):", mean(abs(x - f32(x)) < 1e-12 | x == 0), "\n")
  cat("distinct birth dates / rows:", length(unique(L[,1])), "/", nrow(L), "\n")
}
# nearest-value parent matching (simulate.R:428) under augmentation of a big tree: ambiguous when duplicates exist

.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths())); library(emphasis)
ex <- emphasis:::.expand_pars; ct <- emphasis:::.contract_pars
# identity case used by the test
mb <- c(1L,1L,1L); compact <- c(0.5,0.02,0.03,0.01,0.1,-0.01,0.02,0.005)
cat("mb=111 expand is identity:", identical(ex(compact, mb), compact), "\n")
for (mb in list(c(1L,0L,0L), c(0L,0L,1L), c(1L,0L,1L), c(0L,0L,0L), c(0L,1L,0L))) {
  n <- 2L + 2L*sum(mb); compact <- seq_len(n) * 0.1
  f8 <- ex(compact, mb); rt <- ct(f8, mb)
  cat(sprintf("mb=%s  full8=%s  roundtrip_ok=%s\n", paste(mb,collapse=""),
      paste(f8,collapse=","), isTRUE(all.equal(rt, compact))))
}

x <- readLines("/Users/pancho/Code/emphasis/tests/testthat/test-mcem-thinning.R")
st <- grep("^test_that", x)
hdr <- x[1:(st[1]-1)]
ends <- c(st[-1]-1, length(x))
d <- "/Users/pancho/Code/emphasis/dev/audit/review"
for (i in seq_along(st)) writeLines(c(hdr, x[st[i]:ends[i]]), file.path(d, sprintf("16_block%d.R", i)))
cat(length(st), "blocks\n")

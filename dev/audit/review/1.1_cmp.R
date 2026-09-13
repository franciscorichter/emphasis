a <- readRDS("dev/audit/review/1.1_adv_rlib.rds"); b <- readRDS("dev/audit/review/1.1_adv_rlib-wave1.rds")
B <- merge(a$battery, b$battery, by = c("link","mb","tree","r","rho"), suffixes = c(".pre",".post"))
cat("battery rows:", nrow(B), " finite logf pre:", sum(is.finite(B$logf.pre)), " identical logf:", sum(identical(B$logf.pre, B$logf.post)), 
    " all logf identical bitwise:", isTRUE(all.equal(B$logf.pre, B$logf.post, tolerance = 0)),
    " logg identical:", isTRUE(all.equal(B$logg.pre, B$logg.post, tolerance = 0)), "\n")
print(B[!is.finite(B$logf.pre) | !is.finite(B$logf.post) | B$logf.pre != B$logf.post, ])
cat("\n--- degenerate cases (pre | post) ---\n")
for (n in names(a$zero)) cat(sprintf("%-28s pre = %-28s post = %s\n", n, paste(format(a$zero[[n]], digits = 8), collapse = " "), paste(format(b$zero[[n]], digits = 8), collapse = " ")))

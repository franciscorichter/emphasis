b <- readRDS("dev/pendant_before.rds")
a <- readRDS("dev/pendant_after.rds")
stopifnot(identical(names(b), names(a)))
for (nm in names(b)) {
  cat(sprintf("%-12s identical: %-5s  max|diff|: %s\n", nm,
              identical(b[[nm]], a[[nm]]),
              format(max(abs(unlist(b[[nm]][sapply(b[[nm]], is.numeric)]) -
                             unlist(a[[nm]][sapply(a[[nm]], is.numeric)])),
                         na.rm = TRUE), digits = 3)))
}
cat("\nALL IDENTICAL:", identical(b, a), "\n")
cat("fit_cr before:", format(b$fit_cr$pars, digits = 12), " loglik", format(b$fit_cr$loglik, digits = 12), "\n")
cat("fit_cr after :", format(a$fit_cr$pars, digits = 12), " loglik", format(a$fit_cr$loglik, digits = 12), "\n")
cat("fit_dd before:", format(b$fit_dd$pars, digits = 12), " loglik", format(b$fit_dd$loglik, digits = 12), "\n")
cat("fit_dd after :", format(a$fit_dd$pars, digits = 12), " loglik", format(a$fit_dd$loglik, digits = 12), "\n")
cat("m_dd   before:", format(b$m_dd, digits = 12), "\n")
cat("m_dd   after :", format(a$m_dd, digits = 12), "\n")

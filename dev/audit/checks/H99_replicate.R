# H99 replication: vary what the verifier did not.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({ library(emphasis); library(ape) })
pkg <- "/Users/pancho/Code/emphasis"
TIP <- 10e10
g <- function(b0, eta) b0 * exp(-(eta - 1)^2 / 2)
ok <- function(lbl, a, b, tol = 1e-9) { d <- max(abs(a - b)); cat(sprintf("  %-60s |diff| = %.2e %s\n", lbl, d, if (d < tol) "PASS" else "FAIL")) }

cat("== A'. how dead is the C++ path in the live test suite? ==\n")
for (f in list.files(file.path(pkg, "tests/testthat"), pattern = "\\.R$", full.names = TRUE)) {
  l <- readLines(f, warn = FALSE)
  cat(sprintf("  %-22s test_that: %2d  skip(): %2d  link-mentions: %d\n", basename(f),
              sum(grepl("^test_that\\(", l)), sum(grepl("^\\s*skip\\(", l)), sum(grepl("link", l))))
}
cat("  default link in estimate_rates / emphasis_pipeline / auto_bounds / simulate_tree:\n")
for (fn in c("estimate_rates", "emphasis_pipeline", "auto_bounds", "simulate_tree"))
  cat("   ", fn, "->", deparse(formals(get(fn, asNamespace("emphasis")))$link), "\n")

cat("\n== B'. different tree: 5-tip hand-built, dd gaussian, vs README formula ==\n")
bt <- c(5, 3.2, 2.1, 0.7)          # backward branching times (crown 5)
df5 <- data.frame(brts = rev(5 - bt)[-1], n = 2:4, t_ext = TIP, pd = 0)  # forward times of internal events
df5 <- rbind(df5, data.frame(brts = 5, n = 5, t_ext = TIP, pd = 0))
df5$pd <- cumsum(c(df5$n[1] * df5$brts[1], diff(df5$brts) * df5$n[-1]))  # any pd; M unused (beta_M=0)
p5 <- c(0.7, 0.2, 0, 0, 0.3, -0.15, 0, 0)
lam <- function(n) g(p5[1], p5[2] * n); mu <- function(n) g(p5[5], p5[6] * n)
dt <- diff(c(0, df5$brts))
ref5 <- sum(log(lam(df5$n[-nrow(df5)]))) - sum(dt * df5$n * (lam(df5$n) + mu(df5$n)))
lf5 <- emphasis:::eval_logf(p5, list(df5), c(1L, 0L, 0L), 2L)$logf
ok("B' 5-tip dd gaussian vs README formula", lf5, ref5)
# D-path fallback: model c(1,0,1) with beta_D = gamma_D = 0 must equal the no-D path
lfD <- emphasis:::eval_logf(p5, list(df5), c(1L, 0L, 1L), 2L)$logf
ok("B'' gaussian D-path with beta_D=gamma_D=0 == no-D path", lfD, lf5)
# peak: over beta_N grid on 2-tip tree with n=2, lambda max at beta_N=0.5 (already shown); here check n=5 row: peak at beta_N = 0.2
df1 <- data.frame(brts = 2, n = 5, t_ext = TIP, pd = 0)
bN <- seq(0, 0.5, by = 0.01)
sc <- sapply(bN, function(b) emphasis:::eval_logf(c(1, b, 0, 0, 0, 0, 0, 0), list(df1), c(1L, 0L, 0L), 2L)$logf)
cat(sprintf("  B''' n=5: argmax lambda over beta_N = %.2f (expected 0.20), lambda = %.4f (beta_0 = 1)\n", bN[which.min(sc)], -min(sc) / (5 * 2)))

cat("\n== C'. the other sampler: CEM + gaussian, and a dd gaussian MCEM smoke ==\n")
set.seed(11)
tr <- NULL
while (is.null(tr)) { x <- ape::rlineage(0.6, 0.15, Tmax = 5); x <- tryCatch(ape::drop.fossil(x), error = function(e) NULL)
  if (!is.null(x) && Ntip(x) >= 10 && Ntip(x) <= 35) tr <- x }
cat("  tips:", Ntip(tr), "\n")
t0 <- Sys.time()
fit_cem <- tryCatch(estimate_rates(tr, method = "cem", model = "cr", link = "gaussian",
                    control = list(lower_bound = c(0.05, 0.001), upper_bound = c(3, 1.5), max_iter = 3L, max_time = 60,
                                   num_points = 12L, sample_size = 3L, maxN = 200L, num_threads = 1L)),
                    error = function(e) { cat("  CEM ERROR:", conditionMessage(e), "\n"); NULL })
cat(sprintf("  CEM gaussian: %.1fs; pars = %s; loglik = %s\n", as.numeric(Sys.time() - t0, units = "secs"),
            if (is.null(fit_cem)) "NA" else paste(format(fit_cem$pars, digits = 3), collapse = " "),
            if (is.null(fit_cem)) "NA" else format(fit_cem$loglik, digits = 4)))
t0 <- Sys.time()
fit_dd <- tryCatch(estimate_rates(tr, method = "mcem", model = "dd", link = "gaussian",
                    control = list(lower_bound = c(0.05, 0, 0.001, 0), upper_bound = c(3, 0.2, 1.5, 0.2), max_iter = 3L,
                                   max_time = 60, sample_size = 20L, maxN = 300L, num_threads = 1L)),
                    error = function(e) { cat("  MCEM dd ERROR:", conditionMessage(e), "\n"); NULL })
cat(sprintf("  MCEM dd gaussian: %.1fs; pars = %s; loglik = %s; iters = %s\n", as.numeric(Sys.time() - t0, units = "secs"),
            if (is.null(fit_dd)) "NA" else paste(format(fit_dd$pars, digits = 3), collapse = " "),
            if (is.null(fit_dd)) "NA" else format(fit_dd$loglik, digits = 4),
            if (is.null(fit_dd)) "NA" else fit_dd$details$iterations))

cat("\n== E'. auto_bounds gaussian on a different tree (bird.orders, dd) and exponential for contrast ==\n")
data(bird.orders)
for (lk in c("gaussian", "exponential", "linear")) {
  ab <- tryCatch(suppressMessages(auto_bounds(bird.orders, model = "dd", link = lk, verbose = FALSE)), error = function(e) NULL)
  if (!is.null(ab)) cat(sprintf("  %-12s lower = %s | upper = %s\n", lk, paste(format(ab$lower_bound, digits = 3), collapse = " "),
                                paste(format(ab$upper_bound, digits = 3), collapse = " ")))
}
cat("done\n")

.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages(library(emphasis))
g <- function(b0, eta) b0 * exp(-(eta - 1)^2 / 2)
TIP <- 10e10
bt <- c(5, 3.2, 2.1, 0.7)
df5 <- data.frame(brts = c(rev(5 - bt)[-1], 5), n = 2:5, t_ext = TIP, pd = 0)
p5 <- c(0.7, 0.2, 0, 0, 0.3, -0.15, 0, 0)
lam <- g(p5[1], p5[2] * df5$n); mu <- g(p5[5], p5[6] * df5$n); dt <- diff(c(0, df5$brts))
a <- emphasis:::eval_logf(p5, list(df5), c(1L, 0L, 0L), 2L)$logf
b <- emphasis:::eval_logf(p5, list(df5), c(1L, 0L, 1L), 2L)$logf
cat(sprintf("hand-built: D - noD = %.4f ; 2*sum(dt*(lam+mu)) = %.4f ; full compensator sum(dt*n*(lam+mu)) = %.4f\n", b - a, 2 * sum(dt * (lam + mu)), sum(dt * df5$n * (lam + mu))))
# rows alive throughout each segment under the C++ test (s.brts <= prev && s.t_ext >= brts_i), excluding self
alive <- sapply(seq_len(nrow(df5)), function(i) { prev <- if (i == 1) 0 else df5$brts[i - 1]; sum(df5$brts[-i] <= prev & df5$t_ext[-i] >= df5$brts[i]) })
print(cbind(df5[, 1:2], alive))

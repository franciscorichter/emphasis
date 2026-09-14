a <- readRDS("dev/stopping-rule/final-out.rds")
b <- readRDS("dev/stopping-rule/nogrow-out.rds")
out <- rbind(a, b)
out <- out[out$ok, ]
out$key <- paste(out$config, out$cell)
sm <- do.call(rbind, lapply(split(out, list(out$key, out$rule), drop = TRUE),
  function(d) data.frame(key = d$key[1], rule = d$rule[1], n = nrow(d),
    conv = mean(d$converged), med_dell = median(d$delta_ell, na.rm = TRUE),
    within01 = mean(d$delta_ell > -0.1, na.rm = TRUE),
    med_iter = median(d$iterations), draws = sum(d$draws),
    med_S = median(d$S_final), secs = round(sum(d$secs)))))
print(sm[order(sm$key, sm$rule), ], row.names = FALSE, digits = 4)
cat("\n--- totals, all 10 cells ---\n")
print(do.call(rbind, lapply(split(out, out$rule), function(d) data.frame(
  rule = d$rule[1], n = nrow(d), conv = mean(d$converged),
  med_dell = median(d$delta_ell, na.rm = TRUE),
  within01 = mean(d$delta_ell > -0.1, na.rm = TRUE),
  med_iter = median(d$iterations), draws = sum(d$draws),
  secs = round(sum(d$secs))))), row.names = FALSE, digits = 4)
cat("\n--- totals, excluding the IS-collapse cells (both rules > 1 nat out) ---\n")
bad <- unique(out$tree_id[out$delta_ell < -1 & out$rule == "step"])
o2 <- out[!(out$tree_id %in% bad), ]
print(do.call(rbind, lapply(split(o2, o2$rule), function(d) data.frame(
  rule = d$rule[1], n = nrow(d), conv = mean(d$converged),
  med_dell = median(d$delta_ell, na.rm = TRUE),
  within01 = mean(d$delta_ell > -0.1, na.rm = TRUE),
  med_iter = median(d$iterations), draws = sum(d$draws),
  secs = round(sum(d$secs))))), row.names = FALSE, digits = 4)
cat("\n--- stop reasons ---\n")
print(table(out$rule, out$stop_reason))

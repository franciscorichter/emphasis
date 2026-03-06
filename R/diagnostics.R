# --------------------------------------------------------------------------- #
#  Diagnostics for the Monte Carlo Cross-Entropy Method                        #
# --------------------------------------------------------------------------- #

#' Diagnostics for a CEM fit
#'
#' Computes convergence, importance-sampling quality, and population-diversity
#' diagnostics from the output of \code{\link{emphasis_cem}} or
#' \code{\link{estimate_rates}} (method \code{"cem"}).
#'
#' The function works on the \code{details} component of an
#' \code{emphasis_fit} object or directly on the list returned by
#' \code{emphasis_cem}.  It requires that the fit was produced with the
#' current version of the package (which stores \code{history},
#' \code{final_pop}, and \code{best_IS} in the output).
#'
#' @section Convergence (\code{$convergence}):
#' A data frame with one row per iteration containing:
#' \describe{
#'   \item{\code{iteration}}{Iteration index.}
#'   \item{\code{best_loglik}}{IS log-likelihood of the best particle.}
#'   \item{\code{n_valid}}{Number of particles with a finite \code{fhat}.}
#'   \item{\code{rej_lambda}}{Trees rejected due to \code{max_lambda} exceedance.}
#'   \item{\code{rej_overruns}}{Trees rejected due to \code{max_missing}
#'     exceedance.}
#' }
#'
#' @section Population spread (\code{$population_spread}):
#' Per-iteration summary of the \code{fhat} distribution across all valid
#' particles: median, interquartile range, and total range.  A narrowing
#' range indicates population convergence.
#'
#' @section IS quality (\code{$IS_quality}):
#' Importance-sampling diagnostics computed from the final evaluation of the
#' globally best parameter vector (with \code{max(sample_size, n_boot, 1)}
#' augmented trees):
#' \describe{
#'   \item{\code{ESS}}{Effective sample size:
#'     \eqn{(\sum w_i)^2 / \sum w_i^2}.}
#'   \item{\code{ESS_fraction}}{ESS as a fraction of \code{n_trees}.}
#'   \item{\code{mean_lw}, \code{sd_lw}}{Mean and SD of IS log-weights.}
#'   \item{\code{fhat}}{IS log-likelihood estimate.}
#'   \item{\code{n_trees}, \code{n_rejected}}{Trees simulated and rejected.}
#' }
#'
#' @section Final population (\code{$final_pop}):
#' The full particle population at the last completed iteration (before
#' resampling), with columns for parameters and \code{fhat} per particle.
#' Useful for inspecting the distribution of surviving particles.
#'
#' @section Raw IS data (\code{$best_IS}):
#' The raw IS output for the best particle: \code{logf}, \code{log_q},
#' \code{lw} (log-weights), and \code{trees} (augmented tree data frames in
#' C++ format).
#'
#' @param x An \code{emphasis_fit} object (from \code{\link{estimate_rates}}
#'   with \code{method = "cem"}) or the raw list from \code{\link{emphasis_cem}}.
#' @param plot Logical; if \code{TRUE} (default), produce three diagnostic
#'   plots: convergence, population spread, and IS weight distribution.
#' @return A list (invisibly) with components \code{convergence},
#'   \code{population_spread}, \code{IS_quality}, \code{final_pop}, and
#'   \code{best_IS}.
#' @seealso \code{\link{emphasis_cem}}, \code{\link{estimate_rates}}
#' @examples
#' \dontrun{
#' sim <- simulate_tree(c(0.5, 0.1), max_t = 5, model = "cr")
#' fit <- estimate_rates(sim, method = "cem", model = "cr",
#'   lower_bound = c(0, 0), upper_bound = c(2, 1),
#'   control = list(max_iter = 15, n_boot = 100))
#' diag <- diagnose_cem(fit)
#' diag$IS_quality$ESS          # effective sample size
#' diag$IS_quality$ESS_fraction # ESS / n_trees
#' head(diag$final_pop)         # last-iteration particle cloud
#' }
#' @export
diagnose_cem <- function(x, plot = TRUE) {

  # Accept either emphasis_fit or raw emphasis_cem output
  details <- if (inherits(x, "emphasis_fit")) x$details else x

  if (is.null(details$history))
    stop("No history found in this fit object. ",
         "Re-run with the current version of the package (emphasis_cem).")

  n_iter <- length(details$best_loglik)

  # ── 1. Convergence table ──────────────────────────────────────────────────
  convergence <- data.frame(
    iteration    = seq_len(n_iter),
    best_loglik  = details$best_loglik,
    n_valid      = details$history$n_valid,
    rej_lambda   = details$history$rej_lambda,
    rej_overruns = details$history$rej_overruns
  )

  # ── 2. Population spread per iteration ───────────────────────────────────
  pop_spread <- do.call(rbind, lapply(seq_len(n_iter), function(k) {
    fv <- details$history$fhat_all[[k]]
    fv <- fv[is.finite(fv)]
    if (length(fv) == 0L) {
      return(data.frame(iteration = k, median = NA_real_,
                        iqr = NA_real_, range = NA_real_))
    }
    data.frame(
      iteration = k,
      median    = stats::median(fv),
      iqr       = stats::IQR(fv),
      range     = diff(range(fv))
    )
  }))

  # ── 3. IS quality at best particle ───────────────────────────────────────
  IS_quality <- NULL
  if (!is.null(details$best_IS)) {
    bi <- details$best_IS
    lw_fin <- bi$lw[is.finite(bi$lw)]
    IS_quality <- list(
      ESS          = bi$ESS,
      ESS_fraction = if (!is.na(bi$ESS)) bi$ESS / bi$n_trees else NA_real_,
      mean_lw      = if (length(lw_fin) > 0) mean(lw_fin)   else NA_real_,
      sd_lw        = if (length(lw_fin) > 1) stats::sd(lw_fin) else NA_real_,
      fhat         = bi$fhat,
      n_trees      = bi$n_trees,
      n_rejected   = bi$n_rejected
    )
  }

  # ── 4. Final population data frame ───────────────────────────────────────
  final_pop_df <- NULL
  if (!is.null(details$final_pop)) {
    fp <- details$final_pop
    final_pop_df <- cbind(fp$pars, fhat = fp$fhat)
  }

  # ── 5. Plots ──────────────────────────────────────────────────────────────
  if (plot) {
    n_plots  <- 2L + (!is.null(IS_quality))
    old_par  <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par))
    graphics::par(mfrow = c(1L, n_plots), mar = c(4, 4, 3, 1))

    # Plot 1: convergence
    stop_label <- if (!is.null(details$converged)) details$converged else "?"
    graphics::plot(convergence$iteration, convergence$best_loglik,
                   type = "b", pch = 16, col = "steelblue",
                   xlab = "Iteration", ylab = expression(hat(f)(theta^"*")),
                   main = paste0("CEM convergence\n(stop: ", stop_label, ")"))

    # Plot 2: population spread (median ± IQR band)
    valid_rows <- !is.na(pop_spread$median)
    if (any(valid_rows)) {
      ps <- pop_spread[valid_rows, ]
      graphics::plot(ps$iteration, ps$median,
                     type = "b", pch = 16, col = "darkorange",
                     ylim = range(c(ps$median - ps$iqr / 2,
                                    ps$median + ps$iqr / 2), na.rm = TRUE),
                     xlab = "Iteration",
                     ylab = expression(hat(f) ~ "(valid particles)"),
                     main = "Population spread\n(median \u00b1 IQR)")
      graphics::polygon(
        c(ps$iteration, rev(ps$iteration)),
        c(ps$median - ps$iqr / 2, rev(ps$median + ps$iqr / 2)),
        col  = grDevices::adjustcolor("darkorange", alpha.f = 0.25),
        border = NA
      )
    }

    # Plot 3: IS log-weight distribution at best particle
    if (!is.null(IS_quality) && !is.null(details$best_IS)) {
      lw_fin <- details$best_IS$lw[is.finite(details$best_IS$lw)]
      if (length(lw_fin) > 1L) {
        graphics::hist(lw_fin, breaks = min(30L, length(lw_fin)),
                       freq = FALSE, col = "grey85", border = "white",
                       xlab = "IS log-weight  (logf \u2212 log_q)",
                       main = sprintf(
                         "IS weight distribution\nESS = %.1f / %d  (%.0f%%)",
                         IS_quality$ESS, IS_quality$n_trees,
                         100 * IS_quality$ESS_fraction))
        graphics::abline(v = IS_quality$mean_lw, col = "red", lty = 2, lwd = 2)
        graphics::legend("topright", legend = "mean lw",
                         col = "red", lty = 2, lwd = 2, bty = "n")
      }
    }
  }

  invisible(list(
    convergence      = convergence,
    population_spread = pop_spread,
    IS_quality       = IS_quality,
    final_pop        = final_pop_df,
    best_IS          = details$best_IS
  ))
}


#' Print method for CEM diagnostics
#'
#' @param x Output of \code{\link{diagnose_cem}}.
#' @param ... Additional arguments (ignored).
#' @export
print.cem_diagnostics <- function(x, ...) {
  cat("CEM diagnostics\n")
  cat("===============\n\n")

  cat(sprintf("Iterations run:  %d\n", nrow(x$convergence)))
  cat(sprintf("Best fhat:       %.4f\n", max(x$convergence$best_loglik, na.rm = TRUE)))

  valid_last <- utils::tail(x$convergence$n_valid, 1L)
  cat(sprintf("Valid particles (last iter): %d\n\n", valid_last))

  if (!is.null(x$IS_quality)) {
    iq <- x$IS_quality
    cat("IS quality at best particle:\n")
    cat(sprintf("  n_trees:      %d  (rejected: %d)\n",
                iq$n_trees, iq$n_rejected))
    cat(sprintf("  ESS:          %.1f  (%.0f%% of n_trees)\n",
                iq$ESS, 100 * iq$ESS_fraction))
    cat(sprintf("  mean lw:      %.4f\n", iq$mean_lw))
    cat(sprintf("  sd lw:        %.4f\n", iq$sd_lw))
    cat(sprintf("  fhat:         %.4f\n", iq$fhat))
  }

  invisible(x)
}

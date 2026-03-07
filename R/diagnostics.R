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
#' @param plot Logical; if \code{TRUE} (default), produce diagnostic plots.
#' @param lower_bound,upper_bound Optional numeric vectors of the same length as
#'   the parameter vector.  When provided, dashed horizontal lines are drawn at
#'   the lower (orange) and upper (purple) bounds in each parameter-trace panel,
#'   making it easy to see whether the estimate has converged away from the
#'   boundary.
#' @param true_pars Optional named or unnamed numeric vector of true parameter
#'   values.  Useful in simulation studies where the generative parameters are
#'   known.  When provided, a red dashed line is added to each parameter-trace
#'   panel at the true value.  In real-data analyses this argument should be
#'   omitted (the true parameters are unknown).
#' @return A list of class \code{"cem_diagnostics"} (invisibly) with components
#'   \code{convergence}, \code{population_spread}, \code{IS_quality},
#'   \code{final_pop}, and \code{best_IS}.
#' @seealso \code{\link{emphasis_cem}}, \code{\link{estimate_rates}}
#' @examples
#' \dontrun{
#' lb  <- c(0.01, -0.5, 0, -0.5)
#' ub  <- c(2,     0.5, 1,  0.5)
#' sim <- simulate_tree(c(0.6, -0.05, 0.15, -0.03), max_t = 10, model = "pd")
#' fit <- estimate_rates(sim, method = "cem", model = "pd",
#'   lower_bound = lb, upper_bound = ub,
#'   control = list(max_iter = 50, num_points = 80, n_boot = 100))
#'
#' # Basic diagnostics -- bounds shown as dashed lines in parameter panels
#' diag <- diagnose_cem(fit, lower_bound = lb, upper_bound = ub)
#'
#' # Simulation study: also show the true generating parameters
#' diag <- diagnose_cem(fit, lower_bound = lb, upper_bound = ub,
#'                      true_pars = c(0.6, -0.05, 0.15, -0.03))
#'
#' diag$IS_quality$ESS_fraction  # IS efficiency at the best particle
#' head(diag$final_pop)          # last-iteration particle cloud
#' diag$convergence              # per-iteration best fhat and rejection counts
#' }
#' @export
diagnose_cem <- function(x, plot = TRUE,
                          lower_bound = NULL,
                          upper_bound = NULL,
                          true_pars   = NULL) {

  # Accept either emphasis_fit or raw emphasis_cem output
  details <- if (inherits(x, "emphasis_fit")) x$details else x

  if (is.null(details$history))
    stop("No history found in this fit object. ",
         "Re-run with the current version of the package (emphasis_cem).")

  n_iter <- length(details$best_loglik)

  # -- 1. Convergence table --------------------------------------------------
  convergence <- data.frame(
    iteration    = seq_len(n_iter),
    best_loglik  = details$best_loglik,
    n_valid      = details$history$n_valid,
    rej_lambda   = details$history$rej_lambda,
    rej_overruns = details$history$rej_overruns
  )
  convergence$rej_total <- convergence$rej_lambda + convergence$rej_overruns

  # -- 2. Population spread per iteration -----------------------------------
  pop_spread <- do.call(rbind, lapply(seq_len(n_iter), function(k) {
    fv <- details$history$fhat_all[[k]]
    fv <- fv[is.finite(fv)]
    if (length(fv) == 0L) {
      return(data.frame(iteration = k, median = NA_real_,
                        q25 = NA_real_, q75 = NA_real_, range = NA_real_))
    }
    data.frame(
      iteration = k,
      median    = stats::median(fv),
      q25       = stats::quantile(fv, 0.25),
      q75       = stats::quantile(fv, 0.75),
      range     = diff(range(fv))
    )
  }))

  # -- 3. IS quality at best particle ---------------------------------------
  IS_quality <- NULL
  if (!is.null(details$best_IS)) {
    bi     <- details$best_IS
    lw_fin <- bi$lw[is.finite(bi$lw)]
    IS_quality <- list(
      ESS          = bi$ESS,
      ESS_fraction = if (!is.na(bi$ESS)) bi$ESS / bi$n_trees else NA_real_,
      mean_lw      = if (length(lw_fin) > 0) mean(lw_fin)      else NA_real_,
      sd_lw        = if (length(lw_fin) > 1) stats::sd(lw_fin) else NA_real_,
      fhat         = bi$fhat,
      n_trees      = bi$n_trees,
      n_rejected   = bi$n_rejected
    )
  }

  # -- 4. Final population data frame ---------------------------------------
  final_pop_df <- NULL
  if (!is.null(details$final_pop)) {
    fp <- details$final_pop
    final_pop_df <- cbind(fp$pars, fhat = fp$fhat)
  }

  # -- 5. Parameter names ----------------------------------------------------
  par_names <- if (inherits(x, "emphasis_fit") && !is.null(names(x$pars))) {
    names(x$pars)
  } else if (!is.null(details$best_pars) && !is.null(colnames(details$best_pars))) {
    colnames(details$best_pars)
  } else if (!is.null(details$best_pars)) {
    paste0("par", seq_len(ncol(details$best_pars)))
  } else {
    character(0)
  }
  n_par <- length(par_names)

  # -- 6. Plots --------------------------------------------------------------
  if (plot) {
    has_IS  <- !is.null(IS_quality) && !is.null(details$best_IS) &&
               length(details$best_IS$lw[is.finite(details$best_IS$lw)]) > 1L
    has_rej <- any(convergence$rej_total > 0L)

    n_plots <- 2L + has_IS + has_rej + n_par
    nc      <- min(n_plots, 3L)
    nr      <- ceiling(n_plots / nc)

    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par))
    graphics::par(mfrow = c(nr, nc), mar = c(4, 4, 3, 1))

    iter <- convergence$iteration

    # Plot 1: log-likelihood trace
    stop_label <- if (!is.null(details$converged)) details$converged else "?"
    best_k <- which.max(convergence$best_loglik)
    graphics::plot(iter, convergence$best_loglik,
                   type = "b", pch = 16, col = "steelblue", lwd = 1.5,
                   xlab = "Iteration",
                   ylab = expression(hat(f)^"*" ~ "(best particle)"),
                   main = paste0("Log-likelihood convergence\nstop: ", stop_label))
    graphics::abline(v = best_k, col = "red", lty = 3, lwd = 1)

    # Plots 2..(1+n_par): per-parameter traces with optional bound/true lines
    if (n_par > 0L && !is.null(details$best_pars)) {
      bp <- details$best_pars
      for (j in seq_len(n_par)) {
        nm  <- par_names[j]
        lb_j <- if (!is.null(lower_bound) && j <= length(lower_bound))
                  lower_bound[j] else NULL
        ub_j <- if (!is.null(upper_bound) && j <= length(upper_bound))
                  upper_bound[j] else NULL
        tp_j <- if (!is.null(true_pars)) {
          if (!is.null(names(true_pars)) && nm %in% names(true_pars))
            true_pars[[nm]]
          else if (is.numeric(true_pars) && j <= length(true_pars))
            true_pars[[j]]
          else NULL
        } else NULL

        # y-range: include bounds and true value so lines are always visible
        y_vals <- c(bp[, j], lb_j, ub_j, tp_j)
        y_rng  <- range(y_vals, na.rm = TRUE)
        y_pad  <- max(diff(y_rng) * 0.08, 0.05)

        graphics::plot(iter, bp[, j],
                       type = "b", pch = 16, col = "forestgreen", lwd = 1.5,
                       ylim = c(y_rng[1] - y_pad, y_rng[2] + y_pad),
                       xlab = "Iteration", ylab = nm, main = nm)

        leg_lab <- character(0); leg_col <- character(0)
        leg_lty <- integer(0);   leg_lwd <- numeric(0)

        if (!is.null(lb_j)) {
          graphics::abline(h = lb_j, col = "darkorange", lty = 2, lwd = 1.5)
          leg_lab <- c(leg_lab, "lower bound")
          leg_col <- c(leg_col, "darkorange")
          leg_lty <- c(leg_lty, 2L); leg_lwd <- c(leg_lwd, 1.5)
        }
        if (!is.null(ub_j)) {
          graphics::abline(h = ub_j, col = "purple", lty = 2, lwd = 1.5)
          leg_lab <- c(leg_lab, "upper bound")
          leg_col <- c(leg_col, "purple")
          leg_lty <- c(leg_lty, 2L); leg_lwd <- c(leg_lwd, 1.5)
        }
        if (!is.null(tp_j)) {
          graphics::abline(h = tp_j, col = "red", lty = 2, lwd = 2)
          leg_lab <- c(leg_lab, "true value")
          leg_col <- c(leg_col, "red")
          leg_lty <- c(leg_lty, 2L); leg_lwd <- c(leg_lwd, 2)
        }
        if (length(leg_lab) > 0L)
          graphics::legend("topright", legend = leg_lab, col = leg_col,
                           lty = leg_lty, lwd = leg_lwd, bty = "n", cex = 0.75)
      }
    }

    # Population spread plot (median + IQR band)
    valid_rows <- !is.na(pop_spread$median)
    if (any(valid_rows)) {
      ps    <- pop_spread[valid_rows, ]
      y_rng <- range(c(ps$q25, ps$q75), na.rm = TRUE)
      graphics::plot(ps$iteration, ps$median,
                     type = "b", pch = 16, col = "darkorange", lwd = 1.5,
                     ylim = y_rng,
                     xlab = "Iteration",
                     ylab = expression(hat(f) ~ "(valid particles)"),
                     main = "Population spread (median + IQR)")
      graphics::polygon(
        c(ps$iteration, rev(ps$iteration)),
        c(ps$q25, rev(ps$q75)),
        col    = grDevices::adjustcolor("darkorange", alpha.f = 0.25),
        border = NA
      )
    }

    # Rejection rate plot (only when rejections occurred)
    if (has_rej) {
      graphics::plot(iter, convergence$rej_total,
                     type = "b", pch = 16, col = "tomato", lwd = 1.5,
                     xlab = "Iteration", ylab = "Rejected trees",
                     main = "Tree rejections per iteration")
      if (any(convergence$rej_lambda > 0L))
        graphics::lines(iter, convergence$rej_lambda,
                        col = "darkorange", lty = 2, lwd = 1.2)
      if (any(convergence$rej_overruns > 0L))
        graphics::lines(iter, convergence$rej_overruns,
                        col = "steelblue", lty = 2, lwd = 1.2)
      graphics::legend("topright",
                       legend = c("total", "max_lambda", "max_missing"),
                       col = c("tomato", "darkorange", "steelblue"),
                       lty = c(1L, 2L, 2L), lwd = c(1.5, 1.2, 1.2),
                       bty = "n", cex = 0.75)
    }

    # IS log-weight histogram at best particle
    if (has_IS) {
      lw_fin <- details$best_IS$lw[is.finite(details$best_IS$lw)]
      graphics::hist(lw_fin,
                     breaks = min(30L, max(10L, length(lw_fin))),
                     freq = FALSE, col = "grey85", border = "white",
                     xlab = expression("IS log-weight" ~ (log*f - log*q)),
                     main = sprintf(
                       "IS weights at best particle\nESS = %.1f / %d  (%.0f%%)",
                       IS_quality$ESS, IS_quality$n_trees,
                       100 * IS_quality$ESS_fraction))
      graphics::abline(v = IS_quality$mean_lw, col = "red", lty = 2, lwd = 2)
      graphics::legend("topright", "mean lw",
                       col = "red", lty = 2, lwd = 2, bty = "n", cex = 0.75)
    }
  }

  structure(
    list(
      convergence       = convergence,
      population_spread = pop_spread,
      IS_quality        = IS_quality,
      final_pop         = final_pop_df,
      best_IS           = details$best_IS
    ),
    class = c("cem_diagnostics", "list")
  )
}


# --------------------------------------------------------------------------- #
#  Diagnostics for MCEM                                                        #
# --------------------------------------------------------------------------- #

#' Diagnostics for an MCEM fit
#'
#' Computes convergence and importance-sampling quality diagnostics from
#' the output of \code{\link{estimate_rates}} with \code{method = "mcem"}.
#'
#' @section Convergence (\code{$convergence}):
#' A data frame with one row per EM iteration containing:
#' \describe{
#'   \item{\code{iteration}}{Iteration index.}
#'   \item{\code{fhat}}{IS log-likelihood at the M-step estimate.}
#'   \item{\code{delta_max}}{Maximum relative parameter change (fraction of
#'     search range) from the previous iteration.}
#'   \item{\code{num_trees}}{Number of augmented trees used.}
#'   \item{\code{time}}{Wall-clock time for the iteration (ms).}
#' }
#'
#' @section IS quality (\code{$IS_quality}):
#' Importance-sampling diagnostics from the final EM iteration:
#' \describe{
#'   \item{\code{ESS}}{Effective sample size.}
#'   \item{\code{ESS_fraction}}{ESS as a fraction of \code{num_trees}.}
#'   \item{\code{mean_lw}, \code{sd_lw}}{Mean and SD of IS log-weights.}
#'   \item{\code{fhat}}{IS log-likelihood at the final parameters.}
#'   \item{\code{loglik_var}}{Bootstrap variance of \code{fhat}.}
#' }
#'
#' @param x An \code{emphasis_fit} object (from \code{\link{estimate_rates}}
#'   with \code{method = "mcem"}) or the raw list from \code{.mcem_dynamic_fresh}.
#' @param plot Logical; if \code{TRUE} (default), produce diagnostic plots.
#' @param true_pars Optional numeric vector of true parameter values
#'   (for simulation studies).  Red dashed lines are drawn at these values.
#' @param lower_bound,upper_bound Optional numeric vectors for bound lines.
#' @return A list of class \code{"mcem_diagnostics"} (invisibly) with
#'   components \code{convergence}, \code{IS_quality},
#'   and \code{final_IS}.
#' @seealso \code{\link{diagnose_cem}}, \code{\link{estimate_rates}}
#' @examples
#' \dontrun{
#' sim <- simulate_tree(c(0.5, 0.1), max_t = 5, model = "cr")
#' fit <- estimate_rates(sim, method = "mcem", model = "cr",
#'   control = list(lower_bound = c(0, 0), upper_bound = c(2, 1)))
#' diagnose_mcem(fit, lower_bound = c(0, 0), upper_bound = c(2, 1))
#' }
#' @export
diagnose_mcem <- function(x, plot = TRUE,
                           lower_bound = NULL,
                           upper_bound = NULL,
                           true_pars   = NULL) {

  # Accept either emphasis_fit or raw .mcem_dynamic_fresh output
  if (inherits(x, "emphasis_fit")) {
    details  <- x$details
    par_names <- names(x$pars)
  } else {
    details   <- x
    par_names <- NULL
  }

  if (is.null(details$mcem))
    stop("No MCEM trace found. Is this an MCEM fit?")

  mcem_df  <- details$mcem
  n_iter   <- nrow(mcem_df)

  # -- 1. Convergence table --------------------------------------------------
  par_cols <- grep("^par[0-9]+$", names(mcem_df), value = TRUE)
  convergence <- data.frame(
    iteration = seq_len(n_iter),
    fhat      = mcem_df$fhat,
    delta_max = if ("delta_max" %in% names(mcem_df)) mcem_df$delta_max
                else rep(NA_real_, n_iter),
    rejected  = if ("rejected" %in% names(mcem_df)) mcem_df$rejected
                else rep(NA_integer_, n_iter),
    num_trees = mcem_df$num_trees,
    time      = mcem_df$time
  )

  # -- 2. IS quality at final iteration --------------------------------------
  IS_quality <- NULL
  final_IS   <- details$final_IS
  if (!is.null(final_IS)) {
    lw_fin <- final_IS$lw[is.finite(final_IS$lw)]
    n_trees <- length(final_IS$lw)
    IS_quality <- list(
      ESS          = final_IS$ESS,
      ESS_fraction = if (!is.na(final_IS$ESS) && n_trees > 0)
                       final_IS$ESS / n_trees else NA_real_,
      mean_lw      = if (length(lw_fin) > 0) mean(lw_fin) else NA_real_,
      sd_lw        = if (length(lw_fin) > 1) stats::sd(lw_fin) else NA_real_,
      fhat         = final_IS$fhat,
      num_trees    = n_trees,
      loglik_var   = if (!is.null(details$loglik_var)) details$loglik_var else NA_real_
    )
  }

  # -- 3. Parameter names --------------------------------------------------
  if (is.null(par_names) && length(par_cols) > 0L)
    par_names <- paste0("par", seq_along(par_cols))
  n_par <- length(par_cols)

  # -- 4. Plots ------------------------------------------------------------
  if (plot) {
    has_IS    <- !is.null(final_IS) &&
                 length(final_IS$lw[is.finite(final_IS$lw)]) > 1L
    has_delta <- any(!is.na(convergence$delta_max))

    n_plots <- 1L + has_delta + has_IS + n_par
    nc      <- min(n_plots, 3L)
    nr      <- ceiling(n_plots / nc)

    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par))
    graphics::par(mfrow = c(nr, nc), mar = c(4, 4, 3, 1))

    iter <- seq_len(n_iter)

    # Plot 1: fhat trace
    graphics::plot(iter, mcem_df$fhat,
                   type = "b", pch = 16, col = "steelblue", lwd = 1.5,
                   xlab = "EM iteration", ylab = expression(hat(f)),
                   main = sprintf("Log-likelihood trace  (%d iterations)", n_iter))
    if (n_iter > 1L) {
      fhat_fin <- mcem_df$fhat[is.finite(mcem_df$fhat)]
      if (length(fhat_fin) > 0L)
        graphics::abline(h = utils::tail(fhat_fin, 1L),
                         col = "red", lty = 3, lwd = 1)
    }

    # Plots: per-parameter traces
    if (n_par > 0L) {
      for (j in seq_len(n_par)) {
        nm   <- par_names[j]
        vals <- mcem_df[, par_cols[j]]

        lb_j <- if (!is.null(lower_bound) && j <= length(lower_bound))
                  lower_bound[j] else NULL
        ub_j <- if (!is.null(upper_bound) && j <= length(upper_bound))
                  upper_bound[j] else NULL
        tp_j <- if (!is.null(true_pars)) {
          if (!is.null(names(true_pars)) && nm %in% names(true_pars))
            true_pars[[nm]]
          else if (is.numeric(true_pars) && j <= length(true_pars))
            true_pars[[j]]
          else NULL
        } else NULL

        y_vals <- c(vals, lb_j, ub_j, tp_j)
        y_rng  <- range(y_vals, na.rm = TRUE)
        y_pad  <- max(diff(y_rng) * 0.08, 0.05)

        graphics::plot(iter, vals,
                       type = "b", pch = 16, col = "forestgreen", lwd = 1.5,
                       ylim = c(y_rng[1] - y_pad, y_rng[2] + y_pad),
                       xlab = "EM iteration", ylab = nm, main = nm)

        leg_lab <- character(0); leg_col <- character(0)
        leg_lty <- integer(0);   leg_lwd <- numeric(0)

        if (!is.null(lb_j)) {
          graphics::abline(h = lb_j, col = "darkorange", lty = 2, lwd = 1.5)
          leg_lab <- c(leg_lab, "lower"); leg_col <- c(leg_col, "darkorange")
          leg_lty <- c(leg_lty, 2L);     leg_lwd <- c(leg_lwd, 1.5)
        }
        if (!is.null(ub_j)) {
          graphics::abline(h = ub_j, col = "purple", lty = 2, lwd = 1.5)
          leg_lab <- c(leg_lab, "upper"); leg_col <- c(leg_col, "purple")
          leg_lty <- c(leg_lty, 2L);     leg_lwd <- c(leg_lwd, 1.5)
        }
        if (!is.null(tp_j)) {
          graphics::abline(h = tp_j, col = "red", lty = 2, lwd = 2)
          leg_lab <- c(leg_lab, "true"); leg_col <- c(leg_col, "red")
          leg_lty <- c(leg_lty, 2L);    leg_lwd <- c(leg_lwd, 2)
        }
        if (length(leg_lab) > 0L)
          graphics::legend("topright", legend = leg_lab, col = leg_col,
                           lty = leg_lty, lwd = leg_lwd, bty = "n", cex = 0.75)
      }
    }

    # Parameter stability trace
    if (has_delta) {
      dm <- convergence$delta_max[!is.na(convergence$delta_max)]
      dm_iter <- convergence$iteration[!is.na(convergence$delta_max)]
      graphics::plot(dm_iter, dm,
                     type = "b", pch = 16, col = "darkorange", lwd = 1.5,
                     log = "y",
                     xlab = "EM iteration",
                     ylab = expression(Delta[max]),
                     main = "Parameter stability")
      graphics::abline(h = dm[length(dm)], col = "grey50", lty = 3, lwd = 1)
    }

    # IS log-weight histogram at final iteration
    if (has_IS) {
      lw_fin <- final_IS$lw[is.finite(final_IS$lw)]
      graphics::hist(lw_fin,
                     breaks = min(30L, max(10L, length(lw_fin))),
                     freq = FALSE, col = "grey85", border = "white",
                     xlab = expression("IS log-weight" ~ (log*f - log*q)),
                     main = sprintf(
                       "IS weights (final E-step)\nESS = %.1f / %d  (%.0f%%)",
                       IS_quality$ESS, IS_quality$num_trees,
                       100 * IS_quality$ESS_fraction))
      graphics::abline(v = IS_quality$mean_lw, col = "red", lty = 2, lwd = 2)
    }
  }

  structure(
    list(
      convergence = convergence,
      stop_reason = details$stop_reason,
      IS_quality  = IS_quality,
      final_IS    = final_IS
    ),
    class = c("mcem_diagnostics", "list")
  )
}


#' Print method for MCEM diagnostics
#'
#' @param x Output of \code{\link{diagnose_mcem}}.
#' @param ... Additional arguments (ignored).
#' @export
print.mcem_diagnostics <- function(x, ...) {
  cat("MCEM diagnostics\n")
  cat("================\n\n")

  cat(sprintf("Iterations:      %d\n", nrow(x$convergence)))
  if (!is.null(x$stop_reason))
    cat(sprintf("Stop reason:     %s\n", x$stop_reason))

  fhat_fin <- x$convergence$fhat[is.finite(x$convergence$fhat)]
  if (length(fhat_fin) > 0L)
    cat(sprintf("Final fhat:      %.4f\n", utils::tail(fhat_fin, 1L)))

  dm <- x$convergence$delta_max
  dm <- dm[!is.na(dm)]
  if (length(dm) > 0L)
    cat(sprintf("Final delta_max: %.2e\n", utils::tail(dm, 1L)))

  if (!is.null(x$IS_quality)) {
    iq <- x$IS_quality
    cat("\nIS quality (final E-step):\n")
    cat(sprintf("  num_trees:    %d\n", iq$num_trees))
    cat(sprintf("  ESS:          %.1f  (%.0f%% of num_trees)\n",
                iq$ESS, 100 * iq$ESS_fraction))
    cat(sprintf("  mean lw:      %.4f\n", iq$mean_lw))
    cat(sprintf("  sd lw:        %.4f\n", iq$sd_lw))
    cat(sprintf("  fhat:         %.4f\n", iq$fhat))
    if (!is.na(iq$loglik_var))
      cat(sprintf("  loglik_se:    %.4f  (bootstrap, B=200)\n",
                  sqrt(iq$loglik_var)))
  }

  invisible(x)
}


# --------------------------------------------------------------------------- #
#  Diagnostics for GAM-based MLE                                               #
# --------------------------------------------------------------------------- #

#' Diagnostics for a GAM-based MLE fit
#'
#' Produces diagnostic plots and returns summary statistics from a
#' \code{\link{estimate_rates}} fit obtained with \code{method = "gam"}.
#'
#' @section Surface coverage (\code{$surface_summary}):
#' Summary of the IS log-likelihood grid: total grid points, number with
#' finite fhat, fhat range, and mean trees per grid point.
#'
#' @section GAM quality (\code{$gam_summary}):
#' Key statistics from the fitted GAM: deviance explained, GCV score,
#' and effective degrees of freedom.
#'
#' @section MLE (\code{$mle}):
#' The MLE parameter estimates and predicted log-likelihood at the optimum.
#'
#' @param x An \code{emphasis_fit} object from \code{\link{estimate_rates}}
#'   with \code{method = "gam"}.
#' @param plot Logical; if \code{TRUE} (default), produce diagnostic plots.
#' @param true_pars Optional numeric vector of true parameter values
#'   (for simulation studies).
#' @return A list of class \code{"gam_diagnostics"} (invisibly) with
#'   components \code{surface_summary}, \code{gam_summary}, and \code{mle}.
#' @seealso \code{\link{diagnose_cem}}, \code{\link{diagnose_mcem}},
#'   \code{\link{estimate_rates}}
#' @examples
#' \dontrun{
#' sim <- simulate_tree(c(0.5, 0.1), max_t = 5, model = "cr")
#' fit <- estimate_rates(sim, method = "gam", model = "cr",
#'   control = list(lower_bound = c(0, 0), upper_bound = c(2, 1)))
#' diagnose_gam(fit)
#' }
#' @export
diagnose_gam <- function(x, plot = TRUE, true_pars = NULL) {

  if (inherits(x, "emphasis_fit")) {
    details <- x$details
    par_names <- names(x$pars)
  } else {
    details   <- x
    par_names <- NULL
  }

  if (is.null(details$surface) || is.null(details$gam_fit))
    stop("No GAM data found. Is this a method=\"gam\" fit?")

  surface   <- details$surface
  gam_fit   <- details$gam_fit
  if (is.null(par_names))
    par_names <- details$par_names

  n_par     <- length(par_names)
  valid     <- is.finite(surface$fhat)
  n_valid   <- sum(valid)
  n_total   <- nrow(surface)

  # -- 1. Surface summary ---------------------------------------------------
  surface_summary <- list(
    n_grid_points = n_total,
    n_valid       = n_valid,
    coverage      = n_valid / n_total,
    fhat_range    = if (n_valid > 0) range(surface$fhat[valid]) else c(NA, NA),
    mean_trees    = if (n_valid > 0) mean(surface$n_trees[valid]) else NA_real_
  )

  # -- 2. GAM summary -------------------------------------------------------
  gs <- summary(gam_fit)
  gam_summary <- list(
    deviance_explained = gs$dev.expl,
    gcv_score          = gam_fit$gcv.ubre,
    edf_total          = sum(gs$edf),
    r_squared          = gs$r.sq
  )

  # -- 3. MLE ---------------------------------------------------------------
  mle_pars   <- if (inherits(x, "emphasis_fit")) x$pars else NULL
  mle_loglik <- if (inherits(x, "emphasis_fit")) x$loglik else NULL

  # -- 4. Plots --------------------------------------------------------------
  if (plot && n_valid > 0L) {
    n_plots <- n_par + 2L  # 1D slices + residuals + predicted vs observed
    nc      <- min(n_plots, 3L)
    nr      <- ceiling(n_plots / nc)

    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par))
    graphics::par(mfrow = c(nr, nc), mar = c(4, 4, 3, 1))

    surf_valid <- surface[valid, , drop = FALSE]

    # Per-parameter: fhat vs parameter value with GAM prediction curve
    for (j in seq_len(n_par)) {
      nm <- par_names[j]
      x_vals <- surf_valid[[nm]]
      y_vals <- surf_valid$fhat

      graphics::plot(x_vals, y_vals,
                     pch = 16, col = grDevices::adjustcolor("steelblue", 0.5),
                     cex = 0.7,
                     xlab = nm, ylab = expression(hat(f)),
                     main = paste("fhat vs", nm))

      # GAM prediction curve along this parameter (others at MLE)
      x_seq <- seq(min(x_vals), max(x_vals), length.out = 100L)
      if (!is.null(mle_pars)) {
        nd <- as.data.frame(as.list(mle_pars))
        nd <- nd[rep(1L, 100L), , drop = FALSE]
        nd[[nm]] <- x_seq
        pred <- stats::predict(gam_fit, newdata = nd, type = "response")
        graphics::lines(x_seq, pred, col = "red", lwd = 2)
      }

      if (!is.null(true_pars)) {
        tp_j <- if (!is.null(names(true_pars)) && nm %in% names(true_pars))
          true_pars[[nm]]
        else if (j <= length(true_pars)) true_pars[[j]]
        else NULL
        if (!is.null(tp_j))
          graphics::abline(v = tp_j, col = "red", lty = 2, lwd = 2)
      }
      if (!is.null(mle_pars))
        graphics::abline(v = mle_pars[[nm]], col = "forestgreen", lty = 2, lwd = 1.5)
    }

    # Predicted vs observed fhat
    fitted_vals <- stats::fitted(gam_fit)
    obs_vals    <- gam_fit$y
    graphics::plot(obs_vals, fitted_vals,
                   pch = 16, col = grDevices::adjustcolor("steelblue", 0.5),
                   cex = 0.7,
                   xlab = "Observed fhat", ylab = "GAM predicted fhat",
                   main = sprintf("Predicted vs observed\n(R2=%.3f, dev.expl=%.1f%%)",
                                  gam_summary$r_squared,
                                  100 * gam_summary$deviance_explained))
    graphics::abline(0, 1, col = "red", lty = 2, lwd = 1.5)

    # Residuals
    resid_vals <- stats::residuals(gam_fit)
    graphics::hist(resid_vals,
                   breaks = min(30L, max(10L, length(resid_vals))),
                   freq = FALSE, col = "grey85", border = "white",
                   xlab = "Residual", main = "GAM residuals")
    graphics::abline(v = 0, col = "red", lty = 2, lwd = 1.5)
  }

  result <- structure(
    list(
      surface_summary = surface_summary,
      gam_summary     = gam_summary,
      mle             = list(pars = mle_pars, loglik = mle_loglik),
      surface         = surface,
      gam_fit         = gam_fit
    ),
    class = c("gam_diagnostics", "list")
  )
  invisible(result)
}


#' Print method for GAM diagnostics
#'
#' @param x Output of \code{\link{diagnose_gam}}.
#' @param ... Additional arguments (ignored).
#' @export
print.gam_diagnostics <- function(x, ...) {
  cat("GAM-based MLE diagnostics\n")
  cat("=========================\n\n")

  ss <- x$surface_summary
  cat(sprintf("Grid coverage:    %d/%d valid (%.0f%%)\n",
              ss$n_valid, ss$n_grid_points, 100 * ss$coverage))
  cat(sprintf("fhat range:       [%.2f, %.2f]\n",
              ss$fhat_range[1], ss$fhat_range[2]))
  cat(sprintf("Mean trees/point: %.0f\n\n", ss$mean_trees))

  gs <- x$gam_summary
  cat(sprintf("GAM fit:\n"))
  cat(sprintf("  R-squared:        %.4f\n", gs$r_squared))
  cat(sprintf("  Deviance expl:    %.1f%%\n", 100 * gs$deviance_explained))
  cat(sprintf("  Effective df:     %.1f\n\n", gs$edf_total))

  if (!is.null(x$mle$pars)) {
    cat("MLE estimates:\n")
    print(round(x$mle$pars, 6))
    if (!is.null(x$mle$loglik))
      cat(sprintf("\nlog-likelihood:   %.4f\n", x$mle$loglik))
  }

  invisible(x)
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

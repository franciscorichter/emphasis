# ═══════════════════════════════════════════════════════════════
# BDI exact augmentation sampler
#
# Draws augmented trees from the exact conditional distribution
# p(A | T, theta) via an aggregate Gillespie algorithm.
# Under constant rates: ESS = S (zero-variance proposal).
# Under DD/endogenous: self-consistent backward-forward iteration.
# ═══════════════════════════════════════════════════════════════


# --------------------------------------------------------------------------- #
#  Survival probability helpers                                                #
# --------------------------------------------------------------------------- #

# CR analytical solution
.bdi_p_cr <- function(t, lam0, mu0, tp) {
  E0 <- exp(-(lam0 - mu0) * (tp - t))
  (lam0 - mu0) / (lam0 - mu0 * E0)
}


# --------------------------------------------------------------------------- #
#  Analytical integrals for exact BDI logg (CR only)                           #
# --------------------------------------------------------------------------- #

#' Exact integral of total BDI rate over [t1, t2] for CR with constant n, k.
#'
#' total_rate(t) = (n+2k)*lam0*(1-p(t)) + n*mu0/(1-p(t))
#'
#' Uses the identities:
#'   int lam0*(1-p) dt = mu0*(t2-t1) + ln[p(t1)/p(t2)]
#'   int mu0/(1-p) dt  = lam0*(t2-t1) + ln[(1-p(t1))/(1-p(t2))]
#'
#' @return The integral value (may be Inf when n>0 and t2==tp).
#' @keywords internal
.bdi_integral_cr <- function(t1, t2, n, k, lam0, mu0, tp) {
  if (t2 - t1 < 1e-15) return(0)
  d <- lam0 - mu0
  if (abs(d) < 1e-15) {
    # Critical case: lam0 = mu0
    # p(t) = lam0*(tp-t)/(1+lam0*(tp-t)), 1-p = 1/(1+lam0*(tp-t))
    a1 <- 1 + lam0 * (tp - t1)
    a2 <- 1 + lam0 * (tp - t2)
    I_lam <- log(a1 / a2)
    I_mu  <- lam0 * (t2 - t1) + log(a2 / a1)  # mu/(1-p) = mu*(1+lam*(tp-t))
    return((n + 2L * k) * I_lam + n * I_mu)
  }

  E1 <- exp(-d * (tp - t1))
  E2 <- exp(-d * (tp - t2))

  # I_lam = int lam0*(1-p) = mu0*dt + ln(p1/p2)
  #       = mu0*(t2-t1) + ln[(lam0-mu0*E2)/(lam0-mu0*E1)]
  I_lam <- mu0 * (t2 - t1) + log((lam0 - mu0 * E2) / (lam0 - mu0 * E1))

  if (n > 0L) {
    oE1 <- 1 - E1
    oE2 <- 1 - E2
    if (oE2 < 1e-300) return(Inf)   # t2 at tp: divergent
    I_mu <- lam0 * (t2 - t1) + log(oE1 / oE2)
  } else {
    I_mu <- 0
  }

  (n + 2L * k) * I_lam + n * I_mu
}


#' Find next event time using exact time-change for CR.
#'
#' Given cumulative hazard H(s) = .bdi_integral_cr(t_cur, s, ...),
#' find t* such that H(t*) = U where U ~ Exp(1).
#'
#' @return Event time, or value > t_max if no event before boundary.
#' @keywords internal
.bdi_find_event_time_cr <- function(t_cur, U, n, k, lam0, mu0, tp, t_max) {
  # Quick check: will H reach U before t_max?
  # When n>0 and t_max==tp, H→Inf so event always exists.
  is_last_seg <- (tp - t_max < 1e-10) && (n > 0L)
  if (!is_last_seg) {
    H_max <- .bdi_integral_cr(t_cur, t_max, n, k, lam0, mu0, tp)
    if (U >= H_max) return(t_max + 1)  # no event
  }

  # Upper bound for root search
  t_hi <- if (is_last_seg) tp - 1e-14 else t_max
  f <- function(s) .bdi_integral_cr(t_cur, s, n, k, lam0, mu0, tp) - U
  stats::uniroot(f, c(t_cur + 1e-15, t_hi), tol = 1e-13)$root
}

#' Can the BDI sampler be used for this model/link?
#'
#' The sampler's rate functions and mean-field ODEs are written in the
#' \code{(N, P, E)} covariate layout (slots 2-4 of \code{pars8}), and its
#' link handling covers only \code{linear} (0) and \code{exponential} (1).
#' Under the package's current \code{(N, M, D)} basis and with the
#' \code{gaussian} link (2), only N-only models on links 0/1 are exact, so
#' every other case is routed to the thinning proposal by the callers.
#' @keywords internal
.bdi_supported <- function(model_bin, link) {
  model_bin <- as.integer(model_bin)
  link      <- as.integer(link)
  length(model_bin) == 3L && model_bin[2L] == 0L && model_bin[3L] == 0L &&
    link %in% c(0L, 1L)
}

# Compute speciation rate from 8-param vector + model
.bdi_lam <- function(pars8, N, P, E, model_bin, link) {
  eta <- pars8[1] + pars8[2] * N + pars8[3] * P + pars8[4] * E
  if (link == 0L) max(0, eta) else exp(eta)
}

.bdi_mu <- function(pars8, N, P, E, model_bin, link) {
  eta <- pars8[5] + pars8[6] * N + pars8[7] * P + pars8[8] * E
  if (link == 0L) max(0, eta) else exp(eta)
}


# --------------------------------------------------------------------------- #
#  Backward-forward iteration (DD / general endogenous)                        #
# --------------------------------------------------------------------------- #

#' Solve survival probability p(t) backward from tp to 0.
#' ODE: dp/dt = lam(t)*p^2 - (lam(t)-mu(t))*p, p(tp) = 1.
#' Rates use the mean-field covariates (N̂, P̂, Ê) — supports DD/PD/EP and
#' mixed models. For CR (all model_bin = 0) the covariate arguments are
#' ignored inside .bdi_lam/.bdi_mu.
#' @keywords internal
.bdi_solve_p_backward <- function(pars8, model_bin, link, bt, tp,
                                  Nhat_fun, Phat_fun, Ehat_fun, t_grid) {
  n_grid <- length(t_grid)
  t_rev  <- rev(t_grid)
  p_vals <- numeric(n_grid)
  p      <- 1.0

  for (i in seq_along(t_rev)) {
    ti <- t_rev[i]
    p_vals[n_grid - i + 1L] <- p
    if (i < length(t_rev)) {
      dt <- t_rev[i] - t_rev[i + 1L]
      f <- function(t, p) {
        Nh <- Nhat_fun(t); Ph <- Phat_fun(t); Eh <- Ehat_fun(t)
        la <- .bdi_lam(pars8, Nh, Ph, Eh, model_bin, link)
        mu <- .bdi_mu(pars8, Nh, Ph, Eh, model_bin, link)
        la * p^2 - (la - mu) * p
      }
      # RK4 backward step
      k1 <- f(ti, p)
      k2 <- f(ti - dt / 2, p - dt / 2 * k1)
      k3 <- f(ti - dt / 2, p - dt / 2 * k2)
      k4 <- f(ti - dt, p - dt * k3)
      p  <- p - dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
      p  <- min(max(p, 0), 1)
    }
  }
  p_vals
}


#' Solve the forward mean ODE jointly for (m, P_miss).
#'
#' State: s = (m, P_miss) with
#'   dm/dt       = (lam_bdi - mu_bdi) * m + nu      (self-consistent BDI mean)
#'   dP_miss/dt  = m                                (mean-field: each alive
#'                                                   missing lineage extends
#'                                                   its pendant edge at rate 1)
#'
#' Rates lam, mu are evaluated at (N̂, P̂, Ê) to support DD/PD/EP + mixed.
#' @return list(m_vals, pm_vals) on t_grid.
#' @keywords internal
.bdi_solve_mean_forward <- function(pars8, model_bin, link, bt, tp,
                                    p_fun, Nhat_fun, Phat_fun, Ehat_fun,
                                    t_grid) {
  bt <- sort(bt)
  m  <- 0; pm <- 0
  m_vals  <- numeric(length(t_grid))
  pm_vals <- numeric(length(t_grid))

  k_of_t <- function(t) 2L + sum(bt <= t)

  for (i in seq_along(t_grid)) {
    ti <- t_grid[i]
    m_vals[i]  <- m
    pm_vals[i] <- pm
    if (i < length(t_grid)) {
      dt <- t_grid[i + 1L] - t_grid[i]
      f <- function(t, state) {
        mi <- state[1L]; pmi <- state[2L]
        Nh <- Nhat_fun(t); Ph <- Phat_fun(t); Eh <- Ehat_fun(t)
        la <- .bdi_lam(pars8, Nh, Ph, Eh, model_bin, link)
        mu <- .bdi_mu(pars8, Nh, Ph, Eh, model_bin, link)
        p   <- p_fun(t); omp <- 1 - p
        if (omp < 1e-14) return(c(-mi * 1e6, mi))
        la_bdi <- la * omp
        mu_bdi <- mu / omp
        k      <- k_of_t(t)
        nu     <- 2 * k * la * omp
        c((la_bdi - mu_bdi) * mi + nu, mi)
      }
      s  <- c(m, pm)
      k1 <- f(ti,            s)
      k2 <- f(ti + dt / 2,   s + dt / 2 * k1)
      k3 <- f(ti + dt / 2,   s + dt / 2 * k2)
      k4 <- f(ti + dt,       s + dt     * k3)
      s  <- s + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
      m  <- max(s[1L], 0)
      pm <- max(s[2L], 0)
    }
  }
  list(m_vals = m_vals, pm_vals = pm_vals)
}


#' Solve the forward mean + second-moment ODE (Gaussian closure).
#'
#' Augments `.bdi_solve_mean_forward` with variance state needed to cancel
#' the process-level Jensen residual E[P/N] − P̂/N̂ ≈ −Cov(P,N)/N̂² + P̂ Var(N)/N̂³.
#' State: s = (m, P_miss, Var(N), Cov(N,P)) with
#'   dm/dt         = (lam_bdi − mu_bdi) m + nu
#'   dP_miss/dt    = m
#'   dVar(N)/dt    = 2(lam_bdi − mu_bdi) Var(N) + (lam_bdi + mu_bdi) m + nu
#'   dCov(N,P)/dt  = (lam_bdi − mu_bdi) Cov(N,P) + Var(N)
#'
#' The Var(N) ODE is the standard jump-diffusion variance infusion for a BDI
#' process; the Cov(N,P) ODE comes from the deterministic coupling dP = N dt
#' between events. Cross-jumps at deaths (ΔN = −1, ΔP = −age) are dropped
#' at this order — they are second-order corrections to cNP and would require
#' tracking the age distribution.
#' @return list(m_vals, pm_vals, vN_vals, cNP_vals) on t_grid.
#' @keywords internal
.bdi_solve_meanvar_forward <- function(pars8, model_bin, link, bt, tp,
                                       p_fun, Nhat_fun, Phat_fun, Ehat_fun,
                                       t_grid) {
  bt <- sort(bt)
  m <- 0; pm <- 0; vN <- 0; cNP <- 0
  m_vals   <- numeric(length(t_grid))
  pm_vals  <- numeric(length(t_grid))
  vN_vals  <- numeric(length(t_grid))
  cNP_vals <- numeric(length(t_grid))

  k_of_t <- function(t) 2L + sum(bt <= t)

  for (i in seq_along(t_grid)) {
    ti <- t_grid[i]
    m_vals[i]   <- m
    pm_vals[i]  <- pm
    vN_vals[i]  <- vN
    cNP_vals[i] <- cNP
    if (i < length(t_grid)) {
      dt <- t_grid[i + 1L] - t_grid[i]
      f <- function(t, state) {
        mi   <- state[1L]; pmi  <- state[2L]
        vNi  <- state[3L]; cNPi <- state[4L]
        Nh <- Nhat_fun(t); Ph <- Phat_fun(t); Eh <- Ehat_fun(t)
        la <- .bdi_lam(pars8, Nh, Ph, Eh, model_bin, link)
        mu <- .bdi_mu(pars8, Nh, Ph, Eh, model_bin, link)
        p   <- p_fun(t); omp <- 1 - p
        if (omp < 1e-14) return(c(-mi * 1e6, mi, 0, 0))
        la_bdi <- la * omp
        mu_bdi <- mu / omp
        k      <- k_of_t(t)
        nu     <- 2 * k * la * omp
        drift  <- la_bdi - mu_bdi
        infuse <- la_bdi + mu_bdi
        c(drift * mi + nu,                           # dm
          mi,                                        # dP_miss
          2 * drift * vNi + infuse * mi + nu,        # dVar(N)
          drift * cNPi + vNi)                        # dCov(N,P)
      }
      s  <- c(m, pm, vN, cNP)
      k1 <- f(ti,          s)
      k2 <- f(ti + dt / 2, s + dt / 2 * k1)
      k3 <- f(ti + dt / 2, s + dt / 2 * k2)
      k4 <- f(ti + dt,     s + dt     * k3)
      s   <- s + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
      m   <- max(s[1L], 0)
      pm  <- max(s[2L], 0)
      vN  <- max(s[3L], 0)
      cNP <- s[4L]
    }
  }
  list(m_vals = m_vals, pm_vals = pm_vals,
       vN_vals = vN_vals, cNP_vals = cNP_vals)
}


#' Observed-tree phylogenetic diversity under the emphasis tip_start=0
#' convention.
#'
#' The BDI tree data frame (.bdi_to_tree_df) sets tip_start = 0 for every
#' observed lineage, so calculate_pendant_pd(t) = k(t) * t.  Missing species
#' carry their actual birth time.  P(t) = P_obs(t) + P_miss(t).
#' @keywords internal
.bdi_P_obs <- function(t, bt) {
  k <- 2L + sum(bt <= t)
  k * t
}

#' Backward-forward iteration for self-consistent BDI rates.
#'
#' Mean-field state: N̂(t) = k(t) + m̂(t), P̂(t) = P_obs(t) + P̂_miss(t),
#' Ê(t) = P̂(t) / N̂(t).  All three enter the rate function; inactive
#' covariates (pars8 slopes pinned to 0) make the extension transparent
#' for CR/DD and activate PD/EP/mixed models.
#' @return List with p_fun, Nhat_fun, Phat_fun, Ehat_fun.
#' @keywords internal
.bdi_iterate <- function(pars8, model_bin, link, bt, tp,
                         max_iter = 20, tol = 1e-4, n_grid = 500,
                         use_gaussian_closure = TRUE) {
  bt     <- sort(bt)
  t_grid <- seq(0, tp, length.out = n_grid)

  k_of_t <- function(t) 2L + sum(bt <= t)
  k_vals     <- sapply(t_grid, k_of_t)
  P_obs_vals <- k_vals * t_grid

  # Initial guess: m = 0, P_miss = 0 ⇒ N̂ = k, P̂ = P_obs, Ê = P̂/N̂.
  Nhat_vals <- k_vals
  Phat_vals <- P_obs_vals
  Ehat_vals <- ifelse(Nhat_vals > 0, Phat_vals / Nhat_vals, 0)
  vN_vals   <- numeric(length(t_grid))
  cNP_vals  <- numeric(length(t_grid))

  # Scale for convergence check (normalize P̂ by a stable magnitude).
  Pscale <- max(max(Phat_vals), 1)

  p_vals <- rep(0, length(t_grid))

  for (iter in seq_len(max_iter)) {
    Nhat_fun <- stats::approxfun(t_grid, Nhat_vals, rule = 2)
    Phat_fun <- stats::approxfun(t_grid, Phat_vals, rule = 2)
    Ehat_fun <- stats::approxfun(t_grid, Ehat_vals, rule = 2)

    p_vals <- .bdi_solve_p_backward(pars8, model_bin, link, bt, tp,
                                    Nhat_fun, Phat_fun, Ehat_fun, t_grid)
    p_fun  <- stats::approxfun(t_grid, p_vals, rule = 2)

    if (use_gaussian_closure) {
      mean_res <- .bdi_solve_meanvar_forward(pars8, model_bin, link, bt, tp,
                                             p_fun, Nhat_fun, Phat_fun, Ehat_fun,
                                             t_grid)
      m_vals   <- mean_res$m_vals
      pm_vals  <- mean_res$pm_vals
      vN_vals  <- mean_res$vN_vals
      cNP_vals <- mean_res$cNP_vals
    } else {
      mean_res <- .bdi_solve_mean_forward(pars8, model_bin, link, bt, tp,
                                          p_fun, Nhat_fun, Phat_fun, Ehat_fun,
                                          t_grid)
      m_vals   <- mean_res$m_vals
      pm_vals  <- mean_res$pm_vals
      vN_vals  <- rep(0, length(t_grid))
      cNP_vals <- rep(0, length(t_grid))
    }

    Nhat_new <- k_vals + m_vals
    Phat_new <- P_obs_vals + pm_vals
    # Second-order Taylor correction: E[P/N] ≈ P̂/N̂ − Cov(P,N)/N̂² + P̂·Var(N)/N̂³
    # (zero under plain mean-field closure ⇒ backwards-compatible).
    Ehat_raw <- ifelse(Nhat_new > 0,
                       Phat_new / Nhat_new
                       - cNP_vals / pmax(Nhat_new, 1e-12)^2
                       + Phat_new * vN_vals / pmax(Nhat_new, 1e-12)^3,
                       0)
    # Clip to [0, tp]: E_s is a pendant age.
    Ehat_new <- pmin(pmax(Ehat_raw, 0), tp)

    delta_N <- max(abs(Nhat_new - Nhat_vals))
    delta_P <- max(abs(Phat_new - Phat_vals)) / Pscale
    delta   <- max(delta_N, delta_P)

    Nhat_vals <- Nhat_new
    Phat_vals <- Phat_new
    Ehat_vals <- Ehat_new
    Pscale    <- max(max(Phat_vals), 1)

    if (delta < tol) break
  }

  Nhat_fun <- stats::approxfun(t_grid, Nhat_vals, rule = 2)
  Phat_fun <- stats::approxfun(t_grid, Phat_vals, rule = 2)
  Ehat_fun <- stats::approxfun(t_grid, Ehat_vals, rule = 2)
  p_fun    <- stats::approxfun(t_grid, p_vals, rule = 2)

  list(p_fun    = p_fun,
       Nhat_fun = Nhat_fun,
       Phat_fun = Phat_fun,
       Ehat_fun = Ehat_fun,
       vN_vals  = vN_vals,
       cNP_vals = cNP_vals,
       t_grid   = t_grid)
}


# --------------------------------------------------------------------------- #
#  BDI Gillespie augmentation (one tree)                                       #
# --------------------------------------------------------------------------- #

#' Draw one augmented tree via the BDI process.
#'
#' Under CR: uses exact time-change method with analytical cumulative hazard.
#'   The event times are drawn from the correct inhomogeneous BDI process,
#'   giving zero-variance IS weights (logf - logg = constant).
#' Under DD: uses approximate Gillespie with piecewise-constant rates.
#'   IS weights have nonzero variance (importance sampling, not exact).
#'
#' @return List with $species, $n_alive_at_tp, $logg; or NULL on overflow.
#' @keywords internal
.bdi_augment_one <- function(bt, pars8, model_bin, link, tp,
                             p_fun = NULL, Nhat_fun = NULL,
                             Phat_fun = NULL, Ehat_fun = NULL,
                             max_missing = 1e4L) {
  bt <- sort(bt)
  is_cr <- all(model_bin == 0L)

  lam0 <- .bdi_lam(pars8, 0, 0, 0, model_bin, link)
  mu0  <- .bdi_mu(pars8, 0, 0, 0, model_bin, link)

  boundaries <- c(0, bt, tp)
  alive   <- numeric(0)
  n_alive <- 0L
  species <- list()
  n_total <- 0L
  logg    <- 0

  k_of_t <- function(t) 2L + sum(bt <= t)

  for (seg_idx in seq_len(length(boundaries) - 1L)) {
    t0 <- boundaries[seg_idx]
    t1 <- boundaries[seg_idx + 1L]
    k  <- k_of_t(t0 + 1e-12)
    t  <- t0

    while (t < t1) {
      if (is_cr) {
        # ── Exact time-change method for CR ──
        U <- stats::rexp(1)
        t_star <- .bdi_find_event_time_cr(t, U, n_alive, k, lam0, mu0, tp, t1)

        if (t_star >= t1) {
          # No event before boundary
          logg <- logg - .bdi_integral_cr(t, t1, n_alive, k, lam0, mu0, tp)
          break
        }

        # Event at t_star: survival integral contribution = -U
        logg <- logg - U
        t    <- t_star

        # Compute exact rates at event time for type selection
        p    <- .bdi_p_cr(t, lam0, mu0, tp)
        omp  <- 1 - p
        la   <- lam0 * omp
        mu_r <- mu0 / max(omp, 1e-15)
        nu   <- 2 * k * lam0 * omp
        total <- n_alive * (la + mu_r) + nu

      } else {
        # ── Approximate Gillespie for DD / PD / EP / mixed ──
        # Rates use the mean-field covariates from the backward-forward
        # iteration. Under DD only N̂ is active; PD activates P̂;
        # EP activates Ê = P̂/N̂ (mean isolation time).
        p   <- p_fun(t)
        omp <- 1 - p
        if (omp < 1e-12) break
        Nh <- Nhat_fun(t); Ph <- Phat_fun(t); Eh <- Ehat_fun(t)
        la_raw <- .bdi_lam(pars8, Nh, Ph, Eh, model_bin, link)
        mu_raw <- .bdi_mu(pars8, Nh, Ph, Eh, model_bin, link)
        la   <- la_raw * omp
        mu_r <- mu_raw / max(omp, 1e-15)
        nu   <- 2 * k * la_raw * omp
        total <- n_alive * (la + mu_r) + nu

        if (total < 1e-15) break
        dt <- stats::rexp(1, total)
        tn <- t + dt
        if (tn >= t1) {
          logg <- logg - total * (t1 - t)
          break
        }
        t <- tn
        logg <- logg - total * dt
      }

      # ── Event type selection (shared) ──
      # logg uses per-lineage rates (labeled density) to match C++ logf:
      #   logf has log(λ) per speciation = per-lineage rate
      #   so logg must use per-lineage BDI rates, not total rates.
      r <- stats::runif(1) * total
      if (r < n_alive * la) {
        logg    <- logg + log(la)           # per-species birth rate
        alive   <- c(alive, t)
        n_alive <- n_alive + 1L
        n_total <- n_total + 1L
      } else if (r < n_alive * la + nu) {
        logg    <- logg + log(la)           # per-lineage immigration rate
        alive   <- c(alive, t)
        n_alive <- n_alive + 1L
        n_total <- n_total + 1L
      } else if (n_alive > 0L) {
        logg <- logg + log(mu_r)            # per-species death rate
        idx  <- sample.int(n_alive, 1L)
        species[[length(species) + 1L]] <- c(alive[idx], t)
        alive   <- alive[-idx]
        n_alive <- n_alive - 1L
      }

      if (n_total > max_missing) return(NULL)
    }
  }

  # Under CR: exact process ensures all species die before tp.
  # Under DD: approximate Gillespie may leave survivors — reject.
  if (n_alive > 0L) return(NULL)

  list(species = species, n_alive_at_tp = 0L, logg = logg)
}


# --------------------------------------------------------------------------- #
#  Convert BDI output to emphasis tree data frame                              #
# --------------------------------------------------------------------------- #

#' Convert BDI species list to the tree data frame format used by emphasis.
#' Columns: brts, n, t_ext, pd, tip_start, id, parent_id
#' @keywords internal
.bdi_to_tree_df <- function(species, bt, tp) {
  bt_sorted <- sort(bt)
  n_obs     <- length(bt_sorted)
  n_aug     <- length(species)
  # obs + (spec + ext) per missing + closing node at tp
  n_total   <- n_obs + 2L * n_aug + 1L

  # Pre-allocate vectors
  v_brts      <- numeric(n_total)
  v_t_ext     <- numeric(n_total)
  v_tip_start <- numeric(n_total)
  v_id        <- integer(n_total)
  v_parent_id <- integer(n_total)

  # Observed speciation nodes
  idx <- seq_len(n_obs)
  v_brts[idx]      <- bt_sorted
  v_t_ext[idx]     <- 1e11   # t_ext_tip
  v_tip_start[idx] <- 0
  v_id[idx]        <- seq(0L, n_obs - 1L)
  v_parent_id[idx] <- -1L

  # Augmented species: speciation + extinction nodes
  if (n_aug > 0L) {
    off <- n_obs
    for (i in seq_along(species)) {
      sp    <- species[[i]]
      birth <- sp[1]; death <- sp[2]
      sid   <- as.integer(n_obs + i - 1L)
      parent <- max(0L, which(bt_sorted <= birth) - 1L)
      if (length(parent) == 0L) parent <- 0L else parent <- max(parent)
      j <- off + 2L * (i - 1L)
      # Speciation node
      v_brts[j + 1L]      <- birth
      v_t_ext[j + 1L]     <- death
      v_tip_start[j + 1L] <- birth
      v_id[j + 1L]        <- sid
      v_parent_id[j + 1L] <- parent
      # Extinction node
      v_brts[j + 2L]      <- death
      v_t_ext[j + 2L]     <- 0     # t_ext_extinct
      v_tip_start[j + 2L] <- birth
      v_id[j + 2L]        <- sid
      v_parent_id[j + 2L] <- parent
    }
  }

  # Closing node at tp (present-day marker, matches C++ convention)
  v_brts[n_total]      <- tp
  v_t_ext[n_total]     <- 1e11
  v_tip_start[n_total] <- 0
  v_id[n_total]        <- -1L
  v_parent_id[n_total] <- -1L

  # Sort by brts
  ord <- order(v_brts)
  v_brts      <- v_brts[ord]
  v_t_ext     <- v_t_ext[ord]
  v_tip_start <- v_tip_start[ord]
  v_id        <- v_id[ord]
  v_parent_id <- v_parent_id[ord]

  # Compute n: n[i] = lineage count during [brts[i-1], brts[i]).
  # C++ convention: n[i] = n_after(event i-1) = n[i-1] + type(event i-1).
  v_n    <- numeric(n_total)
  n_cur  <- 2
  v_n[1] <- n_cur
  if (n_total >= 2L) {
    for (i in 2L:n_total) {
      n_cur  <- n_cur + if (v_t_ext[i - 1L] == 0) -1 else 1
      v_n[i] <- n_cur
    }
  }

  data.frame(
    brts      = v_brts,
    n         = v_n,
    t_ext     = v_t_ext,
    pd        = rep(0, n_total),
    tip_start = v_tip_start,
    id        = v_id,
    parent_id = v_parent_id,
    stringsAsFactors = FALSE
  )
}


# --------------------------------------------------------------------------- #
#  Main BDI augmentation wrapper (replaces .augment_tree_internal)             #
# --------------------------------------------------------------------------- #

#' BDI augmentation: draw sample_size augmented trees.
#'
#' Returns same structure as augment_trees() C++ function:
#' list(trees, logf, logg, weights, fhat)
#'
#' @keywords internal
.augment_tree_bdi <- function(tree,
                              pars,
                              model_bin   = c(0L, 0L, 0L),
                              sample_size = 1L,
                              max_missing = 1e4,
                              link        = 0L,
                              rho         = 1.0,
                              use_gaussian_closure = TRUE) {
  brts  <- .extract_brts(tree)
  # Accept either compact or 8-element pars
  if (length(pars) == 8L) {
    pars8 <- pars
  } else {
    pars8 <- .expand_pars(pars, model_bin)
  }
  tp    <- brts[1L]

  # Convert emphasis brts (crown-age first, decreasing) to forward-time bt
  # emphasis brts: tp, t_{n-1}, ..., t_1 (decreasing, present = 0)
  # forward time bt: tp - brts (increasing from 0)
  bt <- sort(tp - brts[-1L])

  is_cr <- all(model_bin == 0L)

  # Solve BDI rates
  p_fun <- Nhat_fun <- Phat_fun <- Ehat_fun <- NULL
  if (!is_cr) {
    sol <- .bdi_iterate(pars8, model_bin, link, bt, tp,
                        use_gaussian_closure = use_gaussian_closure)
    p_fun    <- sol$p_fun
    Nhat_fun <- sol$Nhat_fun
    Phat_fun <- sol$Phat_fun
    Ehat_fun <- sol$Ehat_fun
  }

  # Draw augmented trees.
  # For DD (approximate Gillespie), some draws are rejected (survivors at tp),
  # so we oversample and collect until we have sample_size valid trees.
  trees      <- vector("list", sample_size)
  logg       <- numeric(sample_size)
  n_valid    <- 0L
  max_tries  <- if (is_cr) sample_size else 5L * sample_size

  for (attempt in seq_len(max_tries)) {
    if (n_valid >= sample_size) break
    aug <- .bdi_augment_one(bt, pars8, model_bin, link, tp,
                            p_fun, Nhat_fun, Phat_fun, Ehat_fun,
                            max_missing)
    if (is.null(aug)) next

    n_valid <- n_valid + 1L
    trees[[n_valid]] <- .bdi_to_tree_df(aug$species, bt, tp)
    logg[n_valid]     <- aug$logg
  }

  # Trim to actual count
  trees <- trees[seq_len(n_valid)]
  logg  <- logg[seq_len(n_valid)]

  # Compute logf (model log-likelihood) via C++ eval_logf.
  # eval_logf also returns a thinning-based logg — we discard it
  # and use the Gillespie-accumulated logg from above instead.
  if (length(trees) > 0L) {
    ev   <- eval_logf(pars8, trees,
                      model = as.integer(model_bin),
                      link  = as.integer(link),
                      rho   = as.numeric(rho))
    logf <- ev$logf
  } else {
    logf <- numeric(0)
  }

  weights <- logf - logg
  max_lw  <- if (length(weights) > 0) max(weights) else 0
  sum_w   <- sum(exp(weights - max_lw))
  fhat    <- if (length(weights) > 0) log(sum_w / length(weights)) + max_lw else -Inf

  list(trees   = trees,
       logf    = logf,
       logg    = logg,
       weights = weights,
       fhat    = fhat)
}


# --------------------------------------------------------------------------- #
#  BDI MCEM loop: E-step (BDI) + M-step (C++ nlopt)                           #
# --------------------------------------------------------------------------- #

#' MCEM with BDI augmentation.
#'
#' Drop-in replacement for .mcem_dynamic_fresh(). Uses BDI exact sampling
#' for the E-step and m_cpp() for the M-step.
#'
#' Under CR: zero-variance IS, so every E-step gives the same fhat regardless
#' of the augmented trees.  Convergence depends only on the M-step.
#'
#' Stopping rule.  With theta_k the M-step output of iteration k,
#' \deqn{\delta_k = \max_j |\theta_k[j] - \theta_{k-1}[j]| / \max(|\theta_{k-1}[j]|, \epsilon)}
#' with \eqn{\epsilon = 10^{-2}}, and the run stops as "converged" once
#' \eqn{\delta_k < tol} for \code{patience} consecutive iterations.  The
#' scale is the parameter itself, not the search box, so the same fit stops
#' at the same point whatever bounds the user passes.  An iteration whose
#' M-step returned its starting point unchanged carries no information about
#' stability and resets the streak.  The trace also records the absolute
#' step \code{abs_step} and \code{drift}, the same relative displacement
#' measured over the last \code{patience} iterations: under pure Monte Carlo
#' noise \code{drift} is of the order of \code{sqrt(patience)} steps, under a
#' deterministic EM drift it is \code{patience} steps.
#'
#' The reported likelihood is a separate E-step at the final iterate, so
#' \code{fhat}, ESS, \code{final_IS} and \code{loglik_var} all describe
#' \code{pars}; that E-step is the last row of \code{mcem} (\code{m_step =
#' FALSE}, no step columns) and is not counted in \code{iterations}.
#'
#' @return A list: \code{mcem} (trace, one row per iteration plus the final
#'   E-step), \code{pars}, \code{iterations} (number of completed E+M
#'   iterations), \code{stop_reason}, \code{loglik} (fhat at \code{pars},
#'   \code{NA} when the final E-step failed), \code{loglik_var},
#'   \code{final_IS}, \code{n_failed} (E- or M-step failures).
#' @keywords internal
.mcem_bdi <- function(brts,
                      pars,
                      sample_size,
                      max_missing,
                      lower_bound,
                      upper_bound,
                      max_iter,
                      xtol,
                      tol = 1e-2,
                      patience,
                      num_threads,
                      verbose    = FALSE,
                      conditional = NULL,
                      model      = c(0L, 0L, 0L),
                      link       = 0L,
                      max_time   = NULL,
                      rho        = 1.0) {

  if (inherits(brts, "phylo"))
    brts <- sort(ape::branching.times(brts), decreasing = TRUE)
  if (!is.numeric(brts)) stop("`brts` must be numeric or a `phylo` object.")
  if (!is.numeric(pars)) stop("`pars` must be numeric.")

  model_bin <- as.integer(model)
  link_int  <- as.integer(link)
  patience  <- max(1L, as.integer(patience))

  # Floor of the relative-change denominator: a parameter below eps in
  # magnitude is measured against eps, so a coordinate sitting at zero does
  # not turn every step into an infinite relative change.
  eps <- 1e-2
  rel_change <- function(new, old)
    max(abs(new - old) / pmax(abs(old), eps))

  # One BDI E-step at `theta`; NULL on error or when no tree was drawn.
  e_step_at <- function(theta) {
    e <- tryCatch(
      .augment_tree_bdi(
        tree        = brts,
        pars        = theta,
        model_bin   = model_bin,
        sample_size = as.integer(sample_size),
        max_missing = as.integer(max_missing),
        link        = link_int,
        rho         = as.numeric(rho)
      ),
      error = function(e) NULL
    )
    if (is.null(e) || length(e$trees) == 0L) return(NULL)
    e
  }

  # Trace row.  Step columns are NA for the final E-step (no M-step).
  trace_row <- function(theta, e, elapsed_e, m_time = 0, m_step = TRUE,
                        delta_max = NA_real_, abs_step = NA_real_,
                        drift = NA_real_, m_moved = NA) {
    par_df <- as.data.frame(as.list(
      stats::setNames(theta, paste0("par", seq_along(theta)))))
    cbind(par_df, data.frame(
      fhat        = e$fhat,
      delta_max   = delta_max,
      abs_step    = abs_step,
      drift       = drift,
      m_step      = m_step,
      m_moved     = m_moved,
      rejected    = .n0(e$n_rejected),
      n_nonfinite = sum(!is.finite(e$logf)),
      num_trees   = length(e$trees),
      ESS         = .ess_from_lw(e$weights),
      time        = elapsed_e * 1000 + m_time
    ))
  }

  streak      <- 0L
  fail_streak <- 0L
  n_failed    <- 0L
  prev_pars   <- pars          # last iterate whose E-step succeeded
  history     <- list(pars)    # theta_0, theta_1, ... for the drift column
  mcem        <- NULL
  n_iter      <- 0L
  stop_reason <- "max_iter"
  t0_mcem     <- proc.time()[3]

  for (i in seq_len(max_iter)) {
    # ── E-step: BDI augmentation ──
    t0_e  <- proc.time()[3]
    e_raw <- e_step_at(pars)
    elapsed_e <- as.numeric(proc.time()[3] - t0_e)

    if (is.null(e_raw)) {
      fail_streak <- fail_streak + 1L
      n_failed    <- n_failed + 1L
      streak      <- 0L
      # Restart from the last iterate whose E-step succeeded; the box centre
      # is not a known-good point and may itself be where the sampler fails.
      pars <- prev_pars
      if (verbose) message(sprintf(
        "Iteration %d: E-step failed (%d consecutive) - restarting from the last successful iterate",
        i, fail_streak))
      if (fail_streak >= 8L) { stop_reason <- "e_step_failure"; break }
      next
    }

    # M-step set: draws with a finite log-weight.  Non-finite draws stay in
    # the E-step's fhat denominator (handled by .augment_tree_bdi) and are
    # counted in the trace; they carry no information for the M-step.
    lw     <- e_raw$weights
    finite <- is.finite(lw) & is.finite(e_raw$logf)
    if (!any(finite)) {
      fail_streak <- fail_streak + 1L
      n_failed    <- n_failed + 1L
      streak      <- 0L
      pars <- prev_pars
      if (verbose) message(sprintf(
        "Iteration %d: E-step returned no finite log-weight (%d consecutive)",
        i, fail_streak))
      if (fail_streak >= 8L) { stop_reason <- "e_step_failure"; break }
      next
    }
    fail_streak <- 0L
    prev_pars   <- pars

    # m_cpp objective: sum loglik(theta, tree_i) * w[i] with w as direct
    # multipliers.  BDI log-weights are constant and negative under CR, so
    # convert to self-normalised IS weights (mean 1, all positive).
    lw_f   <- lw[finite]
    w_norm <- exp(lw_f - max(lw_f))
    w_norm <- w_norm / sum(w_norm) * length(w_norm)

    e_step <- list(
      trees                 = e_raw$trees[finite],
      weights               = w_norm,
      rejected              = .n0(e_raw$n_rejected),
      rejected_overruns     = 0L,
      rejected_lambda       = 0L,
      rejected_zero_weights = sum(!finite),
      time                  = elapsed_e * 1000,
      fhat                  = e_raw$fhat
    )

    # ── M-step: C++ optimisation ──
    m_result <- tryCatch(
      m_cpp(e_step     = e_step,
            init_pars  = pars,
            plugin     = "rpd1",
            lower_bound = lower_bound,
            upper_bound = upper_bound,
            xtol_rel   = xtol,
            num_threads = as.integer(num_threads),
            model      = model_bin,
            link       = link_int,
            rho        = as.numeric(rho),
            rconditional = conditional),
      error = function(e) NULL
    )

    if (is.null(m_result)) {
      fail_streak <- fail_streak + 1L
      n_failed    <- n_failed + 1L
      streak      <- 0L
      if (verbose) message(sprintf("Iteration %d: M-step failed", i))
      if (fail_streak >= 8L) { stop_reason <- "m_step_failure"; break }
      next
    }

    new_pars  <- as.numeric(m_result$estimates)
    m_moved   <- any(new_pars != pars)
    abs_step  <- max(abs(new_pars - pars))
    delta_max <- rel_change(new_pars, pars)
    n_iter    <- n_iter + 1L
    history[[n_iter + 1L]] <- new_pars
    drift <- if (n_iter > patience)
      rel_change(new_pars, history[[n_iter + 1L - patience]]) else NA_real_
    pars      <- new_pars

    mcem <- rbind(mcem, trace_row(
      pars, e_raw, elapsed_e, m_time = m_result$time,
      delta_max = delta_max, abs_step = abs_step, drift = drift,
      m_moved = m_moved))

    if (verbose) {
      message(sprintf(
        "Iteration %d: fhat=%.4f  delta=%.2e  step=%.2e  streak=%d/%d  ESS=%.0f/%d%s",
        i, e_raw$fhat, delta_max, abs_step, streak, patience,
        .ess_from_lw(e_raw$weights), length(e_raw$trees),
        if (m_moved) "" else "  (M-step returned its start)"))
    }

    # Check convergence.  An M-step that returned its start unchanged says
    # nothing about stability (the objective may have been undefined there),
    # so it does not count toward patience.
    if (!m_moved) {
      streak <- 0L
    } else if (delta_max < tol) {
      streak <- streak + 1L
      if (streak >= patience) { stop_reason <- "converged"; break }
    } else {
      streak <- 0L
    }

    # Time budget
    if (!is.null(max_time)) {
      elapsed <- proc.time()[3] - t0_mcem
      if (elapsed > max_time) {
        stop_reason <- "time_budget"
        if (verbose) message(sprintf("Time budget reached (%.0fs)", elapsed))
        break
      }
    }
  }

  # ── Final E-step at the returned iterate ──
  # The trace rows pair theta_k with fhat(theta_{k-1}); this E-step gives
  # fhat, ESS and the IS diagnostics at pars itself.
  t0_e    <- proc.time()[3]
  e_final <- e_step_at(pars)
  elapsed_e <- as.numeric(proc.time()[3] - t0_e)
  if (!is.null(e_final)) {
    mcem <- rbind(mcem, trace_row(pars, e_final, elapsed_e, m_step = FALSE))
  } else if (verbose) {
    message("Final E-step at the returned parameters failed; loglik is NA.")
  }

  loglik     <- NA_real_
  loglik_var <- NA_real_
  final_IS   <- NULL
  if (!is.null(e_final)) {
    loglik <- e_final$fhat
    lw     <- e_final$logf - e_final$logg
    finite <- is.finite(e_final$logf)
    if (sum(finite) >= 2L) {
      loglik_var <- .bootstrap_fhat_var(e_final$logf[finite], e_final$logg[finite],
                                        K = 2L, B = 200L)
    }
    final_IS <- list(
      logf  = e_final$logf,
      logg  = e_final$logg,
      lw    = lw,
      fhat  = e_final$fhat,
      ESS   = .ess_from_lw(lw),
      n_rejected = .n0(e_final$n_rejected),
      rejected_zero_weights = sum(!finite)
    )
  }

  list(
    mcem        = mcem,
    pars        = pars,
    iterations  = n_iter,
    stop_reason = stop_reason,
    loglik      = loglik,
    loglik_var  = loglik_var,
    final_IS    = final_IS,
    n_failed    = n_failed
  )
}

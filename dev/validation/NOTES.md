# Validation study — reference-side facts (established 2026-09-13, `ref_check.R`)

Tree: `TreeSim::sim.bd.taxa(n = 30, lambda = 0.6, mu = 0.2, complete = FALSE)`, seed 7, crown age 11.528.

| Reference | MLE (λ, μ) | Note |
|---|---|---|
| `ape::birthdeath` | (0.8225, 0.7416) | Nee et al. 1994, conditioned on both crown lineages surviving |
| `DDD::bd_ML(cond = 1, soc = 2)` | (0.8225, 0.7416) | identical to `ape::birthdeath` |
| `DDD::bd_ML(cond = 0, soc = 2)` | (0.5436, 0.2758) | unconditioned crown likelihood — the surface emphasis targets with `cond = NULL` (map §5.1) |

- `btorph` (0/1) shifts the log-likelihood by a θ-independent constant; MLEs unchanged.
- `DDD:::lambdamu`, `ddmodel = 1`: `λ(N) = pmax(0, λ0 − (λ0 − μ0)·N/K)`, `μ = μ0` — the same truncation as emphasis's linear link, so the mapping `β0 = λ0, βN = −(λ0−μ0)/K, γ0 = μ0, γN = 0` is exact.
- **Closed (2026-09-13).** The apparent theta-dependent offset between `dd_loglik(K -> inf)` and
  `bd_loglik` was a `pars2` slot error in this file: DDD 5.2.4 takes
  `bd_loglik: c(tdmodel, cond, btorph, verbose, soc)` and
  `dd_loglik: c(lx, ddmodel, cond, btorph, verbose, soc)`, and the earlier calls ended
  `(verbose, soc) = (2, 0)`, so the "reference" was a stem-age likelihood. With the documented
  layouts, over 4 trees x 5 theta: `max|dd_loglik(K = 1e6) - bd_loglik| = 3.4e-04`,
  `max|bd_loglik(btorph = 1) - Nee closed form| = 8.5e-14`, and
  `max|fhat_BDI - bd_loglik| = 1.6e-11` — the constant between emphasis's `f` and DDD at
  `btorph = 1` is zero. `dd_loglik(ddmodel = 1)` returns `-Inf` for `mu0 >= lam0`, so the DD
  reference exists only on `lam0 > mu0`.
- DDD's `bd_loglik`/`dd_loglik` print "Parameters: … Loglikelihood: …" unconditionally; wrap calls in `capture.output()`.
- `timeout` is not on macOS; use Rscript's own `setTimeLimit()` or the emphasis `max_time` control.

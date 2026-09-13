# Validation study — reference-side facts (established 2026-09-13, `ref_check.R`)

Tree: `TreeSim::sim.bd.taxa(n = 30, lambda = 0.6, mu = 0.2, complete = FALSE)`, seed 7, crown age 11.528.

| Reference | MLE (λ, μ) | Note |
|---|---|---|
| `ape::birthdeath` | (0.8225, 0.7416) | Nee et al. 1994, conditioned on both crown lineages surviving |
| `DDD::bd_ML(cond = 1, soc = 2)` | (0.8225, 0.7416) | identical to `ape::birthdeath` |
| `DDD::bd_ML(cond = 0, soc = 2)` | (0.5436, 0.2758) | unconditioned crown likelihood — the surface emphasis targets with `cond = NULL` (map §5.1) |

- `btorph` (0/1) shifts the log-likelihood by a θ-independent constant; MLEs unchanged.
- `DDD:::lambdamu`, `ddmodel = 1`: `λ(N) = pmax(0, λ0 − (λ0 − μ0)·N/K)`, `μ = μ0` — the same truncation as emphasis's linear link, so the mapping `β0 = λ0, βN = −(λ0−μ0)/K, γ0 = μ0, γN = 0` is exact.
- **Open:** `dd_loglik(K = 10⁶, ddmodel = 1, cond = 0)` − `bd_loglik(cond = 0)` is **not** constant in θ (5.42, 4.66, 4.02, 2.90 at four points; unchanged by `lx` 500→2000 and by `K` 10⁴→10⁶). Before `dd_ML` is used as the DD reference, find which `pars2` encoding of the two functions describes the same likelihood (check `?bd_loglik` / `?dd_loglik` positional meanings of `pars2`, and `cond` semantics), or establish the DD reference independently (e.g. `dd_ML` vs the generating parameters on `dd_sim` trees).
- DDD's `bd_loglik`/`dd_loglik` print "Parameters: … Loglikelihood: …" unconditionally; wrap calls in `capture.output()`.
- `timeout` is not on macOS; use Rscript's own `setTimeLimit()` or the emphasis `max_time` control.

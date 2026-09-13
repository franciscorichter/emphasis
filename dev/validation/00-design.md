# emphasis validation study — design

Written 2026-09-13. Runs against a post-wave-1 build (see §11). Every number in
§4 and §10 was measured today on this Mac (10 cores, R 4.4.0, DDD 5.2.4, ape
5.8.1, TreeSim 2.4) against the pre-fix scratch build at
`/Users/pancho/.claude/jobs/867af780/tmp/rlib`; the scripts that produced them
are `R/00-selfcheck.R` and `scratch/dd_ref_cost.R`.

---

## 1. The question and the estimand

The question is how well `estimate_rates(method = "mcem")` recovers **the exact
maximum-likelihood estimate of the tree in front of it**, in the constant-rates
model and in the diversity-dependent model.

The estimand is therefore the **per-tree exact MLE** `θ_MLE`, not the
generating parameters. The generating parameters are a second, separate
quantity: the exact MLE misses them too, by an amount that is a property of the
data and not of emphasis. Both are reported; they are never pooled.

Primary scalar: the **log-likelihood deficit**

```
Δℓ = ℓ_exact(θ̂_E) − ℓ_exact(θ_MLE)     ≤ 0
```

`Δℓ` is invariant to reparameterisation, is defined when the MLE sits on the
`μ = 0` boundary (which happens on a large fraction of small trees), and has a
scale: 0.5 nats is about one standard error in one parameter.

Secondary scalars: the per-parameter error in SE units `z_j = (θ̂_j − θ_MLE,j)/SE_j`,
and the reported-log-likelihood error `e = fit$loglik − ℓ_exact(θ̂_E)`.

Three error sources are separated by design rather than pooled:

| source | how it is isolated |
|---|---|
| Monte Carlo | replicate fits of the same tree (seeded for BDI; independent processes for thinning, which cannot be seeded) |
| optimisation / stopping | fits started **at** the exact MLE vs started far from it, paired per tree; plus `N` over a 20-fold range |
| MLE sampling error | `θ_MLE − θ_gen`, reported beside `θ̂_E − θ_gen`, never mixed in |

### Aims

- **A1 (cr, primary)** distribution of `θ̂_E − θ_MLE` and of `Δℓ` for
  `method = "mcem"` alone, by tree size `n`, turnover `ε = μ/λ`, sampler, Monte
  Carlo sample size `N`, and starting point; decomposed into Monte Carlo and
  systematic components.
- **A2 (cr, primary)** accuracy of the reported log-likelihood,
  `e = fit$loglik − ℓ_exact(θ̂_E)`. The constant between emphasis's `f` and DDD's
  `btorph = 1` convention is **measured to be zero** (§4.3), so no calibration
  constant enters `e`.
- **A3 (dd, primary)** A1 and A2 against `DDD::dd_ML` / `dd_loglik`
  (`ddmodel = 1`, `γ_N` fixed at 0).
- **A4 (secondary)** `θ̂_E − θ_gen` beside `θ_MLE − θ_gen`.
- **A5 (secondary, descriptive)** the full `emphasis_pipeline`
  (`bounds → gam → cem → mcem`) against `bd_ML(cond = 0)` and `bd_ML(cond = 1)`,
  with a record of whether `auto_bounds`'s box contains either.

---

## 2. Scope

- Models: `cr` and `dd` with `link = "linear"` and `γ_N` pinned at 0 by
  `lower_bound = upper_bound = 0`, which is exactly DDD's `ddmodel = 1`.
- `rho = 1`, `cond = NULL`, `num_threads = 1` everywhere.
- Samplers: `sampling = "bdi"` (default) and `sampling = "dynamic_fresh"`
  (thinning).
- `n ≤ 100` in the main tier, `n = 200` in an optional extension tier.

Out of scope, and why, in §12.

---

## 3. Data-generating mechanisms

### 3.1 CR

The CR likelihood is scale invariant: multiplying the branching times by `c` and
dividing the rates by `c` leaves every log-likelihood **difference** unchanged.
`00-selfcheck.R` assertion (g) verifies this to 0.00e+00 on a real tree. So
`λ_gen = 1` fixes the time unit and only turnover varies.

```
TreeSim::sim.bd.taxa(n = n, numbsim = 1, lambda = 1, mu = eps, complete = FALSE)
```

| factor | levels |
|---|---|
| `n` | 20, 50, 100 (main); 200 (extension) |
| `ε = μ/λ` | 0, 0.3, 0.6, 0.9 |
| trees per cell | 20 at `n` = 20 and 50, 12 at `n` = 100; set by the smoke-tier calibration rule `trees_per_cell = min(20, floor(cell_budget / measured_cost_per_tree))`, floor 10 |

Plus one **time-scale cell**: `λ_gen = 0.1`, `ε = 0.5`, `n = 50`, 10 trees. Same
process in different units; it probes the underflow cut (H1) and any residual
scale dependence of the stopping rule (H21).

`sim.bd.taxa` conditions on `n`, so the crown age `T` is random and is recorded.
That conditioning affects A4 only, and is handled there by also computing
`bd_ML(cond = 2)`, the matching conditional.

Seeds: `seed = digest::digest2int(tree_id)`, `tree_id` = `cr-n050-e06-l100-t07`.
Trees, seeds and generator calls are written to `data/trees-<tier>.rds` before
any fit runs.

### 3.2 DD

```
DDD::dd_sim(pars = c(lambda0, mu0, K), age = age, ddmodel = 1)$tes
```

| regime | `(λ0, μ0, K)` | age | note |
|---|---|---|---|
| A | (0.8, 0.1, 40) | 10 | strong DD, low turnover |
| B | (0.8, 0.2, 60) | 12 | strong DD, moderate turnover |
| C | (0.8, 0.4, 40) | 12 | high turnover; H58-exposed stress cell, reported separately |

A draw is accepted only when `n ∈ [0.5 K, min(1.2 K, 80)]`. Both bounds matter:
below `0.5 K` the clade has not approached the carrying capacity and `K` is not
identifiable, which makes `dd_ML` a boundary search rather than a reference;
above 80 the reference itself becomes the cost driver (§10). Rejected draws are
counted and the seed advances. 10 trees per regime in the main tier.

---

## 4. Exact references

All references use `btorph = 1`, `soc = 2` (crown). Every `DDD` call is wrapped
in `capture.output()` because `bd_loglik`/`dd_loglik` print unconditionally.

### 4.1 The `pars2` layout — what the open question in NOTES.md was

`dev/validation/NOTES.md` records an open discrepancy: `dd_loglik(K → ∞)` minus
`bd_loglik` was not constant in θ (5.42, 4.66, 4.02, 2.90 at four points). It is
a slot error in `ref_check.R`, not a property of DDD.

DDD 5.2.4 puts `verbose` **before** `soc`:

```
bd_loglik  pars2 = c(tdmodel, cond, btorph, verbose, soc)
dd_loglik  pars2 = c(lx, ddmodel, cond, btorph, verbose, soc)
```

`ref_check.R` lines 25 and 28 pass `c(0, cond, 1, 2, 0)` and
`c(500, 1, 0, 1, 2, 0)` — that is `verbose = 2` and `soc = 0` in both calls.
`soc = 0` is neither stem nor crown, and the two functions treat it differently.
Re-running that layout today reproduces the θ-dependent gap (1.26 to 6.03 nats
across 4 trees × 4 θ points); this is assertion (a2), which the study keeps as a
**guard**: the wrong layout must still differ, or the test is not testing
anything.

With the documented layout, measured today across trees of 10, 30, 30 and 60
tips and five θ points with `λ > μ`:

```
max | dd_loglik(pars1 = c(λ, μ, 1e6), pars2 = c(lx, 1, 0, 1, 0, 2))
    − bd_loglik(pars1 = c(λ, μ, 0, 0), pars2 = c(0, 0, 1, 0, 2)) |  =  3.38e-04
```

`lx` at `10(n+1)` versus `30(n+1)`, and `methode = "analytical"` versus
`"odeint::runge_kutta_fehlberg78"`, move `dd_loglik` by at most 6.23e-04
(assertion (e)). That 1e-3 tolerance is why the dd decision thresholds are
doubled (§8).

**`ref_check.R` and `NOTES.md` still carry the wrong layout.** This study does
not edit them; the correction is for the package owner.

### 4.2 CR

Closed form, per-lineage (labelled-history) crown likelihood, implemented in
`R/00-common.R::ll_cr_nee` and used as the reference surface:

```
ℓ(λ, μ; brts) = 2 log p1(T) + Σ_{i=2}^{n-1} [ log λ + log p1(t_i) ]
p1(t) = r² e^{−rt} / (λ − μ e^{−rt})²,   r = λ − μ
```

Two facts about this formula the implementation has to respect, both found by
running the check:

- For `r < 0` the denominator is negative for all `t > 0` (its zero lies at
  `t < 0`), so the squared form is well defined; the code takes `log|·|`. A
  guard that rejects a negative denominator makes the reference `-Inf` at
  `μ > λ`, which the study visits.
- At `λ = μ` the expression is 0/0; the limit is `p1(t) = 1/(1 + λt)²`.
- For `r < 0`, `exp(−rt)` overflows to `Inf` at moderate `|r|·T` — it does so on
  `n = 200`, `ε = 0.9` trees — and the naive expression then returns `−Inf` over
  a region the optimiser walks through. The implementation factors the growing
  term out, `log|den| = log μ − rt + log|1 − (λ/μ)e^{rt}|`, where the `−rt`
  cancels against the numerator's.

The optimiser side needs two guards for the same reason: the objective is capped
rather than `+Inf` (L-BFGS-B errors on a non-finite value), and each start is
tried with L-BFGS-B and with Nelder-Mead, best kept, then polished. A
first-order check is recorded per tree: `max_j |∂ℓ/∂θ_j| · SE_j < 1e-3`, so the
criterion is "one SE of movement changes `ℓ` by less than 1e-3 nats" and is
comparable across `n` and across time units. On the boundary only the `λ`
component is required, with `∂ℓ/∂μ < 0`.

`bd_ML` itself returns `pars = (−1, −1)` with `conv = −1` when its optimiser
fails, which happened on an `n = 200`, `ε = 0.9` tree. That is recorded as
`ddd_status = "bd_ML_failed"`: the closed form is the reference surface — tier 0
asserts it equals `bd_loglik` — so the `optim` optimum stands and the DDD
cross-check is reported as unavailable rather than as a disagreement. A tree
where **both** optimisers fail is flagged `mle_failed` and gets no fit jobs.

Measured agreement with DDD, over 4 trees × 7 θ points including `μ > λ` and
`λ = μ`:

```
max | bd_loglik(pars2 = c(0, 0, 1, 0, 2)) − ll_cr_nee |  =  8.53e-14
```

MLEs:

- `θ_MLE0` (the target of `cond = NULL`): `optim(method = "L-BFGS-B", lower = c(1e-6, 0))`
  on `ll_cr_nee` from three starts, best kept, cross-checked against
  ```
  DDD::bd_ML(brts, cond = 0, btorph = 1, soc = 2, tdmodel = 0,
             initparsopt = c(λ_gen, max(μ_gen, 0.05)), idparsopt = 1:2)
  ```
  Disagreement above 1e-4 in either parameter flags the tree, and each optimiser
  is re-run from the other's optimum.
- `θ_MLE1 = bd_ML(cond = 1, btorph = 1, soc = 2)`, which equals
  `ape::birthdeath` — measured today at (0.8225, 0.7416) on both, agreement
  < 1e-4 (assertion (d)). This is the survival-conditioned MLE and the target
  the pipeline's GAM stage approximates.
- `θ_MLE2 = bd_ML(cond = 2)`, conditioned on `n`, the conditional that matches
  `sim.bd.taxa`. Used in A4 only.
- `SE = sqrt(diag(solve(−H)))` from `numDeriv::hessian` on the closed form. When
  `μ̂ = 0` the profile SE of `λ` is used and `z_μ` is not defined; those trees
  are summarised by `Δℓ` and by `λ` alone.

### 4.3 The constant between emphasis and DDD is zero

Measured, not assumed. On 4 trees × 5 θ points with `μ < λ`, at `N = 50`:

```
max | fhat_BDI(θ) − bd_loglik(θ; btorph = 1) |  =  1.60e-11
max sd(lw)                                      =  1.56e-11
```

(assertions (c1), (c1b)). `log((n−1)!)` — 71.257 at `n = 30` — is the constant
between `btorph = 0` and `btorph = 1`, not between emphasis and DDD. The wave-1.3
test pin in `dev/audit/02-findings.md` reads `−log((n−1)!)`; that is against
`btorph = 0` and must not be copied into this study. The study keeps a per-tree
assertion `|c_n| < 1e-8` and aborts if it ever fails, rather than carrying a
calibration constant into `e`.

### 4.4 DD

```
DDD::dd_ML(brts, initparsopt = c(λ0, μ0, K)_gen, idparsopt = 1:3, ddmodel = 1,
           cond = 0, btorph = 1, soc = 2, res = max(300, 10(n+1)),
           methode = "analytical", optimmethod = "subplex")
```

run from two starts (generating parameters, and 1.5×), higher log-likelihood
kept. `conv != 0` or `K̂ > 1e4` sets the flag **`K_unidentified`**: the tree is
counted, reported, and analysed as `cr`, never scored as estimator error.

Surface: `dd_loglik(pars1 = c(λ0, μ0, K), pars2 = c(lx, 1, 0, 1, 0, 2))` with
`lx = max(300, 10(n+1))`.

Mapping, exact because `DDD:::lambdamu` for `ddmodel = 1` is
`pmax(0, la − (la − mu)·n/K)`, the same truncation as emphasis's linear link:

```
β0 = λ0 ,  β_N = −(λ0 − μ0)/K ,  γ0 = μ0 ,  γ_N = 0
```

Evaluating the exact likelihood at an emphasis estimate inverts it, with three
cases that the study records rather than silently coerces:

| emphasis estimate | reference |
|---|---|
| `β_N < 0` and `β0 > γ0` | `K = −(β0 − γ0)/β_N`, then `dd_loglik` |
| `β_N ≥ 0` (upper bound) | CR limit: `bd_loglik(cond = 0, btorph = 1, soc = 2)` |
| `β0 ≤ γ0` with `β_N < 0` | `outside_ddd`: nothing is evaluated, the fit is counted under that label |

The third case is not hypothetical. Measured today: **`dd_loglik(ddmodel = 1)`
returns `-Inf` whenever `μ0 ≥ λ0`** (checked at (0.4, 0.6) and (0.5, 0.5),
`K = 1e4` and `1e6`). The DD reference exists only on `λ0 > μ0`, which is exactly
the region where the mapping is defined.

`lx` guard: the largest augmented `N` in every E-step is recorded, and the
reference is valid only while it stays below `lx/2`. Violations are flagged,
never silently accepted.

### 4.5 The tier-0 self-check (`R/00-selfcheck.R`)

Runs first, serially, ~60 s. Aborts on (a), (b), (c1), (c1b), (d), (e), (g);
records three sentinels that the study is designed to *measure* rather than
require:

| sentinel | pins | pre-fix result (measured today) |
|---|---|---|
| c2 | H12 / wave 1.3 — BDI at `μ ≥ λ` | **FAIL**: 0 of 8 draws finite |
| f | H7 / wave 1.8 — thinning envelope | borderline: the sampler is clock-seeded, so this is 3 replicates and the criterion is the mean against its own replicate SE; single draws on the same tree ranged from −0.06 to +0.13 at ESS 39–482 |
| h | H10 / wave 1.4 — dd BDI `acc` term | **FAIL**: no `acc` field; `fhat − dd_loglik` = +0.2866 |

On a post-wave-1 build all three are expected to pass; if they do not, the study
still runs and reports them, because the arms that measure them (C9/D7, C1/C2)
are the point.

---

## 5. Methods — emphasis configurations

Common to every mcem job: `cond = NULL`, `rho = 1`, `link = "linear"`,
`num_threads = 1L`, `max_missing = 1e4`, `max_iter = 400L`, `maxN = 20 N`,
`max_time = timeout_s − 30` (so a clean `time_budget` stop is recorded before
the driver's hard kill), `tol`/`patience`/`xtol` at post-wave-1 defaults.

**Bounds are supplied by hand and contain the exact MLE by construction**, so
`auto_bounds` (H74, unfixed) enters only in the pipeline arm:

```
cr:  lower = c(1e-3, 0)                              upper = c(5 λ_gen, 5 λ_gen)
dd:  lower = c(1e-3, −5(λ0−μ0)/K, 0, 0)              upper = c(5 λ0, 0, 5 λ0, 0)
```

`mle_in_box` is recorded on every tree. The box is asymmetric and the init is
always explicit and off the `λ = μ` diagonal, so the H12/H21 midpoint start is
never used.

Inits: `far` = `c(2 λ_gen, 0.5 λ_gen)` for cr, `c(1.5 λ0, −0.5(λ0−μ0)/K, 0.25 λ0, 0)`
for dd. `mle` = the exact MLE (`μ̂ = 0` passed as 0).

### CR

| id | sampler | `N` | init | reps | trees |
|---|---|---|---|---|---|
| C1 | bdi | 200 | far | 3 | all |
| C2 | dynamic_fresh | 200 | far | 3 | all |
| C3 | bdi | 200 | mle | 1 | all |
| C4 | dynamic_fresh | 200 | mle | 1 | all |
| C5 | bdi | 1000 | far | 1 | subset S |
| C6 | dynamic_fresh | 1000 | far | 1 | subset S |
| C7 | bdi | 50 | far | 1 | subset S |
| C8 | dynamic_fresh | 50 | far | 1 | subset S |
| C9 | both, fixed θ | 2000 | — | 2 (bdi) / 3 (thin) | subset S |
| C10 | pipeline | — | — | 1 | subset S |
| C12a/b | bdi / thinning, box × 10 | 200 | far | 2 | `n = 50`, `ε = 0.6` |

Subset S is the first 5 trees per cell by index (4 at `n ≥ 100`).

C1/C2 × C5..C8 give the `N ∈ {50, 200, 1000}` scaling arm. C1 paired with C3 (and
C2 with C4) separates path and stopping error from the Monte Carlo noise around
the MLE. C12 measures the box dependence of the stopping rule directly rather
than assuming the wave-1.5/1.6 fix.

**C9, the fixed-θ importance-sampling arm.** At each of 9 θ points — the exact
MLE and `±1 SE`, `±2 SE` along each axis, `μ` clipped at 0 — the E-step alone is
run at `N = 2000`: BDI via `emphasis:::.augment_tree_bdi`, thinning via
`augment_trees` + `eval_logf` + `.is_fhat(n_zero_weight = rejected_zero_weights)`.
Recorded per point and replicate: `fhat`, `g(θ) = fhat − ℓ_exact(θ)`, ESS,
`sd(lw)`, non-finite count, `acc` when the build returns it. A θ-dependent slope
in `g` along the `μ` axis is the H7 signature; the intercept is the IS bias. A
single point at the MLE cannot see that slope, which is why the grid is not one
point.

**C10, the pipeline arm** runs all four stages with `num_threads = 1L`, and
records the `auto_bounds` box, the run log, each stage's parameters re-scored
with `ℓ_exact` (stage log-likelihoods are not comparable with each other —
H13/H27 — so they are stored, not interpreted), and whether the box contains
`θ_MLE0` and `θ_MLE1`.

### DD

| id | sampler | `N` | init | reps | trees |
|---|---|---|---|---|---|
| D1 | bdi | 200 | far | 2 | all |
| D2 | dynamic_fresh | 200 | far | 2 | all |
| D3 | bdi | 200 | mle | 1 | all |
| D4 | dynamic_fresh | 200 | mle | 1 | all |
| D5 | bdi | 1000 | far | 1 | subset |
| D6 | dynamic_fresh | 1000 | far | 1 | subset |
| D7 | both, fixed θ | 2000 | — | 2 / 3 | subset |

D7 evaluates `fhat` both raw and, when `acc` is returned, with `+log(acc)`, so
the report states whether wave 1.4 landed instead of assuming it.

---

## 6. Quantities recorded

One RDS per job, written by the child process, under
`results/<tier>/jobs/<job_id>.rds`. Schema:

- job identity and every factor level; `seed`; `timeout_s`; `est_cost_s`
- `outcome ∈ {ok, error, timeout, crash}` with the message — **no exclusions**;
  a failure is a result
- build fingerprint: `packageVersion`, library path, md5 of `DESCRIPTION` and of
  the shared object, R and DDD versions
- from the fit: `pars` (compact), `loglik`, `loglik_var`, `AIC`, `n_pars`,
  `stop_reason`, `iterations`, and the **full `details$mcem` trace**
  (`par1..par8`, `fhat`, `delta_max`, `rejected`, `num_trees`, `time`)
- `final_IS`: ESS, `n_rejected`, `rejected_zero_weights`, `acc`, count of
  non-finite `lw`, `sd(lw)` — each tolerated as absent on a pre-wave-1 build
- derived in the child: `ℓ_exact(θ̂)`, `Δℓ`, `z_j`, `at_bound`,
  `e_at_hat = loglik − ℓ_exact(θ_K)` **and** `e_at_prev = loglik − ℓ_exact(θ_{K−1})`
  from the stored trace, so the H20 conclusion holds whichever E-step the
  package reports; distance of the init from the MLE in SE units; for dd the
  mapping status and `K̂`
- C9/D7: per point and replicate `fhat`, `g`, ESS, `sd(lw)`, non-finite count,
  `acc`
- C10: everything the pipeline returns, plus `stage_ll_exact` and the two
  inside-box flags

---

## 7. Performance measures

Per cell × configuration (`cell = n × ε × λ_gen`, or the dd regime):

1. **Deficit** — median and 10th percentile of `Δℓ`. `Δℓ > 1e-3` is impossible
   if the reference is the maximum; any such tree is re-optimised from `θ̂_E` and
   logged as a reference failure.
2. **Decomposition** — for tree *i* with *R* replicates, `m_i = mean_r θ̂`,
   `s_i² = var_r θ̂`, systematic error `b_i = m_i − θ_MLE,i`. Reported in SE
   units: bias `B = mean_i b_i` with a tree-level bootstrap 95 % CI, the
   between-tree SD of `b_i`, and the pooled Monte Carlo SD `sqrt(mean_i s_i²)`.
   With 3 replicates each tree contributes 2 degrees of freedom, so only the
   pooled cell-level MC SD is used in a decision.
3. **Optimisation vs Monte Carlo** — the paired per-tree difference
   `Δℓ(C1 mean) − Δℓ(C3)`; and the regression of `log(−Δℓ)` on `log N` over
   {50, 200, 1000}. A slope near −1 means Monte Carlo dominates; a flat floor is
   optimisation error (`xtol`, the stopping rule, the M-step start — the arm
   cannot say which of the three).
4. **Reported log-likelihood** — median and 90th percentile of `|e|`. BDI under
   `cr` is expected at `|e| < 1e-6` once the final E-step is in place; thinning
   `e` is plotted against `ESS/N`. `var_r(fit$loglik) / mean_r(loglik_var)` is
   reported as a side check of H14 and enters no decision.
5. **Fixed-θ IS** — `g(θ)` with its replicate SE; a linear model of `g` on the
   SE-scaled coordinates tests θ-dependence; the intercept is the bias.
6. **Failures** — fraction with `stop_reason != "converged"`, timeouts, crashes,
   `NA` log-likelihood, iteration counts, seconds.
7. **A4 recovery** — `θ̂_E − θ_gen` against `θ_MLE0 − θ_gen`, with
   `θ_MLE2 − θ_gen` as the `n`-conditioned comparison.
8. **A5 pipeline** — final parameters against both `bd_ML(cond = 0)` and
   `bd_ML(cond = 1)` in SE units, stratified by the inside-box flags.

---

## 8. Decision rules

Fixed before the main tier runs.

**MCEM is adequate for `cr` at (`n`, `ε`, sampler, `N = 200`)** when all hold:

- (i) ≥ 95 % of C1/C2 fits stop with `"converged"` and none time out;
- (ii) median `Δℓ ≥ −0.1` nats and 10th percentile `≥ −0.5`;
- (iii) `|B| ≤ 0.1` SE for `λ`, and for `μ` on interior trees, with the
  bootstrap CI excluding `|B| ≥ 0.25`;
- (iv) pooled Monte Carlo SD `≤ 0.25` SE;
- (v) median `|e| ≤ 0.1` nats and 90th percentile `≤ 0.5`; for BDI under `cr`,
  `|e| ≤ 1e-6` throughout.

**Adequate with larger `N`** if only (iv) fails and the `N`-scaling slope is
`≤ −0.7`. **Inadequate** otherwise, naming the failing rule.

**The thinning IS estimator is unbiased for `cr`** if the C9 grid slopes are
within 2 replicate SE of 0 and `|intercept| ≤ 0.05` nats on ≥ 90 % of subset
trees.

**Adequate for `dd`** under the same rules with every threshold doubled
(`dd_loglik`'s own numerical tolerance is ~1e-3 nats and `dd_ML`'s optimiser adds
more), plus ESS ≥ 20 in ≥ 90 % of final E-steps, plus max augmented `N < lx/2` in
every E-step. Trees flagged `K_unidentified` are excluded from the dd verdict and
counted. The dd verdict is issued only if the tier-0 reference assertions pass.

**The pipeline is reported descriptively, with no adequacy verdict.** Its target
is identified as `cond = 1` if the paired distance to `bd_ML(cond = 1)` is
smaller than to `bd_ML(cond = 0)` on ≥ 80 % of trees with both inside the box.
No verdict is issued because the survival GAM (H91) and the unfixed
`auto_bounds` criterion (H74) are confounders this study cannot remove.

---

## 9. Figures and tables

One claim per figure.

| id | content | claim it settles |
|---|---|---|
| F1 | `Δℓ` boxplots by `n × ε`, faceted by sampler, at `N = 200`, reference lines at −0.05 and −0.5 | how far MCEM lands from the exact MLE |
| F2 | stacked variance components of `z_λ` (bias / between-tree / Monte Carlo) | whether what remains is noise or bias |
| F3 | `−Δℓ` vs `N`, log-log, slope per sampler | 1/N decay vs an optimiser floor |
| F4 | paired `Δℓ` at `init_far` vs `init_mle` | the cost of starting away from the MLE |
| F5 | `e` vs `ESS/N` (thinning) and the BDI `|e|` histogram | whether the reported log-likelihood is the log-likelihood |
| F6 | `g(θ)` across the 9-point grid with replicate error bars | IS bias and its θ-gradient (H7) |
| F7 | `θ̂ − θ_gen` vs `θ_MLE − θ_gen` | estimator error relative to MLE sampling error |
| F8 | pipeline estimates vs both conditioned MLEs, with the inside-box fraction | which target the pipeline reaches; H74 |
| F9 | dd deficits by regime and sampler | dd behaviour where the BDI proposal holds |
| F10 | outcome counts and iteration counts per cell | failure and cost structure |

Tables: **T1** per-cell summary; **T2** reference reconciliation (tier 0 results
and the sentinels); **T3** decision-rule verdicts; **T4** pipeline; **T5**
fixed-θ surface.

---

## 10. Tiers and compute plan

All timings were measured on the **pre-fix** build. The post-wave-1 stopping
rule changes iteration counts, and the counts already vary by more than a factor
of ten between runs at the same `n` (6 to 104 iterations observed). **The
wall-clock figures below are a forecast; the smoke tier's measured cost table is
what sizes the main tier.**

### Measured costs (1 thread, this Mac)

| operation | `n` = 20 | 50 | 100 | 200 |
|---|---|---|---|---|
| mcem bdi, `N` = 200 | 2.0 s / 16 it | 10.5 s / 26 it | 15.4 s / 10 it | 15.1 s / 6 it |
| mcem thinning, `N` = 200 | 0.8 s / 17 it | 1.8 s / 22 it | 25.7 s / 104 it | 25.8 s / 44 it |
| per iteration, bdi | 0.13 s | 0.40 s | 1.54 s | 2.51 s |
| per iteration, thinning | 0.04 s | 0.08 s | 0.25 s | 0.59 s |
| mcem bdi, `N` = 1000 | — | 69 s / 30 it | 98 s / 12 it | — |
| mcem thinning, `N` = 1000 | — | 9 s / 19 it | 147 s / 96 it | — |
| dd mcem, `N` = 200 (bdi / thin) | 15.6 / 1.7 s | 2.7 / 3.9 s | 22.0 / 19.9 s | — |
| fixed-θ E-step, `N` = 2000 (bdi / thin) | — | 1.0 / 0.8 s | 22.8 / 51 s | — |
| `emphasis_pipeline`, all four stages | 33 s | 66 s | 266 s | — |
| `bd_ML` | < 0.1 s | | | |

`dd_ML` (3 free parameters, `ddmodel = 1`, `cond = 0`, subplex,
`res = 10(n+1)`), measured today:

| `n` | 25 | 37 | 59 | 98 |
|---|---|---|---|---|
| `dd_ML` | 8.6 s | 10.1 s | 35.2 s | 74.5 s |
| one `dd_loglik` | 0.19 s | 0.09 s | 0.32 s | 0.74 s |

This is why the dd arm caps `n` at 80 and why the observed-information matrix for
dd uses 19 central differences rather than `numDeriv`'s Richardson scheme (70+
ODE solves). With the cap and the cheap Hessian, the whole dd reference costs
30–103 s per tree (measured in the smoke run). Without them, an `n = 98` tree
took over 7 CPU-minutes.

`dd_ML` is never run on a CR tree: with `K → ∞` it becomes a boundary search and
took 100–516 s in `scratch/refs_timing.out`.

### Tiers

**Tier 0 — references and build sentinels.** Serial. Hard gate: writes
`results/GATE.ok`; `03-fit.R` refuses `--tier main` without it. About 60 s on a
build where the `μ ≥ λ` BDI branch works; on the pre-fix build the eight `μ ≥ λ`
probes took several minutes between them (H12's error path is retried), so each
draw is time-limited to 120 s and a timeout is recorded as a sentinel failure
rather than a stall.

**Smoke, ~20 min wall end to end (tier 0 included).** 9 CR trees (`n ∈ {20, 50}` × 4 `ε`, plus one at
`n = 100`) and 3 DD trees, one per regime; C1–C4 and D1–D4 at one replicate;
one fixed-θ grid per tree at `N = 500`; one pipeline run per CR tree.
Concludes: the plumbing works, the schema is right, the references hold, and the
per-configuration cost table exists. **No statistical statement.** It produces
`results/smoke/cost-calibration.csv`, which fixes `trees_per_cell` and the
per-job timeouts (`3 ×` the measured cost, floor 120 s), and a list of
configurations that timed out and therefore have no measured cost.

Measured on the pre-fix build, 2026-09-13, with the machine shared (10 workers
plus one unrelated R process on 10 cores, so each child ran at about 60 % of a
core and these figures are roughly 2× the dedicated-machine cost):

```
00-selfcheck  ~5 min     mostly the eight mu >= lambda BDI probes (H12)
01-simulate   0.3 s      12 trees
02-reference  131 s      9 CR (bd_ML < 0.1 s each) + 3 DD (dd_ML 29-99 s each)
03-fit        918 s      81 jobs: 75 ok, 0 error, 6 timeout, 0 crash
04-analyse    0.5 s      48 fit rows, 351 fixed-theta rows, 9 pipeline rows,
                         10 figures, 5 tables
```

About 20 minutes end to end, and a second `./run.sh smoke` finished in 2.5
minutes because `03-fit.R` skipped all 81 finished jobs.

The 6 timeouts are information, not breakage: 3 thinning MCEM fits (`C2` at
`n = 50` and `n = 100`, which ran past 300 s and 600 s — the pre-fix build took
104 iterations at `n = 100`) and 3 thinning fixed-θ grids. The hard kill fired
in each case and a `timeout` row was written, which is what the mechanism is
for. `04-analyse.R` lists the configurations that timed out and therefore have
no measured cost, so the main tier keeps the default timeout for those.

**Main, target 2–3 h wall at 10 workers.** The full CR and DD designs of §3 and
§5. Estimated from the table above at roughly 60 000 core-seconds, so ~1.7 h at
perfect packing and ~2.2 h with tail imbalance and process-spawn overhead
(0.5 s × ~1700 jobs ≈ 15 min of core time).

**Extension, ~30–40 min, only if the main tier finishes.** `n = 200`,
`ε ∈ {0.3, 0.9}`, 5 trees each: thinning ×3 replicates and one BDI fit per tree
under a cap, plus a BDI-only fixed-θ check on 2 trees per cell (thinning at
`n = 200` would exceed the 120 s internal C++ E-step cap). No pipeline arm, no
`init_mle` arm, no `N`-scaling arm: the extension reports BDI's cost scaling and
whether it reaches the MLE within the cap. Dry-run job table (measured cost
model): 44 jobs, ~12 000 core-seconds.

Per-job timeouts are `max(the kind/size default, 3 × the cost model, 120 s)`.
Once `results/smoke/cost-calibration.csv` exists, the cost model is that table,
which is the design's "3 × the measured cost" rule made operational; the
pre-fix formula is the fallback, and it under-predicted the smoke run by about
3×.

### What the dry run already shows (pre-fix build — not results)

These are from the smoke run against the pre-fix scratch build. They are
recorded because they show each arm is wired to the signal it was built to
measure, not because they are findings.

- **H20 lag.** On a `n = 100` C1 fit: `fit$loglik − ℓ_exact(θ_K) = 6.1e-03`
  while `fit$loglik − ℓ_exact(θ_{K−1}) = 6.9e-12`. The reported value is the
  E-step at the previous iterate, exactly. Recording `e` at both iterates makes
  A2's conclusion independent of whether wave 1.5/1.6 landed.
- **C9 arm, cr.** BDI `g(θ) = 0.000` with ESS = `N` at every grid point of every
  CR tree — the zero-variance weights. Thinning `g` between −0.008 and +0.019 at
  `n = 20` and `n = 50` with ESS 390–450, and −0.648 at `n = 50`, `ε = 0.6` with
  ESS 87.
- **C9/D7 arm, dd.** BDI `g = +0.303` in regime A (ESS 175) and `−0.680` in
  regime B (ESS 38); `−58.4` in the regime-C stress cell. The first is the H10
  acceptance term with the wrong sign convention on this build; the last is the
  H11 zero-rate channel. The stress cell is reported separately for this reason.
- **`init_mle` arm.** C3/C4 stop after 3 iterations with `Δℓ = 0.000` and
  `|e| = 0` — the patience counter fires immediately when the iterate does not
  move. That is the reference point the C1/C3 pairing needs.
- **`μ̂ = 0` frequency.** 6 of the 9 smoke CR trees have an exact MLE on the
  `μ = 0` boundary. `z_μ` is undefined for those, which is the reason `Δℓ` and
  not a z-score is the primary scalar.

### Drop order

Applied by the driver when observed throughput after the first ~100 jobs
projects past budget:

1. `N = 1000` at `n = 100`
2. the pipeline arm at `n = 100`
3. the third replicate of C1/C2 at `n = 100`
4. `n = 100` trees per cell, 12 → 8
5. the `N = 50` level

Never dropped: tier 0, the `init_mle` arms, the fixed-θ grid, the CR `N = 200`
mcem arm, and the floor of 10 trees per cell.

### Execution

`03-fit.R` runs **one OS process per fit** (`callr::r_bg`), so a hang or a C++
abort costs one job and not the run. Jobs are ordered longest-first and
dispatched dynamically to 10 slots. Each job carries a hard timeout; emphasis's
own `max_time` is set 30 s below it so a clean `time_budget` stop is recorded
before the kill. Each job writes its own RDS, so a rerun skips finished jobs and
a crash loses nothing. `OMP_NUM_THREADS`, `RCPP_PARALLEL_NUM_THREADS` and
`TBB_NUM_THREADS` are pinned to 1 in every child (H62/H63: the M-step would
otherwise take all cores). Processes are **spawned, not forked**, so no two
thinning jobs share the clock-seeded C++ engine state (H47).

`04-analyse.R` reads job RDS files only and never runs a fit; it refuses to pool
rows from different build fingerprints.

---

## 11. Known distortions and how each is handled

Wave 1 of `dev/audit/02-findings.md` is being implemented now and changes cr/dd
numbers; the design targets post-wave-1 behaviour. Each item below is either
neutralised by the design or measured by an arm.

| finding | status | handling |
|---|---|---|
| H74 `auto_bounds` may exclude the MLE | wave 3, unfixed | mcem arms never use `auto_bounds`; the hand box contains every exact MLE and `mle_in_box` is recorded; the pipeline arm records whether the auto box contains `θ_MLE0` and `θ_MLE1`, so H74 is measured |
| H21 box-scaled stopping rule | wave 1.5/1.6 | the box is on the parameter scale; iterations and the full `delta_max` trace are recorded; the `init_mle` arm exposes stopping displacement directly; the C12 box × 10 sub-arm measures box dependence without relying on the fix |
| H12/H54 BDI fails at `μ ≥ λ`; symmetric-box midpoint init | wave 1.3/1.9 | `init_pars` is always explicit and off the diagonal; tier-0 sentinel (c2) reports the state of the fix; `ε = 0.9` trees whose MLE has `μ̂ ≥ λ̂` are kept and flagged (`mle_r_negative`), so a surviving H12 shows as `e_step_failure` in that stratum |
| H20 reported loglik lags one iteration | wave 1.5/1.6 | `e` is computed at `θ_K` **and** at `θ_{K−1}` from the stored trace; the primary criterion scores the returned `pars` against `ℓ_exact`, which is lag-independent |
| H7 thinning envelope | wave 1.8 | C9's grid slope in the `μ` direction is the signature; sentinel (f) pins it; the deficit and `z` of C2 measure the downstream effect |
| H10/H59 dd BDI omits `−log P(accept)` | wave 1.4 | D7 reports `fhat` raw and `+log(acc)`; sentinel (h) says whether the field exists; estimates are unaffected either way |
| H11 dd/linear `±Inf` weights, frozen M-step | wave 1.1–1.5 | non-finite counts and `stop_reason` recorded; regime A's `K = 40` puts the zero of `λ(N)` inside the range the augmentation reaches; an iteration that returns the init unchanged is visible in the trace |
| H1/H16 underflow cut at `lw = −745` | wave 1.7 | `|ℓ| < 400` for every tree in the design at `λ_gen = 1`; the `λ_gen = 0.1` time-scale cell shifts `ℓ` by about `(n−2)·log 10 ≈ 110` and is the probe; failures are counted regardless |
| H22/H30 `maxN` alternation | wave 1.6 | `maxN = 20 N` explicitly on every call |
| H23 undocumented 120 s E-step cap | wave 2.6 | thinning E-steps at `n ≤ 200`, `N ≤ 1000` are ≤ 1.5 s; the `N = 2000` fixed-θ arm is restricted to `n ≤ 100` (51 s measured); a job stopped by the cap surfaces as `e_step_failure` |
| H62/H63 thread race inflating `fhat` | wave 1.7/2.6 | `num_threads = 1L` in every control list including the pipeline's (whose default is `detectCores() − 1`), and the env vars pinned in the child |
| H47 C++ RNG not seedable | wave 3.4 | BDI fits seeded and reproducible; thinning replicates are the only handle and no thinning result is reported from a single fit; spawned (not forked) processes |
| H4 conditioned M-step penalty diluted | wave 1.2 | every scored arm uses `cond = NULL`; the pipeline arm compares with both `cond = 0` and `cond = 1`, so an unfixed H4 shows as proximity to `cond = 0` |
| H91 survival GAM error | not fixed (author question) | stated as the pipeline arm's floor; no adequacy verdict for that arm |
| H14 `loglik_var` is another estimator's variance | not fixed | stored, compared with the replicate variance, used in no decision |
| H13/H27 stage logliks not comparable | wave 2.8 | only final parameters are scored; stage parameters are re-scored with `ℓ_exact`; stage logliks stored, not interpreted |
| H35 AIC counts bound-fixed parameters | wave 1.9 | AIC is not used |
| H2/H29/H79 `rho` | wave 2.1 | `rho = 1` throughout; the study makes no claim about `rho < 1` |
| H58 Picard non-convergence under strong DD | wave 2.9, unfixed | confined to regime C, labelled a stress cell, reported separately |
| `ref_check.R` `pars2` slot error | not a package defect | closed by tier-0 assertions (a1)/(a2); `NOTES.md` and `ref_check.R` are left for the owner to correct |

---

## 12. What the study cannot conclude

- Nothing about `rho < 1`, the exponential or gaussian links, `d`/`nd`
  (D-covariate) models, the general `dd` model with `γ_N` free, or M-active
  models. Those have no exact reference; H6, H8 and H9 say why.
- Nothing about conditioned fits except through the pipeline arm, where H91 and
  H74 make the comparison descriptive.
- `"adequate for cr"` is bounded at `n ≤ 100`. BDI's per-iteration cost grew from
  0.40 s at `n = 50` to 1.54 s at `n = 100` and 2.51 s at `n = 200`, so size
  scaling for BDI rests on three points plus a capped extension arm.
- The `N`-scaling arm attributes a floor to "the optimiser" collectively —
  `xtol`, the stopping rule and the M-step start — and cannot say which.
- Thinning cannot be seeded, so a thinning result is reproducible only in law;
  its Monte Carlo component includes any run-to-run non-determinism of the C++
  engine. Three replicates give two degrees of freedom per tree, so only the
  pooled cell-level MC SD carries a decision.
- MLEs at `μ̂ = 0` (common at `ε = 0` and at `n = 20`) leave `z_μ` undefined;
  those cells are compared on `Δℓ` and `λ` alone.
- The dd references carry ~1e-3-nat numerical tolerance, so dd deficits below
  that are unresolvable; the dd arm is `ddmodel = 1` with `γ_N` fixed,
  `λ0 > μ0`, `n ≤ 80`, and `K` identifiable by construction. Regime C may fail
  for reasons (H58) the study cannot separate from estimator error.
- `sim.bd.taxa` conditions on `n`, so A4 describes the `n`-conditioned design,
  not fixed-age sampling. `bd_ML(cond = 2)` is reported for that reason.
- The study certifies proximity to the exact MLE **under these control settings
  on these trees**. It does not establish that the shipped defaults are the right
  defaults at other tree sizes or time units.
- Every timing here is from the pre-fix build; wall-clock is a forecast, and the
  drop order is the only guarantee that the main tier finishes.

---

## 13. Files

```
dev/validation/
  00-design.md          this file
  R/00-common.R         exact references, design tables, job helpers, thresholds
  R/00-selfcheck.R      tier 0; writes results/GATE.ok
  R/01-simulate.R       trees        -> data/trees-<tier>.rds
  R/02-reference.R      exact MLEs, SEs, fixed-theta grids -> data/reference-<tier>.rds
  R/03-worker.R         the function each child process runs (one fit per process)
  R/03-fit.R            job table, callr pool, hard timeouts, resumable
  R/04-analyse.R        summaries, figures, tables, report  (reads RDS only)
  run.sh                tiers in order
  results/<tier>/jobs/  one RDS per fit
  figures/  tables/     <tier>-F*.pdf, <tier>-T*.md
```

Run:

```sh
EMPHASIS_LIB=/path/to/post-wave-1/rlib ./run.sh smoke
EMPHASIS_LIB=/path/to/post-wave-1/rlib ./run.sh main
```

`EMPHASIS_LIB` unset uses the normal library. `03-fit.R --dry-run` prints the
job table and its cost estimate without executing anything. Re-running a tier
continues where it stopped.

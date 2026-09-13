# emphasis (v0.4, `/Users/pancho/Code/emphasis`) — the estimator as implemented, and the audit hypotheses that follow

Synthesised 2026-09-13 from seven subsystem maps (inference driver, IS core C++, augmentation C++, BDI sampler, CEM/GAM/pipeline, simulation API, tests/CI/docs). Every claim carries a `file:line` anchor from those maps. Where two maps disagree, the disagreement is stated in place and collected in §1.11.

---

## 1. The estimator, end to end (as implemented)

### 1.1 Entry, resolution, inputs

`estimate_rates(tree, model, link, method, init_pars, cond, control)` (`R/inference.R:672-840`) is the single public inference entry point; `emphasis_pipeline()` (`R/pipeline.R:57-330`) sequences it.

- **Tree → branching times.** `.extract_brts` (`R/inference.R:216-232`): numeric → `sort(decreasing=TRUE)`; `phylo` → `ape::branching.times` sorted decreasing; a `simulate_tree()` list → `$tes`, else `prune_to_extant($tas)`. Result: ages, **crown age first**, length `n_tips - 1`.
- **Model → `model_bin`.** `.resolve_model` (`R/simulate.R:254-275`): `cr=c(0,0,0)`, `dd=c(1,0,0)`, `d=rd=ep=c(0,0,1)`, `nd=c(1,0,1)`; formulas map `N→slot 1`, `D|EP|E→slot 3`, anything else errors (`simulate.R:278-293`); any length-3 0/1 integer vector is accepted, so `c(1,1,0)` (M active) is reachable (`simulate.R:268-274`). Slots are `c(use_N, use_M, use_D)`.
- **Link → integer.** `.resolve_link` (`R/simulate.R:247-251`): linear=0, exponential=1, gaussian=2; a numeric input is `as.integer`-ed with no range check.
- **Control.** `estimate_rates_control(method)` (`R/inference.R:98-144`) supplies per-method defaults; `modifyList` with the user list; `.resolve_control_aliases` (`inference.R:180-213`) syncs `num_points↔num_particles`, `sample_size↔num_trees` (user-supplied name wins; back-ends read the old names). `lower_bound`/`upper_bound` are the only required entries, in the **compact** layout of length `2 + 2*sum(model_bin)` (`inference.R:689-702`). `rho` is taken from control with no validation (`inference.R:107`); C++ resets `rho ≤ 0 or > 1` to 1 (`inst/include/model.hpp:90`).

### 1.2 Parameter layouts

- **Compact (user-facing):** `c(beta_0, beta_active..., gamma_0, gamma_active...)`, active in slot order N, M, D; names from `.par_names` (`inference.R:235-247`).
- **8-element (C++ / internal):** `c(beta_0, beta_N, beta_M, beta_D, gamma_0, gamma_N, gamma_M, gamma_D)` (`R/simulate.R:296`; `model.hpp:58-63`). `.expand_pars` (`simulate.R:298-310`): active covariate `k∈{1,2,3}` → beta slot `k+1`, gamma slot `k+5`; inactive slots 0. `.contract_pars` (`inference.R:347-352`) inverts it. Bounds expand the same way with inactive slots pinned to `[0,0]` (`inference.R:752-753`).
- C++ reads **only** `model_bin[2]` (use_D) to choose an algorithm branch (`model.hpp:245-247, 256-257, 272-274, 360-361, 420-430`; `general_tree.hpp:185, 227`); the N and M terms always enter the linear predictor and are switched off only by R's zero-padding.

### 1.3 Rates and the complete-data log-likelihood `f`

Covariates on a node (`model.hpp:143-190`): `N = node.n` (lineages alive on the segment ending at the node); `M = pd/n`; `E = brts - tip_start` (extinction node), `brts - focal_tip_start` (augmented speciation node with `parent_id ≥ 0`), else `M` (so `D = 0` for observed nodes); `D = E - M`.

Links (`model.hpp:117-212`; identical in `general_tree.hpp:104-140`): linear `max(0, eta)`; exponential `exp(eta)`; gaussian `beta_0 * exp(-(beta_N N + beta_M M [+ beta_D D] - 1)^2/2)` (intercept multiplies, is not inside the exponent). `mu` likewise with `gamma`.

`log f = Model::loglik` (`model.hpp:359-471`), piecewise-constant path:

```
log f = Σ_{extinction nodes} log max(mu_i, 1e-300)
      + Σ_{speciation nodes, i ≠ closing node} log lambda_i          (via log_sum, model_helpers.hpp:126-153)
      - Σ_i Δ_i · n_i · (lambda_i + mu_i)                            (rates frozen at the END node of each segment)
      + n_obs·log rho + n_unsamp·log(1-rho)                          (model.hpp:459-467; no binomial coefficient)
```

For D + exponential the integral is the exact per-lineage sum `S_b·exp_integral(A_lam, beta_D, s_{i-1}, s_i)` (`model.hpp:381-392, 333-336`); for D + gaussian per-lineage erf integrals (`model.hpp:393-417, 348-357`); for D + linear the integral uses `n_i·(lambda_focal + mu_focal)` (`model.hpp:424`). "Speciation node" = tip, missing or unsampled node (`model.hpp:433-438`); `n_obs = 1 + #tip nodes = n_tips`. `f` is a per-lineage (labelled-history) density: `log lambda` per event, no `n·lambda`, no combinatorial factors.

### 1.4 Augmentation proposals

**(a) Thinning proposal (C++, `src/augment_tree.cpp:128-191`, `model.hpp:217-252`).** Used by `augment_trees`/`em_cpp`, i.e. by CEM, GAM, `.mcem_dynamic_fresh`, and by `simulate_tree(method="thinning")`. `create_tree` (`src/E_step.cpp:19-38`) maps `b_0 > … > b_{k-1}` to **forward** times `s_i = b_0 - b_{i+1}` (nodes with `n = 2+i`, `t_ext = 1e11`) plus a closing sentinel at `T = b_0` (`n = n_tips`); the crown (t=0, 2 lineages) is not a node. On each interval `[cbt, next_bt)`: envelope `Λ* = max(nh(cbt), nh(next_bt))` (`augment_tree.cpp:138-140`) with `nh(t) = n·lambda·(1 - rho·exp(-max(mu,1e-10)(T-t)))` evaluated on the node at/after `t` with `pd` recomputed at `t` (`model.hpp:238-252`); throw `augmentation_lambda` if `Λ* > max_lambda` (`augment_tree.cpp:141`); candidate `t* = cbt - log(u)/Λ*`, accepted with probability `nh(t*)/Λ*` (not clipped at 1, lines 148-151). Lifetime from `extinction_time` (`model.hpp:217-235`): with probability `(1-rho)e^{-mu r}/(1 - rho e^{-mu r})` an unsampled extant (`t_ext = 5e10`), else `t* + TruncExp(mu, T - t*)` by rejection (`model_helpers.hpp:157-165`), `mu = extinction_rate` (no D term) on a node whose `pd` is still 0. Parent uniform over alive non-extinction **nodes** with `brts < t* < t_ext` (`augment_tree.cpp:154-166`; crown lineages never eligible; `-1` if none). `num_missing > max_missing` throws `augmentation_overrun`. After the loop `compute_pendant_pd` (`augment_tree.cpp:45-77`) sets `tip_start = 0` for `parent_id == -1` nodes else `brts`, fills `pd = P(brts)` (`model_helpers.hpp:216-229`, `P(t) = Σ_{alive, non-extinction, brts ≤ t} (t - tip_start)`), and `focal_tip_start` from an id→pendant-start map.

`log q = Model::sampling_prob` (`model.hpp:255-328`), with `tips` starting at 2 and incremented per observed node, `Ne` = augmented lineages alive:

```
log q = Σ_{missing i}   [ log(n_i mu_i lambda_i) - mu_i (t_ext,i - s_i) - log(2·tips_i + Ne_i) ]      (line 304)
      + Σ_{unsampled i} [ log(n_i lambda_i (1-rho)) - mu_i (T - s_i)   - log(2·tips_i + Ne_i) ]      (line 308)
      - Σ_i n_i lambda_i [ Δ_i - (rho/mu_i)(e^{-mu_i(T-s_i)} - e^{-mu_i(T-s_{i-1})}) ]              (lines 295-296; mu floored 1e-10)
```

For D + exponential the compensator uses `S_b·exp_integral(...)/Δ_i` instead (`model.hpp:278-293`).

**(b) BDI sampler (pure R, `R/bdi.R`).** Default for MCEM (`inference.R:111`) and for `simulate_tree(tree=…)` (`simulate.R:143`), gated by `.bdi_supported` (`bdi.R:100-105`): `model_bin[2]==0 && model_bin[3]==0 && link ∈ {0,1}`; **rho is not examined**. Simulates the birth–death process conditioned on the reconstructed tree with Nee et al. doomed-lineage rates: per missing lineage birth `λu(t)`, death `μ/u(t)`, immigration from the `k(t)` observed lineages `2kλu(t)`, `u = 1 - p` with `p(t) = (λ0-μ0)/(λ0 - μ0 e^{-(λ0-μ0)(tp-t)})` (`bdi.R:16-19, 451-456`). Under CR the cumulative hazard is analytic (`bdi.R:36-66`) and the event time is a time change solved by `uniroot` (`bdi.R:76-89`); under DD (N-only) a mean-field backward/forward ODE fixed point (`bdi.R:307-392`) feeds a frozen-rate Gillespie (`bdi.R:458-483`) with rejection of draws leaving missing lineages alive at `tp` (`bdi.R:514`). `log g` accumulates `log(la)` per birth and per immigration, `log(mu_r)` per death, minus the integrated total rate (`bdi.R:440-506`). Draws are converted to the tree data frame (`bdi.R:527-611`; `pd ≡ 0`, no `focal_tip_start`, parent = id of the most recent observed branching, `bdi.R:556`) and scored with C++ `eval_logf` (`bdi.R:687-691`); `eval_logf`'s `logg` is discarded.

`.run_mcem` (`inference.R:355-415`) selects `.mcem_bdi` when `ctrl$sampling == "bdi"` and the gate passes, else `.mcem_dynamic_fresh` (thinning), with a message only when `verbose` (`inference.R:359-364`).

### 1.5 Weights and `fhat`

- **Log-weight:** `lw_i = logf(theta, z_i) - logg(theta, z_i)`, both at the theta that generated `z_i` (`src/E_step.cpp:90-92`; `bdi.R:696`).
- **Acceptance (thinning):** keep iff `isfinite(lw) && exp(lw) > 0` (`E_step.cpp:94`); else `++rejected_zero_weights`. Overrun/lambda/other exceptions are counted separately. Stop at `N = sample_size` accepted, or `maxN` attempts, or **120 s** (`E_step.cpp:73`; no caller passes `max_time_seconds`, `src/mcem.cpp:27,46`, `src/rcpp_mce.cpp:61`). Fewer than N → throw (`E_step.cpp:128-130`) → R sees `NULL`.
- **`fhat` (thinning, C++):** `m = max lw`; `sum_w = Σ exp(lw_i - m)` in 100-digit floats (`inst/include/precision_weights.hpp:10-24`); `fhat = log(sum_w / S) + m`, `S = N + rejected_zero_weights` (`E_step.cpp:143-144`). R's `.is_fhat` (`R/de.R:60-76`) reproduces this with `n_zero_weight`. Optional `bias_correct` branch: `mean(lw) + log(1 + Σ_{k=2}^{K} m_k/k!)`, K=2 default, only when `n_zero_weight == 0` (`de.R:82-100`).
- **`fhat` (BDI):** log-mean-exp over accepted trees only, denominator `n_valid` (`bdi.R:696-699`); no finite-weight filter.
- **`fhat` (CEM shared-tree mode):** for particle `j`, `log((1/|Z|) Σ_{z∈Z} exp(logf(z,θ_j) - logg(z,θ_j)))` with both terms at `θ_j` regardless of which particle drew `z` (`de.R:317-357`).
- **Weights handed to the M-step:** thinning path overwrites `E.weights[i]` with `exp(lw_i - m) ∈ (0,1]`, **unnormalised**, `Σw ∈ [1, N]` (`E_step.cpp:133-139`); BDI path normalises to **mean 1**, `Σw = N` (`bdi.R:796-799`).
- **ESS** (R only): `(Σw)^2/Σw^2` on `exp(lw - max)` over finite `lw` (`de.R:111-117`).

### 1.6 M-step

`src/M_step.cpp:41-94`. Objective `Q(θ) = Σ_i w_i · loglik(θ, z_i)` (TBB `parallel_reduce`); NLopt SBPLX minimises `-Q(θ)` or, when a conditional is supplied, `-Q(θ) + log P_tree(θ)` (line 61; the comment at lines 58-60 states the normalised form). Bounds `lb8/ub8`, start = current theta, stop on `xtol_rel = ctrl$xtol` (1e-3) only (`M_step.cpp:81-89`). Any negative `nlopt_result` throws (`src/rsbplx.cpp:144-146`). `m_cpp` (`src/rcpp_mcm.cpp:52-100`, BDI path) rebuilds nodes from `brts/n/t_ext` only (`pd = tip_start = 0`, `id = parent_id = -1`). `cond_fun(pars8) = log max(predict_survival(gam, compact(pars8)), 1e-300)` (`inference.R:333-344`).

### 1.7 MCEM iteration and stopping

**Thinning driver `.mcem_dynamic_fresh`** (`R/emphasis.R:1-171`), for `k = 1..max_iter` (200):
1. `em_cpp(theta_{k-1})` → `fhat_k = fhat(theta_{k-1})` (E-step at the *old* theta), `theta_k` = M-step estimate, per-tree `logf/logg` (`emphasis.R:44-61`). `max_lambda` is hard-coded 1e6 (`emphasis.R:50`).
2. On `NULL` (any error): `fail_streak++`, `maxN ← min(2·maxN, 50000)`, `theta ← 0.8·theta + 0.2·(lb+ub)/2` clamped (`emphasis.R:64-85`); stop `e_step_failure` at 8 consecutive failures (then `.mcem_warn_estep`, `emphasis.R:174-217`); `prev_pars`/streak are not updated; on the next success `maxN` is reset to the original (`emphasis.R:89`).
3. `delta_k = max_j |theta_k[j] - theta_{k-1}[j]| / (ub_j - lb_j)` (range 1 where `ub == lb`) (`emphasis.R:28-29, 93`).
4. Row `(theta_k, fhat(theta_{k-1}), delta, rejected, num_trees, time)` (`emphasis.R:97-105`).
5. Converged when `delta_k < tol` (1e-3) for `patience` (3) consecutive iterations (`emphasis.R:114-122`); `time_budget` when elapsed > `max_time` (3600 s). Sample size `N = 200` is constant; no growth, no averaging.

Reported: `pars = theta_K`; `loglik_var` = variance over 200 bootstrap resamples of the **K=2 bias-corrected** estimator on the last E-step's `lw` (`emphasis.R:136-142`; `de.R:476-489`); `final_IS` with `.is_fhat`, ESS, rejection counts (`emphasis.R:145-161`).

**BDI driver `.mcem_bdi`** (`bdi.R:722-915`): same loop shape with `.augment_tree_bdi` as E-step, `m_cpp` as M-step, same convergence rule, on failure perturbs only (no `maxN` doubling, `bdi.R:778-788`), `rejected` always 0 (`bdi.R:804-807, 851`).

**Init** (`inference.R:714-738`): `init_pars = (lb+ub)/2` compact; for the linear link with covariates all slopes start at `clamp(0)`, then `.validate_linear_init_pars` (`inference.R:256-324`) guards `lambda > 0` / `mu ≥ 0` at the positive covariate extreme only. Ignored for `cem` and `gam`.

**Reported loglik:** `.run_mcem` takes the last finite `fhat` in the table = `fhat(theta_{K-1})` (`inference.R:410-412`); `estimate_rates` subtracts `cond_fun(theta_K)` when `cond` is supplied (`inference.R:819-822`); `AIC = -2·loglik + 2·n_pars` with `n_pars = length(compact)` including bound-fixed parameters (`inference.R:816, 824`).

### 1.8 CEM initialisation (`emphasis_cem`, `R/de.R:567-902`)

Particles are 8-element vectors in `[lb8, ub8]` (`de.R:172-192`); defaults `num_points 50`, `sample_size 1`, `maxN 10`, `max_iter 50`, `disc_prop 0.5`, `tol 1e-4`, `patience 5`, `sd_decay 0.85`, `sd_min_frac 0.01` (`inference.R:121-133`); pipeline overrides `max_iter 20`, `num_particles 50`, `num_trees 5` (`pipeline.R:198`). Per iteration: evaluate all particles (`fhat` via `augment_trees` + `.is_fhat`; NULL → NA, `de.R:213-288`), subtract `cond_fun` (`de.R:651-655`), rescue if all NA (`maxN·10`, `max_missing·10`, abort when `maxN > 10000`, `de.R:659-674`), `max_missing·1.1` on any overrun (`de.R:677`), best tracking, plateau rule `improvement < tol`, stop when `plateau ≥ patience` **and** all `sd_vec ≤ sd_floor` (`de.R:705-708`) or `max_time`, resample elites `ceil(disc_prop·n_valid)` + Gaussian perturbation truncated to bounds (`de.R:407-450`), SD decay on non-improvement (`de.R:733-737`). Final: re-evaluate elites with `max(20, sample_size)` trees, `obtained_estim = Σ_e softmax(fhat_e)·theta_e` (`de.R:793-848`), one more evaluation at that point → `best_IS$fhat - cond_fun` (reported loglik, `inference.R:457-463`), bootstrap `loglik_var` (`de.R:881-885`).

### 1.9 GAM stage and `auto_bounds`

`auto_bounds` (`R/gam.R:133-396`): centre from `r_hat = log(N/2)/T`, `lam_hat = 1.2 r_hat`, `mu_hat = 0.2 r_hat`; log-scale intercepts for any `link_int != 0` (`gam.R:170-179`); feasibility = ≥ half of `n_test` forward simulations end `done` with tips in `[max(2, 0.1N), 10N]` (`gam.R:401-416`); three bisection phases; box = min/max of feasible points ± `margin·span` clamped to `.wide_bounds` (`gam.R:344-352, 661-691`); survival GAM `survived ~ Σ s(p_j)`, binomial, on 500 LHS points (`gam.R:371-388, 34-68`). If no feasible centre: returns the wide box with `survival_gam = NULL` (`gam.R:207-214`).

`.run_gam` (`inference.R:469-571`): LHS grid (`n_grid`, pipeline 150) in compact space; per row `augment_trees` (`sample_size` 200) + `eval_logf` + `.is_fhat` (`gam.R:739-847`); `fhat ← fhat - log P_surv` (no floor, `inference.R:527-530`); `train_likelihood_GAM` additive Gaussian GAM (`gam.R:867-910`); `find_MLE` = L-BFGS-B on `-predict(gam)` from the best grid row (`gam.R:936-971`); **loglik = predicted smooth value at the optimum**; `loglik_var = NA`.

### 1.10 Pipeline choice

`emphasis_pipeline` (`pipeline.R:57-330`): stages `bounds → gam → cem → mcem`, each in `tryCatch`. `control$rho` reaches `auto_bounds` only (`pipeline.R:125`). MCEM init = `cem` pars if finite else `gam` (`pipeline.R:231-239`). Final result = first of `(mcem, cem, gam)` with finite loglik (`pipeline.R:289-302`); no loglik comparison.

### 1.11 Where the maps disagree or stand in tension

1. **Thinning `q` vs measured agreement.** Maps 2 and 3 derive that `sampling_prob`'s `-log(2·tips + Ne)` term (`model.hpp:304, 308`) equals neither `log n_i` nor `log |alive_ids|` and that for M/D models three different `mu`s and two different `lambda`s enter sampling vs scoring, so `q` is not the sampler's density. Map 4 measured on one CR tree that thinning `fhat = -16.89` (ESS 1312/4000) against BDI's `-16.87`, which it verified equals the Nee et al. crown likelihood by `integrate()`. The two are reconcilable only if the extra factor integrates to a θ-independent constant near 1 for CR; that is untested (→ H8, H9).
2. **`log_sum` at `lambda = 0`.** Map 2 states `lambda = 0` gives `-inf` (`model_helpers.hpp:132-135`); Map 4 reports a "sign bug in `log_sum::result`" producing `+Inf` for some trees (`model_helpers.hpp:129-136`) with 56/300 BDI-DD trees at `±Inf`. Direct check needed (→ H31).
3. **Recovery on E-step failure.** Map 1 states the BDI path only perturbs (`bdi.R:778-788`) while the thinning path doubles `maxN` and then resets it (`emphasis.R:68, 89`); README describes one behaviour for "MCEM".
4. **`maxN` defaults.** Three code sites, three values: `.run_mcem` `max(2000, 10·sample_size)` only when `NULL` (`inference.R:365`), `.augment_tree_internal`/GAM `max(2000, 200·n)` (`simulate.R:475`, `inference.R:518`), CEM 10 (`inference.R:126`), BDI `max_tries = sample_size` or `5·sample_size` (`bdi.R:663`). Not a contradiction between maps, but the maps' summaries of "the default maxN" differ because they describe different sites.
5. **Installed vs source.** Map 6 found the installed library 0.4 lacks the `.bdi_supported` fallback at `simulate.R:345`; all maps describe the source.
6. **What `E.weights` holds.** `bdi.R:792-795` (comment) says the thinning M-step receives log-weights; `E_step.cpp:132-139` converts them to linear max-scaled weights first (Maps 2, 4 agree; the disagreement is code-vs-comment).

---

## 2. Conventions table

| Item | Convention as implemented | Anchor |
|---|---|---|
| Time direction, R side | Ages: crown age first, decreasing to 0 at present | `inference.R:216-232` |
| Time direction, C++ augmentation | Forward from crown (0) to present `T = brts[1]`; `create_tree` maps `b_{i+1} → T - b_{i+1}`; closing sentinel node at `T` | `E_step.cpp:19-38` |
| Time direction, forward simulator | Forward `t ∈ [0, max_t]`; L-table output converted to ages (`max_t - date`, `-1` alive) | `general_tree.hpp:152-295`, `div_tree.cpp:45-51` |
| L-table | DDD convention: col1 birth age, col2 parent label, col3 label (one crown clade negative), col4 death age or -1 | `div_tree.cpp:45-51`, `simulate.R:414-453` |
| `brts` ordering | Strictly decreasing, crown first, length `n_tips - 1`; never checked in C++ | `E_step.cpp:19-38` |
| Node `n` | Lineage count on the half-open segment **ending** at the node; first node `n = 2`; closing node `n = n_tips` | `model_helpers.hpp:58` |
| Rates on a segment | Evaluated at the segment's END node (`Δ_i · n_i · (λ_i + μ_i)`) | `model.hpp:295-296, 424` |
| Parameter slots (C++) | 0-based `pars[0..7] = (β0, βN, βM, βD, γ0, γN, γM, γD)`; only `model_bin[2]` read | `model.hpp:58-63` |
| Parameter slots (R compact) | 1-based `(β0, β_active…, γ0, γ_active…)`; `n_lam = 1 + #active` | `simulate.R:298-310` |
| Node ids | 0-based; id `k` ↔ `brts[k+2]` (R 1-based); crown is not a node; closing node id `n-2` (thinning) / `-1` (BDI); augmented species: speciation + extinction rows share one id | `augment_tree.cpp:196-212`, `bdi.R:538-560`, `simulate.R:428` |
| `t_ext` sentinels | `0.0` extinction node, `1e11` observed tip, `5e10` unsampled extant, anything else = augmented speciation node with its own extinction time; duplicated as literals in R | `model_helpers.hpp:48-50, 105-108`; `simulate.R:417-418`; `bdi.R:544, 572` |
| `tip_start` | Observed nodes and augmented nodes with `parent_id == -1`: `0` (pendant age counted since the crown); other augmented nodes: birth time; observed lineages never reset by their own later speciations. Forward simulator: true pendant ages (parent reset at each speciation) | `augment_tree.cpp:50`, `model_helpers.hpp:216-229`, `bdi.R:287-291`, `general_tree.hpp:264-267` |
| `pd`, `M`, `D` | `pd = P(brts)`; `M = pd/n`; `D = E - M` with `E = M` (so `D = 0`) for observed nodes; BDI trees carry `pd ≡ 0` | `model.hpp:177-185`, `bdi.R:605` |
| Log base / sign | Natural logs everywhere; `logf`, `logg`, `fhat` are log-densities; NLopt minimises `-Q` | `M_step.cpp:41-62` |
| `f` normalisation | Labelled per-lineage density: `log λ` per speciation event, no `(n-1)!`, no `2^{n-1}`, no binomial coefficient on `rho` terms | `model.hpp:433-438, 466` |
| `q` normalisation, thinning | Includes `-log(2·tips + Ne)` per augmented lineage; compensator of `n·λ·(1 - ρe^{-μ(T-t)})` | `model.hpp:295-308` |
| `g` normalisation, BDI | Per-lineage: `log(λu)` per birth and per immigration (not `log(nλu)`, not `log(2kλu)`), `log(μ/u)` per death, minus integrated total rate; rejection probability omitted | `bdi.R:440-506, 514` |
| Weight accepted (thinning) | `isfinite(lw) && exp(lw) > 0` — absolute threshold at `lw ≈ -745` | `E_step.cpp:94` |
| Weights passed to M-step | Thinning: `exp(lw - max)`, unnormalised, `Σw ∈ [1, N]`; BDI: normalised to mean 1, `Σw = N`; CEM final: softmax normalised | `E_step.cpp:133-139`; `bdi.R:796-799`; `de.R:837-838` |
| `fhat` denominator | C++: `N + rejected_zero_weights`; BDI: `n_valid`; CEM shared mode: `|Z|`; overrun/lambda rejections never in the denominator | `E_step.cpp:143`; `bdi.R:699`; `de.R:355` |
| "loglik", MCEM | `fhat(θ_{K-1})`: IS estimate with `N = sample_size` (200) trees drawn at the iterate **before** the last M-step, unconditioned in the table; `estimate_rates` subtracts `log P_tree(θ_K)` if `cond` | `inference.R:410-412, 819-822` |
| "loglik", CEM | `best_IS$fhat`: fresh IS estimate with `max(20, sample_size)` trees at the softmax-weighted mean of the final elites, minus `log P_tree` inside `emphasis_cem` | `de.R:852-870`; `inference.R:457-463` |
| "loglik", GAM | Predicted value of the additive Gaussian smooth at the L-BFGS-B optimum; grid points each used 200 trees and had `log P_surv` subtracted; no IS evaluation at the returned pars | `gam.R:967`; `inference.R:527-530, 569` |
| `loglik_var` | Bootstrap (B=200) variance of the K=2 moment-corrected estimator `mean(lw) + log(1 + m₂/2)`, no zero-weight term; GAM: `NA` | `de.R:476-489`; `emphasis.R:136-142` |
| `AIC` | `-2·loglik + 2·length(compact pars)`, bound-fixed parameters counted | `inference.R:816, 824` |
| Survival conditioning | `cond` off by default; `P_surv` from a binomial GAM trained on forward simulations with `status == "done"` (both crown clades alive at `T` and `N < max_lin`) | `inference.R:678`; `gam.R:37`; `general_tree.hpp:175, 290` |
| RNG | C++ engines seeded from clock ⊕ thread id (`std::default_random_engine` in `augment_tree.cpp:41`, `mt19937_64` in `model.hpp:218`, `mt19937` in `general_tree.hpp:34-43`); R's RNG only for rho tip-dropping and `auto_bounds` Phase 3 | `simulate.R:218-221`; `gam.R:318-341` |
| Threads | `num_threads` sets grainsize only in `E_step`/`augment_trees` (arena constructed, not used); `M_step` uses `task_scheduler_init` | `E_step.cpp:62-65`; `augment_tree.cpp:221-224`; `M_step.cpp:77` |

---

## 3. README/docs vs code

Deduplicated across all maps. "README" = `/Users/pancho/Code/emphasis/README.md`.

**Model, covariates, links**
1. README:79-84 and :103, `simulate.R:9-11`: `λ = f(β0, βN N + βD D)`, M "used only internally … not a user covariate". Code: `βM·M` always in the linear predictor (`model.hpp:154,165,200,211`; `general_tree.hpp:126,138`); `model_bin[1] = use_M`; `.resolve_model` accepts `c(1,1,0)`/`c(1,1,1)` (`simulate.R:268-274`); `.par_names`/`.pars_error_msg` expose `beta_M/gamma_M` (`inference.R:240,244`; `simulate.R:313-324`); `test-covariates.R:40-42` asserts those names.
2. README: `D_s = I_s - M`, `I_s` "time since it last speciated". Inference uses `tip_start = 0` for observed lineages and never resets them (`E_step.cpp:34`, `augment_tree.cpp:50`); observed speciation nodes use `D = 0` (`model.hpp:184-185`); the simulator resets the parent (`general_tree.hpp:264-267`).
3. README:86-94 describes three links; `R/RcppExports.R:9,13,36-37,65-66,97-98,131-132`, `src/div_tree.cpp:10-14`, `src/rcpp_mce.cpp:31`, `src/rcpp_mcem.cpp:39`, the generated `man/*_cpp.Rd`, `augment_trees.Rd`, `eval_logf.Rd` document `model = c(use_N, use_P, use_E)` and link `0/1` only. `general_tree.hpp:87-88` comments still say `beta_E`, `{use_N, use_P, use_E}`. Tests' comments carry the same layout (`test-de.R:92-97`, `test-em.R:8`).
4. README says `β0` is the gaussian peak rate; `auto_bounds`/`.wide_bounds` put the gaussian intercept on a log scale (`gam.R:170-176, 434-440, 676-688`).
5. README:114 lists only `"ep"` as legacy alias for `"d"`; `.resolve_model` also accepts `"rd"` (`simulate.R:261`; `test-covariates.R:13`).

**Samplers and IS**
6. README:124 "weight is `f_θ/q`"; evaluated `q` contains `-log(2·tips + Ne)` per augmented lineage (`model.hpp:304,308`), which is not the sampler's uniform parent-choice probability over `tips-2+Ne` nodes (`augment_tree.cpp:154-166`).
7. README: thinning "every model and link" / "Poisson-thinning proposal". Envelope is the max of two endpoint evaluations, not a bound, and acceptance probability is not clipped (`augment_tree.cpp:138-151`); `q` matches the sampled density only when M and D coefficients are 0 (`model.hpp:221-252, 272-296`).
8. README:132, `bdi.R:6`: BDI "exact … ESS = sample size" under constant rates. Verified for `λ0 > μ0` only; for `μ0 > λ0` the hazard returns `Inf` and `uniroot` errors (`bdi.R:59, 88`).
9. README:135: gate excludes D-terms and gaussian; nothing excludes `rho < 1` (`bdi.R:100-105`; `inference.R:359`; `simulate.R:345`).
10. `simulate.R:92-95` (roxygen `rho`): augmentation inserts unsampled extant lineages — only true for `method = "thinning"`; BDI never emits `t_ext == 5e10`.
11. `inference.R:788-790` prints "BDI (exact)" for any supported model, including `dd`, which README:132 and `bdi.R:404-405` describe as approximate.
12. `bdi.R:26-32` comment states `∫μ/(1-p) = λΔt + ln[(1-p1)/(1-p2)]`; code uses `ln[(1-E1)/(1-E2)]` (`bdi.R:60`; code is the valid identity). `bdi.R:40-41` comment gives the extinction probability as `p(t)` for the critical case.
13. `bdi.R:713-715`: `.mcem_bdi` is "a drop-in replacement" for `.mcem_dynamic_fresh`; denominators differ (`bdi.R:699` vs `E_step.cpp:143-144`), rejection counts are always 0 (`bdi.R:804-807, 851`), tree formats differ (`pd ≡ 0`, `bdi.R:605`).
14. `bdi.R:792-795` claims the thinning M-step uses log-weights; `E_step.cpp:132-139` converts to linear max-scaled weights. `rcpp_mcem.cpp:54` documents `weights` without saying they are linear max-scaled; `rcpp_mcem.cpp:50` documents `fhat` as a vector (scalar at `E_step.cpp:144`, `rcpp_mcem.cpp:114`).
15. `rcpp_mcem.cpp:28`: `maxN` = "maximum number of failed trees"; code counts every attempt (`E_step.cpp:65`), matching `inference.R:60` and `rcpp_mce.cpp:26`.
16. `emphasis.hpp:113`: `max_time_seconds = 0` means "no time limit"; `E_step.cpp:51,71-73` implements 0 as 120 s. README and R docs never mention the 120 s cap, the zero-weight-in-denominator convention (`E_step.cpp:141-144`), or that `set.seed()` does not reach any C++ RNG.
17. README:61,133: "parallel augmentation via TBB" with `num_threads`; `E_step` and `augment_trees` never run inside the arena they construct (`E_step.cpp:63-65`; `augment_tree.cpp:221-224`).
18. `rcpp_mce.cpp:28`, `simulate.R:87`: `max_lambda` bounds "the speciation rate"; code bounds the total thinning intensity `n·λ·(1-ρe^{-μ(T-t)})` (`augment_tree.cpp:141`). Three defaults: 500 (`simulate.R:139,468`; `emphasis.hpp:40`), 1e6 (`emphasis.R:50,182`).
19. README:168: rho adds "the binomial sampling factor `n_obs log ρ + n_unsamp log(1-ρ)`"; code has exactly those two terms and no binomial coefficient (`model.hpp:466`).
20. Docs site `index.html` describes only the thinning proposal `ν(t) = Nλ(1-ρe^{-μ(T-t)})`; README makes BDI the default and never appears in the site.

**Inference driver and control**
21. README calls the package "Monte Carlo conditional-likelihood inference", `ℓ_cond = ℓ - log P_θ(survival)`. `cond = NULL` by default (`inference.R:678`); when on, the MCEM M-step adds `log P_tree` with coefficient 1 to an unnormalised weighted sum with `Σw ∈ [1,N]` (`M_step.cpp:45-61`; `E_step.cpp:139`).
22. `inference.R:53-54` roxygen: sampling "Currently only `dynamic_fresh` is implemented"; default is `"bdi"` (`inference.R:111`), README:128-135 says BDI is default.
23. `inference.R:75` roxygen: CEM `max_iter` default 20; code 50 (`inference.R:121`). CEM `tol/patience/sd_decay/sd_min_frac/max_time/rho` and the whole GAM control block are undocumented.
24. README:216-223 MCEM example: `sample_size = 1`, model `dd`, link `gaussian` — not BDI-supported, so it runs thinning MCEM with one tree per M-step and `loglik_var = NA`; README does not say so.
25. README:152: MCEM "doubles augmentation attempts and perturbs" on E-step failure. True for `.mcem_dynamic_fresh` (and the doubled `maxN` is discarded after the next success, `emphasis.R:89`); the default BDI path only perturbs (`bdi.R:778-788`).
26. `estimate_rates` `@return` (`inference.R:650-653`): method `mcem` or `cem`, details from `.mcem_dynamic_fresh` or `emphasis_cem`; `gam` and `.mcem_bdi` also occur.
27. `estimate_rates` `@param init_pars` (`inference.R:616-617`): "ignored for cem"; also ignored for `gam`, and `.validate_linear_init_pars` still runs for any method under the linear link (`inference.R:734-738`).
28. README:212 uses `control = list(n_grid = 200)`; `estimate_rates_control("gam")` documents `grid_points` (`inference.R:137`); `n_grid` is read as an undocumented extra key (`inference.R:484`; `pipeline.R:168`).
29. `.is_fhat` roxygen (`de.R:41-43`): K=2 correction "equivalent to the log-normal correction `μ + σ²/2`"; code computes `ℓ̂ + log(1 + m₂/2)`.
30. `emphasis_cem` `@return` (`de.R:558-559`): `converged ∈ {annealing, plateau, max_iter}`; code also returns `time_budget`, `all_failed` (`de.R:672, 722`). `de.R:556-557` and `inference.R:767-770`: `loglik_var` NA if `sample_size = 1`; the final evaluation always uses `max(20, sample_size)` trees and bootstraps (`de.R:800, 881-885`).
31. `compare_models` roxygen example (`inference.R:886-890`) calls `estimate_rates(..., lower_bound=, upper_bound=)`; no such formals (`inference.R:672-678`).

**Pipeline, bounds, GAM**
32. README "Stage 1": "a second bisection pass tightens bounds where IS augmentation fails". `.diagnose_is_bounds` prints warnings only (`gam.R:487-533`) and runs only when `verbose = TRUE` (`gam.R:363-367`).
33. README "Stage 2": GAM loglik described as an IS log-likelihood; code reports the smooth's predicted value at the optimiser's point (`gam.R:967`; `inference.R:569`).
34. README "The best result across stages is returned": first finite-loglik fit in fixed order `mcem > cem > gam` (`pipeline.R:289-302`); MCEM init is "cem if finite else gam" (`pipeline.R:231-239`).
35. README pipeline sample output (:186-197) has stages `bounds/gam/mcem`, no `cem` row although default `stages` include `cem` (`pipeline.R:60`), and the `*` marker sits on `bounds` while "Best stage: mcem" (`print.emphasis_pipeline`, `pipeline.R:346-354`).
36. README "set rho in the control list": `emphasis_pipeline` forwards `control$rho` to `auto_bounds` only (`pipeline.R:125` vs `168-175, 198-205, 249-257`).

**Simulation API**
37. `simulate.R:46-47`, README: batch `survival_prob` = "fraction of simulations with `status == done`"; code returns `mean(1/n_attempts)` (`simulate.R:168, 211`); measured 0.383 vs 0.467 with `max_tries = 1`.
38. `simulate.R:66-68`, README:75: `tree` may be a numeric branching-time vector; conditional path calls `tree$L` on it (`simulate.R:404` via `343`) and errors.
39. `simulate.R:60-61`: "failed draws are NULL" in `trees`; C++ excludes failed draws (`E_step.cpp:96-104`); NULLs arise only from `useDDD = FALSE` or `L2phylo` failure (`simulate.R:373-378`).
40. README:168 and `simulate.R:92-95`: forward incomplete sampling; implementation marks dropped tips `L[,4] <- max_t` (`simulate.R:215-223`), giving negative edges in `tas` (min edge −4.99 at `rho = 0.5`).
41. `RcppExports.R:12` documents `max_tries` for `simulate_div_tree_cpp`; R always passes 0 (`simulate.R:207`).

**API surface, package metadata, tests, CI**
42. README:71, :291-299 and `index.html` list `compare_models()` and `select_diversification_model()` as public; both are `@keywords internal` (`inference.R:893, 1010`) and absent from `NAMESPACE`; `index.html` calls `select_diversification_model(bird.orders, link = "gaussian")` directly. README:51 says nine exports (NAMESPACE agrees), so the function-reference table is the inconsistent part.
43. README:50: `.expand_pars/.contract_pars` in `R/simulate.R`; `.contract_pars` is in `R/inference.R:347-352`.
44. README:3 badge URL `actions/workflows/R-CMD-check.yaml`; the only workflow is `.github/workflows/main.yml` (`name: R-CMD-check`).
45. README:41 "CI runs on every push"; `main.yml:3-12` triggers on pushes to `main/syst_biol/develop` and PRs to `main/develop`.
46. README:43 "The 55 remaining tests pass": 55 = all `test_that` blocks including 18 unconditionally skipped; 37 run (98 expectations). README:43 names `test-em.R`, `test-inference.R`, `test-simulate.R` as dead and says their bodies "target earlier APIs (mc_loglik(), estimate_rates(lower_bound=), positional simulate_tree())". That holds for `test-augment.R:10`, `test-em.R:12,34-35,47-48,60-61,72-73`, `test-inference.R:74,88,101`; it does not hold for `test-simulate.R:38,83,94,107` or `test-de.R:48-58,77-82,102-107,121-123,137-150`, whose bodies match current signatures and are dead only because of `skip()`. `test-de.R` and `test-augment.R` are omitted from the README's list.
47. README:307 "no vignettes" (true); `man/emphasis-package.Rd:29-30` points to `vignette("Simulation")` and `vignette("documentation_estimation")`, its example (:35) calls `emphasis(bird.orders, model = list(pars = ...))` (no such function), and its `\details` lists `augment_trees()/eval_logf()` as user workflow. `emphasis.Rcheck/00check.log:121-128` is from an older tree with vignettes.
48. `DESCRIPTION:11` `License: file LICENSE`; README:4,320 says MIT.
49. `man/diagnose_mcem.Rd` `\usage` is wrapped/aligned unlike every other generated Rd; README:66 says `man/` is generated.
50. `print.mcem_diagnostics` hard-codes "(bootstrap, B=200)" (`diagnostics.R:579-580`); `diagnose_cem`/`diagnose_mcem` are documented to return invisibly (`diagnostics.R:69, 354`) but return visibly (`305-315, 533-541`); `@return` omits `stop_reason` and the `rejected` column.
51. README:305 "main is identical on the forge and on GitHub"; the checkout is reported as not a git repository by the harness, so not checkable offline.
52. The installed library 0.4 lacks the `.bdi_supported` fallback present at `simulate.R:345`; README describes the source.

---

## 4. Audit hypotheses, ranked

Severity order within each theme: critical → major → minor → question. "Test" gives a concrete procedure. Merged concerns cite every original anchor.

### (a) Likelihood / importance-sampling correctness

**H1 — critical — `src/E_step.cpp:94`.** The acceptance test `isfinite(log_w) && 0 < exp(log_w)` is an absolute threshold: a tree with `lw < ≈ -745` is counted as zero-weight even though `exp(lw - max)` is representable after max-shifting. `lw` contains the observed-speciation `Σ log λ_i` and the `-∫n(λ+μ)` terms that `q` does not cancel, so `|lw| = O(n_tips)` (e.g. 400 tips at `λ ≈ 0.2` give `Σ log λ ≈ -640` before integrals). On large clades every augmentation is rejected, `E_step` throws "maxN exceeded … zero weights", and `.mcem_warn_estep` (`emphasis.R:198-209`) attributes it to `λ = 0`. `+Inf` weights (`logg = -Inf`, proposal support failure) are also swallowed into `rejected_zero_weights`. *Test:* `augment_trees` on a simulated CR tree of 300–500 tips with `λ = 0.2, μ = 0.05`; record `rejected_zero_weights` vs `sample_size` and the `logf/logg` magnitudes of survivors; then compile with the test replaced by `isfinite(log_w)` (keeping `-Inf` as zero weight) and confirm `fhat` becomes finite and matches `DDD::bd_loglik` within IS error.

**H2 — critical — `R/bdi.R:687` (+ `inference.R:359`, `simulate.R:345`).** BDI ignores `rho`: the conditioned process, `p(t)`, the rates and the absence of unsampled-extant lineages assume complete sampling; `rho` is forwarded only to `eval_logf`. The proposal has zero mass on augmentations with unsampled extant lineages, which have positive `f`-mass when `rho < 1`, so the IS estimator targets `f(obs, no unsampled)`. Measured BDI `fhat = -22.41 = C + 8·log 0.5` vs thinning `-15.41` at `rho = 0.5` on the same tree. Neither `.bdi_supported` nor the two dispatch sites gate on `rho`, and `bdi` is the default. *Test:* `estimate_rates(tree, "cr", control = list(rho = 0.7, sampling = "bdi"))` vs `sampling = "dynamic_fresh"`: logliks must agree within IS error; assert `.bdi_supported` returns FALSE (or a rho-aware BDI is implemented) when `rho < 1`.

**H3 — critical — `R/de.R:340` (shared-tree CEM mode).** For a tree `z_i` drawn from `q(·|θ_k)`, the weight at particle `θ_j` is computed as `logf(z_i, θ_j) - logg(z_i, θ_j)`; the density actually used to draw `z_i` (`q(·|θ_k)`, or the pooled mixture) never enters. `n_zero_weight` is dropped (`de.R:355`) and `bias_correct` forced FALSE. The docstring (`de.R:294-300`) asserts the reverse. Default is `shared_trees = FALSE`. *Test:* on a small CR tree run Mode 1 and Mode 2 with the same particle set and large `sample_size`; compare per-particle `fhat` to `DDD::bd_loglik` — Mode 2 should show θ-dependent bias; replace the denominator with the mixture `log((1/n)Σ_k q(z|θ_k))` and re-check.

**H4 — major — `src/M_step.cpp:61` (both MCEM paths).** With a conditional, the objective is `-Σ_i w̃_i log f_i + log P_tree(θ)` where `w̃` are unnormalised (`Σw̃ ∈ [1,N]` thinning, `= N` BDI). The self-normalised conditioned EM objective is `Σ_i (w̃_i/Σw̃) log f_i - log P_tree`, i.e. the penalty needs coefficient `Σw̃`. The comment at lines 58-60 states the normalised form. Consequence: MCEM with `cond` returns approximately the unconditioned maximiser while `estimate_rates` reports a conditioned loglik (`inference.R:820-822`), and the two samplers condition with different effective strength. *Test:* fit CR with `cond` supplied and `N = 200`; compare pars against the maximiser of `DDD::bd_loglik(cond = 1)` (survival-conditioned) and `cond = 0`; the current code should sit near the `cond = 0` optimum. Unit test: with constant `log f`, the M-step argmin must be invariant to rescaling `w`.

**H5 — major — `inst/include/model.hpp:304, 308`.** `q` charges each augmented lineage `-log(2·tips + Ne)` where `tips` starts at 2 and increments per observed node; the sampler's parent choice is uniform over `|alive_ids| = tips - 2 + Ne` nodes (`augment_tree.cpp:154-166`) and the lineage count is `n_i = tips + Ne`. None of the three coincide, and `f` is per-lineage while `q` uses `n·λ`. The cited "tech report" (`E_step.cpp:142`) is absent (`dev/audit`, `dev/validation` empty). Map 4's single-tree check nevertheless found thinning `fhat` within 0.02 nats of the Nee likelihood (§1.11 item 1). *Test:* on a 2-tip and a 3-tip CR tree, enumerate augmentations with exactly one missing lineage and compute `Σ_z q(z)` by numerical integration over `(t*, t_ext)`; it must equal the total probability of the one-lineage event under the thinning process. Then check `E_q[f/q]` against `DDD::bd_loglik(btorph = 1, cond = 0, soc = 2)` for several `(λ, μ)`; a θ-dependent gap falsifies the weight.

**H6 — major — `model.hpp:221, 243, 282-296` (M/D models: three `mu`s, two `lambda`s, two compensators).** For any model with `γ_M, γ_D` or `β_M, β_D ≠ 0`: (a) `extinction_time` draws with `extinction_rate` (no D) on a node whose `pd = 0` during augmentation; (b) `nh_rate` uses `extinction_rate_ep` with `pd(t)`; (c) `sampling_prob` credits `mu_i = extinction_rate_ep` on the finished node with its final `pd/focal_tip_start`. Likewise `nh_rate`'s `lambda` is evaluated on a copy of the proxy node (`focal_tip_start = 0`, `D = 0` if observed) while `sampling_prob` uses the finished node. For D + exponential, `sampling_prob`'s compensator is the exact per-lineage sum `S_b·exp_integral` while the sampler's intensity is `n·λ(proxy)`. The cancellation `f/q` derived for CR holds only when all coincide. *Test:* simulate a `d`/`nd` tree under the exponential link; estimate `p(y|θ)` by (i) thinning IS with large `N` and (ii) an independent brute-force estimate (forward simulation with ABC-style matching of branching times on a 4–5-tip tree, or exact enumeration for 2–3 tips); a θ-dependent discrepancy confirms `q ≠` sampler density. Also assert `nh_rate(t*)` equals the `λ_i μ_i n_i` credited in `logg` for every augmented node by instrumenting the C++.

**H7 — major — `src/augment_tree.cpp:138-140, 148-151` (envelope).** `λ1 = nh(cbt)` uses the pre-event lineage count (lower_bound returns the node at `cbt` whose `n` is the count before the event), so just after a branching the true intensity `(n+1)λs(cbt)` exceeds `Λ* = max(nλs(cbt), (n+1)λs(next_bt))` and `pt > 1` is accepted with probability 1. More generally `Λ*` is the max of two endpoint evaluations, not a dominating rate (M(t), the survival factor, the gaussian link and `max(0,·)` are non-monotone), and after a rejection the envelope drops to `λ2` alone (line 188). The realised intensity is `min(nh, Λ*)` while `log q` integrates unclipped `nh`. *Test:* instrument the loop to count `pt > 1`; on CR trees with small `N` it must be non-zero. Compare the empirical distribution of first-event times on `[0, next_bt)` against `1 - exp(-∫nh)` by KS test. Fix candidates: evaluate the start-of-interval rate with post-event `n`, or reject/clip with a corrected density.

**H8 — major — `model.hpp:424` (D + linear `f`).** For D-models under the linear (and non-exact gaussian) path, `f`'s integral uses `n_i·(λ_focal + μ_focal)` — `n` times the focal lineage's rate — instead of `Σ_s λ_s(t)`; with `max(0,·)` and per-lineage `D` this is not the model's total intensity. `q` and `nh_rate` use the same approximation, so the sampler is self-consistent, but `fhat` estimates a different model from README:79-84, and AIC across links mixes an exact `f` (exponential) with an approximate one (linear). *Test:* on a small augmented tree with two lineages of different `D`, compute `Σ_s λ_s` by hand and compare with `n·λ_focal` used in `loglik`; state which model is being fitted.

**H9 — major — `model_helpers.hpp:225`, `augment_tree.cpp:50`, `E_step.cpp:34`, `model.hpp:320, 412, 443-449`, `bdi.R:315` (pendant-age convention).** Observed nodes carry `tip_start = 0`, so `P(t) = k_obs(t)·t + Σ_missing(t - ts)`: every observed daughter counts as alive since the crown, the two crown lineages count nothing, and observed lineages' pendant starts are never reset by their own speciations; `E = M` (so `D = 0`) for observed speciation events. The forward simulator uses true pendant ages (`general_tree.hpp:178, 264-267`). Within `model.hpp` the D+exp running sums take an observed lineage's `ts` as its birth time (lines 320, 449) while `M` and the gaussian path (line 412) take `ts = 0`. Augmented lineages with `parent_id == -1` are reset to `tip_start = 0` (numerically verified: node at 0.0448 gets `pd = 0.0448`). `bdi.R:315` uses the same `P_obs = k·t` (inert under the gate). Hence the `M`/`D` covariates that generate data in `simulate_tree` are not the covariates the estimator evaluates. *Test:* simulate one `nd` tree, keep its L-table, compute the true `M(t)` and `D_s(t)` at each observed branching from the L-table, and compare with `pd/n` and `e_s - M` returned in the augmented data frame (`augment_trees(..., sample_size = 1)` with no missing lineages, `max_missing = 0`); they must match if the covariate is what README describes.

**H10 — major — `R/bdi.R:514, 508` (DD rejection not in `log g`).** Under DD, draws with missing lineages alive at `tp` or overflowing `max_missing` are rejected, so the effective proposal is `q/P(accept|θ)`, but `-log P(accept)` is not added to `logg`. `fhat` is biased upward by a θ-dependent amount (measured +0.23 and +0.45 nats at two θ). `max_tries = 5·sample_size` then returns fewer trees with no counter (`bdi.R:680-681, 804-807, 851`). *Test:* for a `dd` tree, compare BDI `fhat` with `DDD::dd_loglik` (mapping in §5) across a grid of θ; the gap should track `-log(acceptance rate)` recorded from the sampler.

**H11 — major — `R/bdi.R:696` (+ `model_helpers.hpp:129-136`).** BDI performs no finite-weight filtering: with `dd` on the linear link, mean-field `N̂` proposals reach `N` where `λ(N) = 0`, `loglik` returns `-Inf` (and, per Map 4, sometimes `+Inf` via `log_sum::result`); 56/300 trees had `±Inf`, `fhat = NaN`, and `w_norm` (`bdi.R:796-799`) became `NaN`, passed to `m_cpp`. Maps 2 and 4 disagree on the sign produced by `log_sum` at `λ = 0` (§1.11 item 2). *Test:* call `eval_logf` on a hand-built tree with one speciation node at `λ = 0`; assert `-Inf`. Fit `dd` (linear) with pars `(1.5, -0.12, 0.4, 0)` via BDI and assert all `logf` finite or filtered, `fhat` finite, no `NaN` weights reach `m_cpp`.

**H12 — major — `R/bdi.R:59, 86, 88` (μ > λ region).** For `μ0 > λ0`, `1 - E2 ≤ 0 < 1e-300` makes `.bdi_integral_cr` return `Inf` for every segment with `n > 0`; `uniroot` errors ("f() values at end points not of opposite sign", reproduced); `.mcem_bdi` counts an E-step failure and perturbs toward the box centre, so MCEM cannot visit `μ ≥ λ`. On the final segment `t_hi = tp - 1e-14`, `U > H(t_hi)` (probability `~e^{-30n}`) also errors, and `1 - exp(-d·1e-14)` carries ~2 significant digits. *Test:* `simulate_tree(tree, pars = c(0.3, 0.5), method = "bdi")` must not error; `estimate_rates` CR with bounds spanning `μ > λ` on a tree whose `DDD::bd_ML` MLE has `μ ≈ λ` (small, old clade) should reach it.

**H13 — major — `R/gam.R:967, 888-894`; `inference.R:824`.** The GAM stage loglik is the additive smooth's predicted value at the L-BFGS-B optimum, not an IS estimate; an additive `Σ s(p_j)` cannot represent the `λ–μ` or intercept–slope ridges, so its optimum need not lie on the ridge. This value feeds AIC and README's model-comparison example. *Test:* on a CR tree, compare the GAM-stage pars/loglik to `DDD::bd_ML`; then evaluate `fhat` with 2000 trees at the GAM optimum and compare with the reported smooth value.

**H14 — major — `R/emphasis.R:140`, `R/de.R:883` (+ `de.R:476-489`; `inference.R:855-856, 917-942`).** `loglik_var` bootstraps the K=2 moment-corrected estimator `mean(lw) + log(1 + m₂/2)` without the zero-weight denominator, while the reported loglik is the log-mean-exp with `S = N + n_zero`. `print.emphasis_fit` and `compare_models`' pairwise test present `sqrt(loglik_var)` as the MC s.e. of loglik. *Test:* on one E-step's `lw`, bootstrap both estimators and compare variances under a heavy-tailed weight distribution (e.g. `dd` with low ESS); assert `compare_models` p-values change materially.

**H15 — minor — `R/de.R:819, 866` and `de.R:355`.** `bias_correct` affects only Mode-1 search-time `fhat`; the final elite re-evaluation and final estimate omit it, Mode 2 hard-codes FALSE, and the branch silently falls back when `n_zero_weight > 0` (`de.R:82-86`); untested (`test-de.R`). *Test:* unit test that `bias_correct = TRUE` changes `best_IS$fhat`; property test that `ℓ̂ + log(1 + Σ_{k≤K} m_k/k!) → log mean exp(lw)` as `K → ∞`.

**H16 — minor — `inst/include/model.hpp:434, 275, 250`.** `μ` floored at `1e-300` in `f`'s extinction term but `1e-10` in `q` and `nh_rate`; `λ = 0` in `f` gives `-Inf` (tree rejected) whereas `μ = 0` at a sampled extinction gives finite `f` with relative weight `~e^{-690}` (kept). *Test:* evaluate `eval_logf` on a tree with an extinction node at `μ = 0` and a speciation node at `λ = 0`; document both outcomes; decide one convention.

**H17 — minor — `R/de.R:112` (`.ess_from_lw`).** Non-finite log-weights are dropped before ESS, so `-Inf`-weight trees neither reduce ESS nor count in `n`. *Test:* `ess_from_lw(c(0, -Inf, -Inf))` — decide whether 1 or 1/3 is intended.

**H18 — question — `src/E_step.cpp:143`, `R/de.R:69`.** `S = N + rejected_zero_weights` treats `λ = 0` trees as legitimate zero-weight samples but excludes `augmentation_lambda`/`augmentation_overrun`, which are θ-dependent, so `fhat(θ)` is conditional on completion with a θ-dependent event. Sampling that stops at `N` successes is an inverse-binomial design; `(1/S)Σw` has `O(1/N)` bias, largest at `sample_size = 1` (`lw_1 - log(1 + n_zero)`). The counter also includes attempts finishing after `stop` (`E_step.cpp:104-107`). *Test:* on a CR tree with `max_missing` set low enough to trigger overruns, compare `fhat` with and without the excluded counts against `DDD::bd_loglik`.

**H19 — question — `model.hpp:466`.** `n_obs log ρ + n_unsamp log(1-ρ)` omits `choose(n_obs + n_unsamp, n_obs)`; `n_unsamp` varies across `z`, so the coefficient is not constant across the IS sample. The `(1-ρ)` per unsampled lineage cancels with `q` line 308. *Test:* on a 2-tip tree with `ρ = 0.5` enumerate augmentations with 0/1/2 unsampled lineages and check whether the labelled-history convention makes the coefficient unnecessary.

### (b) MCEM / CEM optimisation logic

**H20 — major — `R/inference.R:411` (+ `emphasis.R:44-61, 92, 99`; `inference.R:821`).** The reported MCEM loglik is `fhat(θ_{K-1})` (E-step before the final M-step) while `pars = θ_K`; `loglik_var` comes from the same E-step; the survival correction is evaluated at `θ_K`. Under `max_iter`/`time_budget`/`e_step_failure` stops (not surfaced) the lag is not small. *Test:* after a fit, run one extra `em_cpp` at `θ_K` with `N = 2000` and compare `fhat(θ_K)` with `fit$loglik`; assert the pipeline's stage log would change ranking.

**H21 — major — `R/emphasis.R:114` (+ `bdi.R` convergence rule).** Constant `N = 200`, no growth or averaging, convergence when `max_j |Δθ_j|/range_j < 1e-3` for 3 consecutive iterations. With `N = 200` the MC noise in `θ_k` is generally larger than `1e-3·range`, so "converged" can fire on three quiet draws or never; the returned `pars` is one draw from the MCEM stationary distribution. *Test:* repeat the same fit 20 times with fixed data; report the sd of `pars` and the distribution of `stop_reason`; compare `pars` spread with the distance to `DDD::bd_ML`.

**H22 — major — `R/emphasis.R:89, 64-85` (+ `inference.R:60-63, 365`).** Recovery doubles `maxN` and pulls `θ` 20 % toward `(lb+ub)/2`, but resets `maxN` to the original after the next success. When the failure cause is structural (`sample_size ≥ maxN`, high rejection, the 120 s cap), the loop alternates fail/success and every failure moves θ toward the box centre; after 8 failures the returned `pars` is the perturbed vector (`emphasis.R:69, 165`). `maxN ≥ sample_size` is documented, not enforced; `ctrl$maxN = 2000` by default so the `max(2000, 10·N)` fallback only applies when the user passes `NULL`, and `num_trees > 2000` cannot succeed. *Test:* fit with `sample_size = 1900, maxN = 2000` on a tree with ~10 % rejection; trace `pars` vs iteration and show drift toward the centre; unit test that `maxN < sample_size` errors.

**H23 — major — `src/E_step.cpp:73` (+ `emphasis.hpp:113`; `mcem.cpp:27,46`; `rcpp_mce.cpp:61`).** A 120 s wall-clock cap is the meaning of `max_time_seconds = 0`; no caller passes the argument; on timeout with `< N` trees the throw says "maxN exceeded" and R responds by doubling `maxN` and perturbing (which cannot help). The check is per attempt, so one augmentation inserting up to `max_missing` lineages is unbounded. Also affects CEM (`de.R:151` → `NA`, indistinguishable from augmentation failure) and GAM. *Test:* on a large tree with `num_threads = 1` and `sample_size = 200`, time one E-step; if > 120 s assert the failure is reported as a timeout, not "zero weights"; thread `ctrl$max_time` through.

**H24 — major — `R/de.R:736, 705-708` (+ `pipeline.R:198`).** With `sd_decay = 0.85`, `sd_min_frac = 0.01` the SD floor is reached only after 29 non-improving iterations, and the plateau stop requires the floor; the pipeline sets `max_iter = 20`, so CEM always stops on `max_iter` with SD ≥ ~4.6 % of `(ub-lb)/4`. *Test:* run the pipeline's CEM stage and assert `converged == "max_iter"` on every run; property test that the stopping rule can fire within `max_iter` for the default control.

**H25 — major — `R/de.R:839, 830-848`.** `obtained_estim` is the softmax-weighted mean of elite parameter vectors; elite `fhat`s are re-drawn with 20 trees and differ by units, so weights are near one-hot on a noisy winner or, on a ridge (`λ–μ`, intercept–slope), the mean lies off the ridge; the reported loglik is evaluated at that averaged point. *Test:* on a CR tree compare `best_IS$fhat` at `obtained_estim` with the best elite's `fhat` re-evaluated with the same tree count; on a ridge (`dd`) compare with `DDD::dd_loglik` at both points.

**H26 — major — `R/inference.R:126` (+ `pipeline.R:198`; `de.R:151, 410, 659`).** CEM default `maxN = 10` with pipeline `num_trees = 5`: a particle needs 5 valid trees in 10 attempts; failures become `NA` and are excluded from elite selection; only the all-NA case escalates. The rejection rate is θ-dependent, so the search is filtered toward easy-to-augment parameters. *Test:* record the fraction of NA particles per iteration as a function of `θ` (e.g. by `μ/λ`); it should be non-uniform; compare CEM estimates with `maxN = 10` vs `maxN = 200` on the same tree.

**H27 — major — `R/pipeline.R:295, 289-302` (+ `pipeline.R:231`).** Stage logliks are not comparable (GAM smoothed surrogate; CEM 20-tree log-mean-exp at an averaged point; MCEM 200-tree `fhat` at `θ_{K-1}`), yet they are printed side by side and the MCEM init is "cem if finite else gam" regardless of value — the justifying comment ("GAM with low sample_size can produce misleading optima") contradicts the defaults (GAM 200 trees/point, CEM 5/particle). *Test:* re-evaluate every stage's `pars` with one common 2000-tree `fhat` and compare with the logged values.

**H28 — major — `R/pipeline.R:145` (+ `gam.R:207-214`).** When `auto_bounds` finds no feasible centre it warns and returns the wide box with `survival_gam = NULL`; the pipeline logs the stage as "ok" and runs later stages unconditioned on the widest box; `result$cond` is FALSE but the log does not flag it. *Test:* force infeasibility (e.g. a 3-tip tree) and assert the log marks the downgrade.

**H29 — major — `R/pipeline.R:125` (+ `168-175, 198-205, 249-257`).** `control$rho` reaches `auto_bounds` only; `gam/cem/mcem` run with `rho = 1` unless nested per stage. *Test:* `emphasis_pipeline(tree, control = list(rho = 0.5))` and inspect `fits$mcem$details` for the `rho` actually used.

**H30 — minor — `R/emphasis.R:93` (+ `64-85`).** On E-step failure `prev_pars` is not updated and the streak is not reset (`next` before lines 93-122), so the next `delta_max` compares against the pre-perturbation iterate and a convergence streak can survive interleaved failures. *Test:* inject an E-step failure in a converging run and assert the streak resets.

**H31 — minor — `R/emphasis.R:166` (+ `inference.R:410-412`).** If every iteration fails, `mcem` is `NULL`, `iterations = nrow(NULL) = NULL`, `loglik = NA`; `stop_reason` and `iterations` live only in `details`, never in the fit or print method. *Test:* unit test on the all-fail path; add `stop_reason` to `print.emphasis_fit`.

**H32 — minor — `R/emphasis.R:178` (`.mcem_warn_estep`).** The diagnostic re-runs `em_cpp` with `rho = 1`, no conditional, `maxN = 200` fixed, and parses "<n> zero weights" thresholding at `n > 50` while the message says "/200". *Test:* trigger the warning under `rho < 1` and compare the diagnosed cause with the failing configuration's rejection breakdown.

**H33 — minor — `R/emphasis.R:50` (+ `emphasis.hpp:40`; `inference.R:140, 518`; `simulate.R:139, 468`; `augment_tree.cpp:141`).** `max_lambda` is 1e6 in MCEM, 500 in GAM/simulate/C++ default; not a control entry for MCEM; it bounds the total thinning intensity `n·λ·(…)`, not the per-lineage rate as documented. *Test:* run the same θ through GAM and MCEM E-steps on a 500-tip tree; `rejected_lambda` differs; document which quantity is bounded.

**H34 — minor — `R/emphasis.R:101`.** The per-iteration `rejected` column and verbose `rej=` use only `results$rejected` (unhandled exceptions), not overruns/lambda/zero-weight; a run with 90 % zero-weight trees shows `rej = 0`. *Test:* assert the trace column equals the sum of all rejection counters.

**H35 — minor — `R/inference.R:816` (+ `emphasis.R:29`; `inference.R:548`).** `n_pars` counts parameters fixed by `lb == ub`; AIC over-penalises by 2 per fixed parameter; `compare_models` rankings shift. *Test:* fit `dd` with `gammaN` fixed via `lb = ub = 0` and check `n_pars`.

**H36 — minor — `R/inference.R:290, 256-324, 740-749`.** `.validate_linear_init_pars` guards only the positive covariate extreme; `D` is centred and can be negative (to about `-crown age`); uses `≤ 0` for `λ` but `< 0` for `μ`; `β0 = 0` yields slope 0 with `λ ≡ 0`. *Test:* init `d` linear with `β_D > 0` at midpoint bounds and check whether the first E-step collapses.

**H37 — minor — `R/de.R:660, 677, 813, 857`.** Rescue escalation multiplies `max_missing` by 10 up to three times and by 1.1 per overrun iteration without bound, persisting into the final evaluation. *Test:* log `max_missing` per iteration in a run that triggers rescue.

**H38 — minor — `R/de.R:745`.** Verbose `sd%` divides by `sd_floor/sd_min_frac`, which is 0 on inactive slots, giving `NaN` for every model. *Test:* run any CEM fit with `verbose = TRUE`.

**H39 — minor — `R/gam.R:836` (+ `inference.R:151-166`; `gam.R:808-823`).** The projected-time check applies only on the sequential path of `estimate_likelihood_surface`; the `mclapply` path has no budget; in the pipeline the error discards the partial surface. *Test:* set a tiny `max_time` with `num_threads = 1` vs `> 1`.

**H40 — question — `src/M_step.cpp:81`.** `xtol` (1e-3) is `xtol_rel` relative to `|x|`; `tol` (1e-3) is relative to the bound range. For slopes near 0 the optimiser stops much later than `tol`; for `β0 ≈ 1` they coincide. Is the M-step tolerance meant to sit below the EM threshold in the same metric? *Test:* rerun MCEM with `xtol = 1e-6` and compare the iterate noise.

**H41 — question — `R/inference.R:714` (+ `model.hpp:127-130`).** For the gaussian link `init_pars` is the bound midpoint with no slope handling; the CR rate is `β0·e^{-1/2}` and slopes act inside a square. Is the midpoint a sensible start? *Test:* count `max_iter` stops for gaussian fits starting at the midpoint vs at the CEM point.

**H42 — question — `R/de.R:698, 418`.** Plateau detection compares maxima of freshly re-drawn noisy `fhat`s (1–5 trees per particle), so `improvement ≥ 1e-4` is dominated by MC noise. *Test:* re-run the same population twice and measure the sd of `best_loglik`.

**H43 — question — `R/inference.R:427, 548` (+ `de.R:382-389`).** Inactive 8-element slots have `lb = ub = 0` and `sd = 0`; are they guaranteed to stay at 0 through SBPLX and CEM perturbation? *Test:* assert all inactive slots are exactly 0 in every particle and every M-step estimate.

### (c) Augmentation samplers (structure, topology, parents)

**H44 — critical — `R/simulate.R:404` (`.extract_Ltable`).** For a `simulate_tree()` result, `tree$L` is the FULL L-table including true extinct lineages, and `.sim_tree_conditional` appends augmentations to it. Measured: 30-tip CR tree with 21 extinct rows → augmented tables averaged 70 rows (49 with `phylo2L(tes)`), 43/50 draws had a parent dying before its child was born, 43/50 `tas` had negative edges; 0/50 with `phylo2L(tes)`. *Test:* assert `min(tas$edge.length) ≥ 0` and `nrow(L_aug) == nrow(phylo2L(tes)) + n_augmented` for simulated inputs.

**H45 — major — `R/bdi.R:556, 551-552` and `R/simulate.R:420, 428, 440` (+ `augment_tree.cpp:153-166`).** Parent assignment: BDI uses the deterministic id of the most recent observed branching (id 0 — an event after the birth — for lineages born before `bt_sorted[1]`; line 557/552 dead code); thinning draws uniformly over event nodes (never the two crown lineages) and returns `-1` when none, and `.aug_to_Ltable` drops `parent_id == -1` rows (33/1265 = 2.6 % of augmented lineages, the oldest ones, absent from `tas` while contributing to `log_q`), attaches the rest to the DAUGHTER row of the parent event via nearest birth age, and falls back to `L_extant[1,3]`. `tas` topology is therefore a convention, not a draw; `tas` and `log_q` describe different `z`. *Test:* count `parent_id == -1` rows per draw; assert every augmented lineage in the C++ data frame appears in `tas`; for BDI assert no lineage is attached to a node younger than itself.

**H46 — major — `R/simulate.R:221-222` (rho forward).** Dropped tips get `L[,4] <- max_t` ("extinct at present" in the comment) which in DDD age units is extinction at the crown; `tas` has negative edges (min −4.99 at `rho = 0.5`), `tes` is unaffected; `n_drop` capped at `n_extant - 2` without guaranteeing a survivor per crown clade; these rows are inherited by H44. *Test:* `simulate_tree(pars = c(0.6, 0.1), max_t = 6, model = "cr", rho = 0.5)$tas` must have non-negative edges (expected fix: `L[drop,4] <- 0`).

**H47 — major — `inst/include/general_tree.hpp:37`, `src/augment_tree.cpp:41`, `model.hpp:218` (+ `de.R:259-263`).** Three independent thread-local engines seeded from the wall clock ⊕ thread id; `std::default_random_engine` is `minstd_rand` (31-bit LCG) on libc++; `set.seed()` reaches only rho tip-dropping and `auto_bounds` Phase 3. Two `set.seed(42)` forward runs gave 29 vs 156 L-table rows. Under `mclapply` forks the process-global engines' post-fork state is unknown. *Test:* two `set.seed(1)` runs of `simulate_tree` and of `augment_trees` must differ (documenting the fact) until a seed argument is added; after adding one, assert equality and that forked children differ.

**H48 — major — `src/augment_tree.cpp:50` (+ `model_helpers.hpp:216-229`; `simulate.R:420`).** Augmented lineages with `parent_id == -1` (born before the first non-crown branching) get `tip_start = 0` and `e_s = M` (`D = 0`), while their extinction node keeps `tip_start = t_spec`; `P(t)` over-counts them by `t_spec`; the D+exp running sums add `e^{-β_D·brts}` but subtract `e^{-β_D·t_spec}`. *Test:* augment a tree with `max_missing = 1` many times, select draws with `parent_id == -1`, and verify `pd` against a hand computation.

**H49 — minor — `src/augment_tree.cpp:167` (+ `model.hpp:216, 231`).** `insert_unsampled_species` fires whenever `ext_time ≥ b` regardless of `rho`; if `extinction_time` can return `≥ T` at `rho = 1` (the truncated exponential should prevent it), `tas` gains an extra extant tip. *Test:* 10⁴ augmentations at `rho = 1`, assert no `t_ext == 5e10`.

**H50 — minor — `inst/include/model_helpers.hpp:161`.** `trunc_exp` is rejection sampling with acceptance `1 - e^{-μr}`; with `μ` clamped to 1e-10 one call can spin ~1e9 iterations; the 120 s cap is the only backstop. *Test:* time `augment_trees` with `γ0 = 0` under the linear link; replace with inverse-CDF sampling.

**H51 — minor — `src/augment_tree.cpp:101`.** The `n`-recomputation loop has no `tree.end()` guard; safe only because the caller routes `ext_time ≥ b` elsewhere. *Test:* add an assert or bound; call `insert_species` directly with `t_ext ≥ T` in a C++ unit test.

**H52 — minor — `inst/include/general_tree.hpp:70-74, 85, 212-216, 250-256`.** Simulator times are `float` (~1e-7 relative), R matches birth ages by nearest value (`simulate.R:428`) and `prune_to_extant` uses `tol = 1e-8` (`inference.R:16-24`); with the linear link a selected lineage can have `λ + μ = 0` through floating fallback, giving `bernoulli(NaN)`. *Test:* simulate with `max_t = 100` and check `branching.times` for duplicates/zero-length edges; store `double`.

**H53 — minor — `src/augment_tree.cpp:57` (+ `model.hpp:242`).** `compute_pendant_pd` is `O(N²)` per tree and `nh_rate` rescans for `P(t)` at every candidate; dominates E-step time at large `max_missing`. *Test:* profile one E-step on a 200-tip tree with `max_missing = 1e4`.

**H54 — minor — `R/bdi.R:39-41`.** Critical-case branch (`|λ0-μ0| < 1e-15`): `I_lam` should be `λΔt - log(a1/a2)` and `I_mu` `λΔt + log((tp-t1)/(tp-t2))` (numeric 4.666 vs code 2.622); practically unreachable, and the general formula loses precision as `d → 0`. *Test:* compare `.bdi_integral_cr` with `integrate()` at `d = 1e-6, 1e-10, 0`.

**H55 — minor — `R/bdi.R:108-116` (+ `bdi.R:171-284, 338-366, 634-635`).** `.bdi_lam/.bdi_mu` read slots 3–4 as coefficients on `P` and `E` (pre-2026-07 layout) while C++ reads them as `β_M, β_D`; an 8-element `pars` with non-zero slots 3/4 is accepted and evaluated inconsistently; ~150 lines of Gaussian-closure code are inert under every reachable configuration. *Test:* pass `pars8` with `β_D ≠ 0` to `.augment_tree_bdi` and compare `logg` conventions with `eval_logf`; or assert the sampler errors.

**H56 — minor — `R/bdi.R:605`.** BDI trees carry `pd ≡ 0` and no `focal_tip_start`; thinning trees carry per-node `pd`; any M/D-aware scoring of BDI trees sees `M = 0`. *Test:* `eval_logf` on a BDI tree with `β_M ≠ 0` vs a thinning tree of the same topology.

**H57 — minor — `R/bdi.R:496, 491, 501-502` (+ `simulate.R:380`; `model.hpp:304`).** BDI `logg` uses per-lineage `log(λu)` per birth and immigration; thinning `log q` uses `log(nμλ) - log(2·tips+Ne)`; both are returned as `log_q` by `simulate_tree`. The zero-variance CR property holds because `f` drops the same multiplicity, but the two "log q" values are not the same density. *Test:* on one fixed augmented tree (convert a BDI draw to the thinning data frame), compute both `logg`s; the difference should be a constant independent of `θ` if they are compatible; if it depends on `θ`, only one can be paired with `f`.

**H58 — minor — `R/bdi.R:377`.** `.bdi_iterate` stops at 20 iterations with no non-convergence report; Cov ODE drops death cross-jumps; `Ê` clipped to `[0, tp]`. *Test:* return the final `delta` and assert it is below tolerance on a strong-DD tree.

**H59 — minor — `R/bdi.R:665, 680-681`.** Under CR `max_tries = sample_size`, so any overflow returns fewer trees with no retry and no counter; `num_trees` reflects it but `rejected` is 0. *Test:* set `max_missing = 2` and assert the shortfall is reported.

**H60 — question — `R/bdi.R:465`.** The DD branch breaks when `1-p < 1e-12` or `total < 1e-15` without the survival term for the remainder; with `n_alive > 0` this is a rejection at 514. Intended? *Test:* count how many rejections originate from this break vs genuine survivors.

**H61 — question — `R/bdi.R:546`.** For a 2-tip tree `seq(0L, -1L)` is assigned to a zero-length index. *Test:* run BDI on a 2-tip tree.

### (d) Numerical stability, concurrency, build

**H62 — major — `src/E_step.cpp:95, 101`.** `if (!stop)` is evaluated outside the mutex and `stop = (E.trees.size() == N)` uses `==` inside it; a second thread can push tree `N+1` and set `stop = false`, after which all threads continue to `maxN`/120 s; `S_completed` then uses `N` rather than `trees.size()`, so `fhat` is biased by `log(N/num_trees)` and the M-step runs on an unplanned sample. *Test:* run `em_cpp` with `num_threads = 8`, `sample_size = 50`, `maxN = 5000` repeatedly and assert `length(logf) == 50`.

**H63 — major — `src/E_step.cpp:63` (+ `augment_tree.cpp:221-224`; `M_step.cpp:77`).** `tbb::task_arena arena(num_threads)` is constructed but the `parallel_for` is not run inside `arena.execute`, so it uses hardware concurrency; `num_threads` only sets grainsize. R callers that set `cpp_threads = 1` for `mclapply` forks (`de.R:223-227, 262`) still get a full TBB pool per child. *Test:* `num_threads = 1` with `sample_size = 2000` and measure CPU utilisation / `top`; assert single-core.

**H64 — major — `src/E_step.cpp:62` (+ `augment_tree.cpp:223`).** `grainsize = maxN / num_threads` (integer division) is 0 when `maxN < num_threads` (CEM default `maxN = 10`, `num_threads` up to hardware concurrency, 32 on whitebox); TBB's contract requires `≥ 1`. *Test:* CEM defaults on a ≥ 12-thread machine; fix with `std::max<size_t>(1, …)`.

**H65 — minor — `src/rsbplx.cpp:144`.** Any negative `nlopt_result` throws, including `NLOPT_ROUNDOFF_LIMITED (-4)`, a benign termination for derivative-free methods; `em_cpp` then errors and the MCEM driver treats it as an E-step failure (doubling `maxN`, perturbing). *Test:* log `nlopt` codes across a fit; count `-4` occurrences; treat `-4` as success.

**H66 — minor — `src/M_step.cpp:84-87` (+ `rsbplx.cpp:108-123`).** No check that `lower/upper` have `pars.size()` elements before passing `.data()` to nlopt; compact bounds with 8-element pars read out of bounds. *Test:* call `m_cpp` with 2-element bounds under ASan.

**H67 — minor — `src/M_step.cpp:77, 76, 61`.** `tbb::task_scheduler_init` is the pre-oneTBB API (compiles only because RcppParallel 5.1.11.2 ships the header); the `is_threadsafe` guard is commented out; the R conditional (`Rcpp::Function`) is called inside the NLopt C call stack so an R error unwinds through C frames. *Test:* make `predict_survival` throw mid-optimisation and observe.

**H68 — minor — `src/E_step.cpp:97, 104-107`.** Each accepted tree is copied under the global mutex (serialising the parallel region); `rejected_zero_weights` increments after `stop` while accepted-after-stop trees are dropped uncounted. *Test:* compare `sum(counters) + N` with the number of attempts made.

**H69 — minor — `src/rcpp_mcem.cpp:11` (+ `rcpp_mcm.cpp:21`).** `rcpp_mcem` uses a private 3-column unpack (`brts, n, t_ext`) instead of `unpack.h`'s 8 columns; `m_cpp` rebuilds nodes with `pd = tip_start = 0`, `id = parent_id = -1`, so trees from `em_cpp(copy_trees = TRUE)` cannot be re-scored for M/D models and `m_cpp` is correct only for N-only models (currently the only caller, `bdi.R:820`). *Test:* `eval_logf` on `em_cpp` trees for an `nd` model vs the E-step's own `logf`.

**H70 — minor — `inst/include/emphasis.hpp:54` (+ `M_step.cpp:29-32`; `mcem.cpp:14-29`; `rcpp_mcm.cpp:54`).** Dead API: `E_step_info_t::logf/logg` scalars, `nlopt_f_data` destructor building an empty tree, `mce()` with ignored arguments and no caller, `sbplx_1`, the `plugin` argument of `m_cpp`, the vectorised `augment_trees` C++ function (not the R symbol). *Test:* remove and rebuild.

**H71 — minor — `inst/include/general_tree.hpp:231, 110-113` (+ `simulate.R` no validation).** Gaussian link with `β0 < 0` gives negative `λ`; if `λ + μ > 0`, `bernoulli(p ∉ [0,1])` is undefined behaviour; if the total is `≤ 0`, `expon()` returns 1e20 and the run ends `done` with 2 rows (confirmed for `pars = c(-0.5, 0.1)`). *Test:* assert `simulate_tree` rejects negative intercepts under the gaussian link.

**H72 — minor — `emphasis.Rcheck/00check.log:108-118`.** Compiled code contains `abort/assert` entry points in `E_step.o`, `loglik.o`, `rcpp_*.o`; an assert can terminate the R session. *Test:* `R CMD check --as-cran`.

### (e) API and cross-component consistency

**H73 — critical — `R/gam.R:170-179, 434-440, 676-688`.** `auto_bounds`, `.find_feasible_center` and `.wide_bounds` branch on `link_int == 0L` else log-scale, so the gaussian link (2) gets `β0 = log(lam_hat)` and bounds `[log lam_lo, log lam_hi]` although `model.hpp:127-130` uses `β0` on the natural scale; for `bird.orders` the centre is `β0 ≈ -2.3`. A zero/negative-rate centre passes `.test_feasibility` because a 2-tip tree satisfies `tip_lo = 2` (`gam.R:161`). The survival GAM is then trained mostly at non-positive rates. *Test:* `auto_bounds(bird.orders, model = "dd", link = "gaussian")` and assert `lower_bound["beta_0"] > 0`.

**H74 — major — `R/gam.R:415, 401-416, 374`.** The feasibility criterion (≥ 50 % of forward simulations `done` with tips in `[0.1N, 10N]`) declares high-turnover regions (`μ ≈ λ`) infeasible before any likelihood is evaluated, although survival conditioning is what makes those regions admissible; the survival GAM is trained only inside the box. *Test:* simulate CR trees with `μ/λ = 0.9`, run `auto_bounds`, and check whether `DDD::bd_ML(cond = 1)` lies inside the returned box.

**H75 — major — `R/simulate.R:343, 404`.** A numeric branching-time vector as `tree` (documented) errors in the conditional path (`tree$L` on an atomic vector). *Test:* `simulate_tree(tree = c(5,3,1), pars = c(0.5,0.1))`.

**H76 — major — `R/simulate.R:168, 211`.** `survival_prob = 1/n_attempts` per tree (values `{1, 0.5, 0}` at `max_tries = 1`), batch = its mean, not the `done` fraction; consumed by `auto_bounds`/`train_GAM`. *Test:* batch of 1000 with `max_tries = 1` vs `0`; compare to the empirical `done` fraction.

**H77 — major — `R/simulate.R:268-274` (+ `inference.R:240,244`; `simulate.R:313-324`; `bdi.R:104`; `test-em.R:23`).** The M slot is reachable through the binary-vector interface though README/roxygen call it internal; the formula path rejects `M`; `.bdi_supported` routes M-active models to thinning; C++ semantics of an M-active model are undocumented and untested. *Test:* decide; either reject `model_bin[2] == 1` in `.resolve_model` or document and test `beta_M`.

**H78 — major — `R/RcppExports.R:9,13,36-37,65-66,97-98,131-132` (+ `div_tree.cpp:10-14`; `general_tree.hpp:87-88`; `rcpp_mce.cpp:31`; `rcpp_mcem.cpp:39`).** Stale `(use_N, use_P, use_E)` / link `0/1` documentation on every C++ wrapper. *Test:* regenerate with corrected roxygen; grep for `use_P`.

**H79 — minor — `R/inference.R:101, 107, 384, 404` (+ `model.hpp:90`; `.mcem_warn_estep`; `select_diversification_model`).** `rho` unvalidated in R; C++ resets out-of-range values to 1 with no message; some callers do not forward it. *Test:* `estimate_rates(..., control = list(rho = 80))` must error.

**H80 — minor — `R/inference.R:529, 532` (+ `gam.R:871`; `inference.R:344`).** `.run_gam` subtracts `log(p_tree)` with no floor; `p = 0` gives `fhat = +Inf`, dropped by `is.finite`, so the lowest-survival grid rows are removed rather than penalised; `cond_fun` floors at 1e-300 so CEM/MCEM differ at the same points. *Test:* grid rows with `predict_survival == 0` and count how many vanish from the smooth.

**H81 — minor — `R/gam.R:561, 553-559`.** `.observed_covariates` sets `mean_pendant = mean(brts)` (all branching times), not a pendant age, for the M/D compensatory diagonals. *Test:* compare with the true mean pendant age from `tes`.

**H82 — minor — `R/simulate.R:248` (+ `general_tree.hpp:104-107`; `model.hpp:87`).** `link = 3` is linear in forward simulation and an out-of-range enum in augmentation. *Test:* `.resolve_link(3)` must error.

**H83 — minor — `R/simulate.R:348-364, 228-231, 377, 407`.** Every sampler and `L2phylo/phylo2L` call is wrapped in `tryCatch(error = NULL)`, so C++ exceptions, the 120 s cap, missing DDD and genuine bugs are indistinguishable from "augmentation failed". *Test:* pass an invalid `model_bin` and observe `tas = NULL, log_q = NA` with no message.

**H84 — minor — `R/simulate.R:374`.** With `useDDD = FALSE` the conditional path returns `tas = NULL` and discards the L-table and data frame. *Test:* request `useDDD = FALSE` and look for the augmentation.

**H85 — minor — `R/simulate.R:158, 168, 172`.** Batch recursion passes 15 positional arguments, forces `num_threads = 1` per row, `vapply(..., "survival_prob")` fails on `try-error` rows from `mclapply`, and the conditional batch return is undocumented. *Test:* make one row error under `mclapply`.

**H86 — minor — `R/diagnostics.R:142, 305-315, 374, 387, 457, 533-541, 579-580`.** `diagnose_cem` errors when `best_IS$ESS` is NULL; ESS taken from the fit rather than recomputed; documented invisible returns are visible; `par_names` from `names(x$pars)` may not match `par[0-9]+` columns; "B=200" hard-coded. *Test:* `diagnose_cem` on a fit lacking `ESS`.

**H87 — minor — `R/inference.R:933, 910`.** `compare_models` indexes AIC by label (two "CR" fits give a length-2 `T_ij`) and reads `n_pars` with `vapply(..., 0L)` (errors if double). *Test:* `compare_models(list(a, b))` with two CR fits.

**H88 — minor — `DESCRIPTION:11`, `.Rbuildignore:8`, `.github/workflows/main.yml:37`.** License field non-standard; `emphasis.Rcheck/` and `emphasis_0.4.tar.gz` not ignored; `checkout@v3`, no coverage job, no `--as-cran`, `syst_biol` trigger. *Test:* `R CMD build` and inspect the tarball; open the Actions page.

**H89 — minor — `R/gam.R:65-66, 907-908`.** `train_GAM`/`train_likelihood_GAM` `cat()` unconditionally, ignoring `verbose`. *Test:* `expect_silent`.

**H90 — minor — `R/simulate.R:475` (+ `inference.R:365`; `bdi.R:663`; `inference.R:126`).** Three attempt budgets for the same operation; `max_tries` ignored in the conditional path; failure statistics not comparable between `simulate_tree()` and `estimate_rates()`. *Test:* document and unify.

**H91 — question — `R/gam.R:37` (+ `general_tree.hpp:175, 290`).** `P_surv` is `P(both crown clades alive at T and N < max_lin = max(20N, 500))`, trained inside the detected box only; whether this is the conditioning event of the reported "conditional likelihood" (and the same event the M-step conditions on) is undetermined. *Test:* compare `predict_survival` with `DDD::bd_loglik`'s `cond = 1` normaliser (`1 - P(extinction)` squared) for CR.

**H92 — question — `src/div_tree.cpp:33`.** The C++ `max_tries` retry loop is never exercised from R. *Test:* remove or document.

**H93 — question — `R/simulate.R:380`.** See H57: two `log_q` definitions under one name.

### (f) Tests, CI, documentation

**H94 — critical — `tests/testthat/test-de.R:43,72,93,116,130`, `test-em.R:5,31,43,55,68`, `test-inference.R:72,86,99`, `test-simulate.R:35,80,91,104`, `test-augment.R:5`.** Every test reaching C++ (`augment_trees`, `eval_logf`, `emphasis_cem` end-to-end, `em_cpp`, `simulate_div_tree_cpp`) is unconditionally `skip()`-ed; the likelihood, proposal, M-step, MCEM, forward simulation, conditional term, multi-threading and the zero-weight denominator have no automated verification; CI green means "compiles and pure-R plumbing works". *Test:* re-enable `test-simulate.R` and `test-de.R` blocks (current signatures) first; add the closed-form checks of §5.

**H95 — critical — `R/bdi.R` (entire file, ~800 lines).** No test touches any `.bdi_*` function though BDI is the default proposal and README:132 claims exactness. *Test:* pin (i) `sd(lw) < 1e-10` and `ESS == N` under CR with `λ > μ`; (ii) `fhat == ` Nee crown likelihood computed by `integrate()` / `DDD::bd_loglik(cond = 0)` up to a constant; (iii) agreement with thinning under `dd`; (iv) behaviour at `μ > λ` and `rho < 1`.

**H96 — major — `tests/testthat/test-augment.R:10`, `test-em.R:12`.** `mc_loglik()` is undefined anywhere in `R/`. *Test:* rewrite against `augment_trees`/`em_cpp`.

**H97 — major — `tests/testthat/test-em.R:34-36, 47-49, 60-62, 72-74`.** `estimate_rates(..., lower_bound=, upper_bound=)` — no such formals. *Test:* move bounds into `control`.

**H98 — major — `tests/testthat/test-inference.R:74, 88, 101`.** `simulate_tree(c(0.5, 0.1), max_t = 5, model = "cr")` binds the pars to `tree`. *Test:* use `pars =`.

**H99 — major — README:92,183,246 (gaussian link untested).** No test in R or C++ exercises `link_int = 2`, the link every README example recommends. *Test:* add gaussian-link fits and a check that `β0` is the peak rate (`eval_logf` on a tree with `η_cov = 1`).

**H100 — major — README:168 (rho untested).** No test passes `rho`; C++ sampling terms (`model.hpp:225-228, 251, 282, 295, 307-308, 456`) are untested. *Test:* CR with `rho < 1` vs a brute-force estimate; assert the two samplers agree (H2).

**H101 — major — `R/pipeline.R:57`, `R/gam.R:133`.** `emphasis_pipeline`, `auto_bounds`, `print.emphasis_pipeline`, `diagnose_mcem`, `diagnose_cem`, `.run_*`, `.resolve_control_aliases`, `.validate_linear_init_pars`, `.build_cond_fun`, `estimate_likelihood_surface`, `.sim_tree_conditional`, `.aug_to_Ltable`, `prune_to_extant`, `select_diversification_model` have no tests. *Test:* smoke tests with `stages = c("bounds")` on a 10-tip tree.

**H102 — minor — `tests/testthat/test-covariates.R:52-58`.** The expand/contract round trip uses `mb = c(1,1,1)`, for which `.expand_pars` is the identity; only `c(1,0,0)` exercises the index arithmetic; `c(0,0,1)` and `c(1,0,1)` are uncovered. *Test:* add `d` and `nd` round trips with distinct values per slot.

**H103 — minor — `tests/testthat/test-gam.R:14` (+ `gam.R:67, 909`).** Tests guard on `mgcv` but package code calls `mgcv::gam` unconditionally with mgcv in Suggests. *Test:* run without mgcv.

**H104 — minor — `tests/testthat/test-inference.R:71-72, 85-86, 98-99`, `test-augment.R:4-5`.** `skip_on_cran()` precedes `skip()`, so the reported reason is "On CRAN" locally. *Test:* remove `skip_on_cran` or set `NOT_CRAN`.

**H105 — minor — documentation set (§3 items 1–52).** Every README/roxygen/Rd discrepancy in §3; in particular the badge URL, "55 remaining tests", the stale `emphasis-package.Rd`, and the two user-facing documents describing different default proposals (§3 item 20). *Test:* `devtools::document()` diff; a docs-vs-code checklist run against §3.

---

## 5. What a validation study can and cannot check

### 5.1 Constant-rates model (`model = "cr"`)

**Exact references.** `ape::birthdeath(phy)` maximises the Nee et al. (1994) likelihood of the branching times conditioned on the crown age, returning `d/b` and `b - d`. `DDD::bd_ML(brts, cond, btorph, soc = 2)` / `DDD::bd_loglik(pars1 = c(λ, μ, 0, 0), pars2 = c(tdmodel = 0, cond, btorph, soc = 2, ...), brts, missnumspec = 0)` gives the same likelihood with explicit control over conditioning and labelling. `DDD::dd_loglik` with `K → ∞` (or `pars1[3]` large) reduces to it.

**Mapping.** emphasis linear or exponential link, `cr`: `λ = max(0, β0)` / `exp(β0)`, `μ = max(0, γ0)` / `exp(γ0)`. Gaussian: `λ = β0·e^{-1/2}`, `μ = γ0·e^{-1/2}` (`model.hpp:127-130`).

**What agreement would establish.** With `cond = NULL` (default), `rho = 1`, `sampling = "bdi"`: Map 4 verified on one tree that BDI `fhat` equals the Nee crown-age likelihood (unconditioned on survival, labelled-history convention) to 1e-13 and that `lw` is constant; so `estimate_rates` reduces to SBPLX on an exact (up to the constant) log-likelihood with mean-1 weights. Agreement of `pars` with `DDD::bd_ML(cond = 0, btorph = 1, soc = 2)` would then establish the M-step, the EM loop, the bound handling and the parameter plumbing for CR — and nothing about the thinning `q`, D-covariates, or `rho`. Running the same with `sampling = "dynamic_fresh"` tests the thinning weight (H1, H5, H7) for CR only: `fhat` at fixed `(λ, μ)` vs `DDD::bd_loglik` across a grid separates weight bias (θ-dependent gap) from IS variance (ESS-controlled scatter). Because `f` has no `(n-1)!`/`2^{n-1}` factors, only differences of loglik across θ are comparable; the offset is constant.

**Confounders.**
- *Conditioning.* emphasis defaults to unconditioned; with `cond` the M-step under-weights the penalty by `Σw` (H4) and `P_surv` is a GAM of `P(both crown clades alive and N < max_lin)` (H91), so a `cond = 1` comparison with DDD tests the GAM and H4 together. Use `cond = NULL` vs DDD `cond = 0` first.
- *Crown vs stem.* emphasis is crown-conditioned (`soc = 2`, two lineages at `t = 0`); the crown is not a node and the crown lineages are never parents (H45).
- *Labelled vs unlabelled.* emphasis `f` is per-lineage labelled; DDD `btorph = 1` gives branching times (unlabelled) and `btorph = 0` the phylogeny; the difference is a θ-independent constant for CR but must be subtracted before comparing logliks (not needed for comparing MLEs).
- *rho.* DDD's `missnumspec` conditions on a known number of missing species; emphasis's `rho` is Bernoulli sampling of extant tips with an unknown number. No exact reference for emphasis's `rho` model in DDD; BDI ignores `rho` (H2). Validate `rho` only by brute force on tiny trees or by the thinning-vs-BDI cross-check after H2 is resolved.
- *IS variance and MCEM noise.* The returned `pars` is one draw around the MLE (H21), loglik lags one iteration (H20), `loglik_var` describes another estimator (H14). Compare distributions over replicate fits, not single fits, and compare `fhat` at the DDD MLE with `N = 2000`.
- *μ > λ.* BDI cannot visit it (H12); DDD can. Trees whose exact MLE has `μ ≈ λ` will show a BDI-specific discrepancy.
- *RNG.* Not seedable (H47); replicate fits are the only handle.
- *Bounds.* `auto_bounds` may exclude the DDD MLE (H74, H73); supply bounds by hand that contain it.

### 5.2 Linear diversity-dependent model (`model = "dd"`, linear link)

**Exact reference.** `DDD::dd_loglik(pars1 = c(λ0, μ0, K), pars2 = c(lx, ddmodel, cond, btorph, soc = 2, ...), brts, missnumspec = 0)` computes the likelihood of the reconstructed branching times under a process whose per-lineage speciation rate depends on the number of lineages alive in the full tree (`N`, including doomed lineages), by integrating the master equation over the number of unobserved lineages — the same `N` that emphasis uses (`node.n` counts all lineages, observed and augmented). `DDD::dd_ML` maximises it.

**Parameter mapping (derivable, with truncation caveat).**
- `ddmodel = 1` (linear in speciation, constant extinction): `λ(N) = λ0 - (λ0 - μ0)·N/K`, `μ(N) = μ0`. In emphasis's linear link `λ = max(0, β0 + βN·N)`, `μ = max(0, γ0 + γN·N)`: `β0 = λ0`, `βN = -(λ0 - μ0)/K`, `γ0 = μ0`, `γN = 0` (fix via `lb = ub = 0`; note H35 for AIC). DDD uses `pmax(0, ·)` for `N > K` as emphasis does, so the truncation matches; confirm in the DDD source for the installed version.
- `ddmodel = 1.3`: `λ(N) = λ0·(1 - N/K)` → `β0 = λ0`, `βN = -λ0/K`, `γN = 0`.
- `ddmodel = 3` (linear in extinction): `λ = λ0`, `μ(N) = μ0 + (λ0 - μ0)·N/K` → `βN = 0`, `γ0 = μ0`, `γN = (λ0 - μ0)/K`.
- `ddmodel = 2` (exponential in speciation): `λ(N) = λ0·(μ0/λ0)^{N/K}` → exponential link with `β0 = log λ0`, `βN = log(μ0/λ0)/K`, `γ0 = log μ0`, `γN = 0`. `ddmodel = 4` likewise for extinction.
- emphasis's general `dd` (both `βN` and `γN` free) has no DDD counterpart; validate on the `γN = 0` or `βN = 0` submodels. DDD's `lx` (lineage truncation) must exceed the largest `N` emphasis augments (`max_missing`).

**What agreement would establish.** `fhat` at fixed `(β0, βN, γ0)` vs `dd_loglik(cond = 0)` across a grid is the direct test of the IS estimator under diversity dependence — for thinning it tests H1/H5/H7 with a θ-dependent `N`; for BDI it tests the mean-field proposal plus the rejection bias (H10, whose signature is a gap tracking `-log P(accept)`) and the `±Inf` weights (H11). Agreement of MLEs with `dd_ML` then establishes the M-step and EM under `dd`. This is the only exact reference for a non-constant `N`-dependent rate; `d`/`nd` (D-covariate) models have no exact likelihood and can only be checked by simulation-based consistency (recover generating parameters, H9) or brute-force enumeration on 2–3-tip trees (H6).

**Confounders (in addition to §5.1).**
- The `max(0, ·)` truncation region: emphasis's augmented trees can reach `N` where `λ = 0`; DDD integrates the same truncation but emphasis rejects those trees (`-Inf` weight) or, under BDI, produces `NaN` (H11), so `fhat` and `dd_loglik` differ exactly where `K` is small relative to the augmented `N`.
- BDI's mean-field `N̂` (`bdi.R:307-392`) is approximate under strong DD; discrepancies there are proposal quality, not estimator bias, and show as low ESS and high rejection rather than a θ-dependent offset — unless H10 applies.
- `cond`: DDD `cond = 1` conditions on survival of both crown lineages under the DD process; emphasis's GAM `P_surv` includes the `max_lin` truncation (H91) and the M-step penalty scale (H4).
- `M` enters emphasis's `dd` linear predictor with coefficient 0 only because of zero-padding (H77); an M-active vector passed by accident changes the model with no message.

---

## 6. Open questions for the author

Only what the code cannot decide.

1. **Thinning `q` derivation.** What is the intended derivation of `-log(2·tips + Ne)` in `sampling_prob` (`model.hpp:304, 308`)? The cited tech report (`E_step.cpp:142`) is not in the repository. Is it a labelled-attachment count from the original emphasis paper, and does it still apply now that the parent is drawn uniformly from `alive_ids` (excluding the crown lineages)?
2. **Exact or heuristic `q` for M/D models.** For `d`/`nd`, is the thinning proposal meant to be exact (`log q` = the sampler's density), which requires closing the three-`mu`/two-`lambda`/two-compensator mismatches, or is an approximate `q` accepted by design (in which case the IS estimator is biased and README should state it)?
3. **Pendant-age convention.** Is `tip_start = 0` for observed lineages (`P(t) = k_obs(t)·t`, `D = 0` at observed speciation events) a deliberate mean-field definition of `M` and `D` for inference given that no topology is passed to C++, or a placeholder until topology is passed? The simulator uses true pendant ages, so "M" and "D" currently name different quantities in simulation and inference.
4. **Labelled-history convention.** Is the per-lineage `f` (no multiplicity factors) with BDI's per-lineage `g` and the thinning `q` intended to be the same target density, and is the omission of the binomial coefficient in the `rho` terms (`model.hpp:466`) intentional under that convention?
5. **Conditioned M-step.** Is the objective meant to be the self-normalised `Q - log P_tree` (penalty coefficient `Σw`), as the comment at `M_step.cpp:58-60` says, or is the unit coefficient deliberate? Should BDI (mean-1 weights) and thinning (max-scaled) share one weight convention?
6. **Which conditioning event.** Is the survival event "both crown lineages alive at `T`" with the `max_lin` truncation (`general_tree.hpp:175, 290`), and is that the event the reported "conditional likelihood" refers to?
7. **Reported loglik.** Should the MCEM loglik be re-evaluated at the returned `pars` (fresh E-step at `θ_K`, as CEM's `best_IS` does) rather than taken from the E-step at `θ_{K-1}`? Should `loglik_var` bootstrap the same estimator?
8. **MCEM design.** Is a constant Monte Carlo sample size with a parameter-stability stopping rule the intended design (with CEM/GAM as global search), or is sample-size growth / iterate averaging planned? What reference run (tree, seed, control) established `tol = 1e-3`, `patience = 3`, `sample_size = 200`, `maxN = 2000`?
9. **Recovery semantics.** What should happen on E-step failure: is the `maxN` reset after success intended, should perturbation toward the box centre apply when the cause is the 120 s cap, and should the cap be `ctrl$max_time` (header says 0 = no limit)?
10. **BDI scope.** Is `μ0 > λ0` deliberately out of scope for BDI (`bdi.R:59`), or should the formulas at lines 54/60 be used for `d < 0`? Should `.bdi_supported` return FALSE for `rho < 1`, or is a rho-aware BDI (with `p_ρ(t)` and an unsampled-extant event) planned? Is BDI-DD `fhat` meant to be a usable likelihood estimate (then the acceptance probability must enter), or only the M-step weights? Why is the `P̂/Ê` Gaussian-closure machinery kept when every reachable configuration zeroes its coefficients?
11. **`tas` from conditional simulation.** Is `tas` meant to be a faithful draw of the augmented tree (then `.extract_Ltable` must use an extant-only table, crown-interval lineages must attach to a crown row, and parent choice must be uniform over lineages) or only a visual aid?
12. **`fhat` denominator.** Are the counting conventions (zero-weight completions in, overrun/lambda out; inverse-binomial stopping) those of the paper's IS estimator definition?
13. **M slot.** Should `model_bin[2]` remain reachable through the binary-vector interface, and if so what are the intended C++ semantics of an M-active model?
14. **Shared-tree CEM mode.** Was it ever validated against independent mode on a known-answer tree, and is it still meant to be offered?
15. **Pipeline semantics.** Are cross-stage logliks in `result$log` meant to be comparable, and is "first finite in fixed order" the intended selection rule? Is CEM `maxN = 10` with `num_trees = 5` intentional?
16. **Reproducibility.** Is the clock-seeded C++ RNG accepted, or should the samplers and simulator take a seed from R (or use `unif_rand`) so `set.seed()` works and the skipped `set.seed(42)` tests can be pinned? Which engine is intended (`default_random_engine` vs `mt19937_64`)?
17. **Installed vs source.** Is the installed 0.4 library expected to match this tree (it lacks the `simulate.R:345` fallback)? Should audit findings be reproduced against source (`load_all`) or the install?
18. **Tests and exports.** Which currently-skipped blocks with current signatures (`test-simulate.R`, `test-de.R`) are meant to be re-enabled? Should `compare_models`/`select_diversification_model` be exported, as README and `index.html` present them?
# emphasis -- Package State

Date: 2026-03-07

## Module Map

### Module 1: Simulation (`R/simulate.R`)
| Function | Exported | Purpose |
|---|---|---|
| `simulate_tree()` | Yes | Unified entry point for forward/conditional tree simulation |
| `.resolve_link()` | No | Parse "linear"/"exponential" to integer |
| `.resolve_model()` | No | Parse "cr"/"dd"/"pd"/"ep" / formula / binary to `c(0,1,0)` |
| `.parse_model_formula()` | No | Parse `~ N + PD` formula |
| `.expand_pars()` | No | Compact pars -> 8-element C++ layout |
| `.pars_error_msg()` | No | Format parameter error message |
| `.sim_tree_conditional()` | No | Conditional simulation / augmentation |
| `.extract_tes()` | No | Extract extant-only phylo |
| `.extract_Ltable()` | No | Extract L-table from tree |
| `.aug_to_Ltable()` | No | Merge augmented branches into L-table |
| `%||%` | No | Null-coalescing |
| `.augment_tree_internal()` | No | Internal C++ wrapper |

Status: Clean. All helpers support `simulate_tree()`. No dead code.

### Module 2: Inference (`R/inference.R`)
| Function | Exported | Purpose |
|---|---|---|
| `prune_to_extant()` | Yes | Drop extinct tips |
| `estimate_rates_control()` | Yes | Default tuning parameters |
| `estimate_rates()` | Yes | Main estimation entry point |
| `print.emphasis_fit()` | Yes | Print method |
| `compare_models()` | Yes | AIC comparison table |
| `select_diversification_model()` | Yes | 3-model comparison pipeline |
| `print.model_selection()` | Yes | Print method |
| `.resolve_control_aliases()` | No | Sync old/new param names |
| `.extract_brts()` | No | Get branching times from tree/list/vector |
| `.par_names()` | No | `c("beta_0", "beta_N", ...)` |
| `.contract_pars()` | No | 8-element -> compact |
| `.run_mcem()` | No | Dispatch to `.mcem_dynamic_fresh` |
| `.run_cem()` | No | Dispatch to `emphasis_cem` |
| `.model_label()` | No | "CR", "N", "PD" label |

Status: Clean. Safe init_pars for DD/PD slopes. Proactive PD/EP warning.

### Module 3: MCEM Engine (`R/emphasis.R`)
| Function | Exported | Purpose |
|---|---|---|
| `.mcem_dynamic_fresh()` | No | MCEM EM loop with convergence |
| `.mcem_warn_estep()` | No | Diagnose E-step failures |

Status: Clean. Dead code removed.

### Module 4: CEM + IS helpers (`R/de.R`, `R/diagnostics.R`)
| Function | Exported | Purpose |
|---|---|---|
| `diagnose_cem()` | Yes | CEM convergence/IS diagnostics |
| `diagnose_mcem()` | Yes | MCEM convergence/IS diagnostics |
| `print.cem_diagnostics()` | Yes | Print method |
| `print.mcem_diagnostics()` | Yes | Print method |
| `emphasis_cem()` | No | Full CEM optimizer |
| `.is_fhat()` | No | IS log-likelihood from logf/logg |
| `.ess_from_lw()` | No | Effective sample size |
| `.simulate_particle()` | No | Augment trees at one parameter set |
| `.init_population()` | No | Initialize CEM particles |
| `.init_particle_grid()` | No | Grid initialization |
| `.eval_independent()` | No | Mode 1 IS evaluation |
| `.eval_shared()` | No | Mode 2 shared-tree evaluation |
| `.eval_particles()` | No | Dispatch evaluation |
| `.perturb_value()` | No | CEM perturbation |
| `.resample_particles()` | No | CEM resampling |
| `.bootstrap_fhat_var()` | No | Bootstrap variance |
| `.total_rejected()` | No | Sum rejection counts |
| `.n0()` | No | Null-safe integer |

Status: Clean. CEM is a substantial optimizer; all helpers are used.

### Module 5: GAM-based methods (`R/gam.R`)
| Function | Exported | Purpose |
|---|---|---|
| `train_GAM()` | Yes | Train survival-probability GAM (use case 1) |
| `predict_survival()` | Yes | Predict from trained GAM |
| `estimate_likelihood_surface()` | Yes | IS log-likelihood at parameter grid (use case 2) |
| `train_likelihood_GAM()` | Yes | Fit GAM to log-likelihood surface |
| `find_MLE()` | Yes | Optimize GAM to find MLE |

Status: Both use cases implemented and tested.
1. Survival conditioning: simulate unconditioned trees, fit GAM on survival
   indicator, use to condition likelihoods.
2. GAM-based MLE: simulate conditioned trees over parameter grid, estimate
   log-likelihood surface via GAM, optimize directly (no iterative EM/CE).

### Module 6: Legacy utilities (`R/generate.R`)
| Function | Exported | Purpose |
|---|---|---|
| `generatePhyloPD()` | Yes | Generate PD-dependent trees (old API) |
| `generateNonHomogeneousExp()` | Yes | Non-homogeneous exponential variates |
| `nhExpRand()` | Yes | Same, different algorithm |
| `rate_t()` | Yes | Compute rate at time t |
| `ExponentialRate()` | Yes | Compute exp rate from covariates |

Status: Likely removable. `generatePhyloPD()` superseded by `simulate_tree(model="pd")`.
The exponential utilities are standalone and unused by the inference pipeline.

### Removed files
- `R/augment.R` -- was empty (comment only)
- `R/utils.R` -- all 13 functions were dead code (none called outside utils.R)
- `get_required_sampling_size()` from `R/emphasis.R` -- dead code, allowed dropping `MASS`

## C++ Functions (src/)

| C++ function | R name | File | Exported | Purpose |
|---|---|---|---|---|
| `simulate_div_tree_cpp()` | `simulate_div_tree_cpp` | div_tree.cpp | No* | Forward tree simulation |
| `generateNonHomogeneousExpCpp()` | `generateNonHomogeneousExpCpp` | generateNonHomogeneousExp.cpp | No* | Non-homogeneous exponential |
| `eval_logf_cpp()` | `eval_logf` | loglik.cpp | No | Evaluate logf and logg |
| `rcpp_mce()` | `augment_trees` | rcpp_mce.cpp | No | E-step: augment trees |
| `rcpp_mcem()` | `em_cpp` | rcpp_mcem.cpp | No | Full E+M step |
| `rcpp_mcm()` | `m_cpp` | rcpp_mcm.cpp | No | M-step only |

\*Available via `emphasis:::` but not user-facing. `augment_trees`, `em_cpp`,
`m_cpp`, `eval_logf` were previously exported; made internal 2026-03-07.

## Augmentation Rejection Types -- Analysis

There are **two fundamentally different** ways an augmented tree can fail.
The IS estimator must handle them differently to remain correct.

### Type 1: Computational failures (overruns / lambda)

**Where**: `augment_tree.cpp:167-168` (overrun), `augment_tree.cpp:138` (lambda)

**What**: The thinning simulator aborted mid-construction because:
- `num_missing > max_missing`: too many extinct lineages were inserted
- `lambda_max > max_lambda`: speciation rate exceeded the thinning bound

**Nature**: These are truncations of the proposal distribution q's support.
The tree was never completed -- we have no logf or logg to compute. The
rejection probability depends on `max_missing` and `max_lambda` (user-chosen
computational limits), NOT on the model parameters theta.

**Handling**: Discarded from both numerator and denominator.

**Correctness**: Justified (see tech report Section 3). Since the rejection
probability is independent of theta, including these as zero-weight samples
would introduce bias proportional to the user-chosen limit, not the model.

### Type 2: Zero IS weights (parameter-dependent)

**Where**: `E_step.cpp:81-94`

**What**: The tree completed augmentation successfully, logf and logg were
computed, but `log_w = logf - logg` is non-finite (-Inf, NaN) or
`exp(log_w) == 0` (underflow).

**Root cause with linear link**: `lambda = max(0, beta_0 + beta_X * X) = 0`
for large X. Then `logf` includes `log(lambda) = -Inf`, so `log_w = -Inf`.

**Nature**: These are NOT computational failures. The proposal q successfully
generated a complete tree, but under the target f that tree has zero
probability. This is a legitimate IS outcome: w_i = f(z_i)/q(z_i) = 0.

**Handling** (FIXED): Zero-weight trees are now included in the IS
denominator. `S_completed = N_valid + rejected_zero_weights` is used in
both C++ (`E_step.cpp`) and R (`.is_fhat()`). They contribute w=0 to the
numerator and 1 to the denominator, which is the correct IS treatment.

**Previous bug**: Denominator used `N_valid` only, inflating fhat by
`log(S_completed / N_valid)` -- significant during IS collapse.

### Previous Recommendations -- Status

1. ~~**Track N_total alongside N_valid**~~ -- **DONE**. `rejected_zero_weights`
   is returned by C++ and passed through R. `S_completed` is computable.

2. ~~**Correct the fhat denominator**~~ -- **DONE**. C++ `E_step.cpp` and
   R `.is_fhat()` both use `S_completed` as denominator.

3. **Report acceptance rate in diagnostics** -- **OPEN**. Neither
   `diagnose_cem()` nor `diagnose_mcem()` reports `N_valid / S_completed`.
   The rejection counts are available but not surfaced as an acceptance rate.

4. ~~**Default to exponential link in examples**~~ -- **DONE**. README updated.

### New Recommendations

5. **Remove dead C++ code in `phylodiv_tree.cpp`**. All 4 functions
   (`simulate_single_pd_tree_cpp`, `simulate_single_ep_tree_cpp`,
   `simulate_pd_trees_cpp`, `explore_grid_cpp`) have no `Rcpp::export`
   annotations and are never called from R. Also `phylodiv_tree.hpp` can
   be removed. This eliminates compilation time and the only compiler warning.

6. **Remove or deprecate Module 6 (`R/generate.R`)**. All 5 exported functions
   (`generatePhyloPD`, `generateNonHomogeneousExp`, `nhExpRand`, `rate_t`,
   `ExponentialRate`) are unused by the inference pipeline. `generatePhyloPD`
   is superseded by `simulate_tree(model="pd")`. Also removes the C++ file
   `generateNonHomogeneousExp.cpp`. Dropping these would reduce the exported
   API from 27 to 22 functions.

7. **Add tests for inference and GAM modules**. The existing test suite
   (41 pass, 15 skip) covers simulation, augmentation, and CEM/MCEM at a
   basic level but has no tests for:
   - `estimate_rates()` end-to-end (CR/DD with small tree)
   - `compare_models()` / `select_diversification_model()`
   - `diagnose_cem()` / `diagnose_mcem()`
   - GAM functions: `train_GAM()`, `estimate_likelihood_surface()`, etc.
   - IS fhat denominator correction (n_zero_weight > 0)

8. **Exponential link for EP model**. Currently EP + exponential is blocked
   ("no closed-form integral"). Investigate whether numerical integration
   or a different parameterization could enable it. EP is the only model
   restricted to linear link, which is prone to IS collapse.

## Recent Changes (2026-03-07)

- Fixed `emphasis.hpp` copy-paste bug: error message printed `rejected_overruns`
  twice instead of `rejected_zero_weights`.
- Fixed `maxN` default from 10 to 2000.
- Added `.mcem_warn_estep()` diagnostic for zero-weight IS failures.
- Fixed all R CMD check issues: non-ASCII chars, C++14 spec, Rd braces,
  LaTeX build artifacts.
- Fixed MCEM `init_pars` for DD/PD/EP models: slope parameters start at 0
  instead of midpoint of bounds.
- Added proactive warning for PD/EP models with linear link.
- Fixed `select_diversification_model` best-model bug: always returned DD
  due to case mismatch in label lookup.
- Added Module 5 use-case-2 functions: `estimate_likelihood_surface()`,
  `train_likelihood_GAM()`, `find_MLE()`.
- Removed dead code: `R/augment.R`, `R/utils.R` (all 13 functions),
  `get_required_sampling_size()`. Dropped `MASS` from Imports.
- **Fixed IS fhat denominator for zero-weight trees**: C++ `E_step.cpp` and
  R `.is_fhat()` now use `S_completed = N_valid + rejected_zero_weights` as
  denominator, not `N_valid` alone. Overrun/lambda rejections still excluded.
  This removes the upward bias of `log(S_completed / N_valid)`.
- Updated technical report (`doc/emphasis_technical.tex`) Section 3: added
  paragraph distinguishing Type 1 (computational) vs Type 2 (zero-weight)
  augmentation failures and their correct IS handling.
- Updated README examples to use `link = "exponential"` for DD/PD/EP models.

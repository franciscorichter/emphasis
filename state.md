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
| `estimate_rates()` | Yes | Main estimation entry point (mcem/cem/gam) |
| `print.emphasis_fit()` | Yes (S3) | Print method |
| `prune_to_extant()` | No | Drop extinct tips |
| `estimate_rates_control()` | No | Default tuning parameters |
| `compare_models()` | No | AIC comparison table |
| `select_diversification_model()` | No | 3-model comparison pipeline |
| `print.model_selection()` | Yes (S3) | Print method |
| `.resolve_control_aliases()` | No | Sync old/new param names |
| `.extract_brts()` | No | Get branching times from tree/list/vector |
| `.par_names()` | No | `c("beta_0", "beta_N", ...)` |
| `.contract_pars()` | No | 8-element -> compact |
| `.run_mcem()` | No | Dispatch to `.mcem_dynamic_fresh` |
| `.run_cem()` | No | Dispatch to `emphasis_cem` |
| `.run_gam()` | No | Dispatch to GAM pipeline |
| `.model_label()` | No | "CR", "N", "PD" label |

Status: Clean. Three methods: mcem, cem, gam. Only `estimate_rates` exported.
Safe init_pars for DD/PD slopes. Proactive PD/EP warning with linear link.

### Module 3: MCEM Engine (`R/emphasis.R`)
| Function | Exported | Purpose |
|---|---|---|
| `.mcem_dynamic_fresh()` | No | MCEM EM loop with convergence |
| `.mcem_warn_estep()` | No | Diagnose E-step failures |

Status: Clean. Stores rejection counts in `final_IS` for diagnostics.

### Module 4: CEM + IS helpers (`R/de.R`, `R/diagnostics.R`)
| Function | Exported | Purpose |
|---|---|---|
| `diagnose_cem()` | Yes | CEM convergence/IS/acceptance diagnostics |
| `diagnose_mcem()` | Yes | MCEM convergence/IS/acceptance diagnostics |
| `diagnose_gam()` | Yes | GAM surface/fit/MLE diagnostics |
| `print.cem_diagnostics()` | Yes (S3) | Print method with acceptance rate |
| `print.mcem_diagnostics()` | Yes (S3) | Print method with acceptance rate |
| `print.gam_diagnostics()` | Yes (S3) | Print method |
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

Status: Clean. Acceptance rate reported in CEM and MCEM diagnostics.

### Module 5: GAM-based methods (`R/gam.R`)
| Function | Exported | Purpose |
|---|---|---|
| `train_GAM()` | Yes | Train survival-probability GAM (use case 1) |
| `predict_survival()` | Yes | Predict from trained GAM |
| `estimate_likelihood_surface()` | No | IS log-likelihood at parameter grid (use case 2) |
| `train_likelihood_GAM()` | No | Fit GAM to log-likelihood surface |
| `find_MLE()` | No | Optimize GAM to find MLE |

Status: Both use cases implemented. Use case 2 integrated into
`estimate_rates(method="gam")`. GAM basis dimension `k` capped to avoid
errors with small grids.

### Removed files
- `R/augment.R` -- was empty (comment only)
- `R/utils.R` -- all 13 functions were dead code (none called outside utils.R)
- `get_required_sampling_size()` from `R/emphasis.R` -- dead code, allowed dropping `MASS`
- `R/generate.R` -- all 5 exported functions (`generatePhyloPD`,
  `generateNonHomogeneousExp`, `nhExpRand`, `rate_t`, `ExponentialRate`)
  were unused by the inference pipeline. `generatePhyloPD` superseded by
  `simulate_tree(model="pd")`.
- `src/generateNonHomogeneousExp.cpp` -- only supported `R/generate.R`
- `src/phylodiv_tree.cpp` + `inst/include/phylodiv_tree.hpp` -- dead C++
  code (4 functions with no `Rcpp::export`, never called from R). Removed
  to eliminate compilation time.

## C++ Functions (src/)

| C++ function | R name | File | Purpose |
|---|---|---|---|
| `simulate_div_tree_cpp()` | `simulate_div_tree_cpp` | div_tree.cpp | Forward tree simulation |
| `eval_logf_cpp()` | `eval_logf` | loglik.cpp | Evaluate logf and logg |
| `rcpp_mce()` | `augment_trees` | rcpp_mce.cpp | E-step: augment trees |
| `rcpp_mcem()` | `em_cpp` | rcpp_mcem.cpp | Full E+M step |
| `rcpp_mcm()` | `m_cpp` | rcpp_mcm.cpp | M-step only |

All C++ wrappers are internal (accessible via `emphasis:::` only).

## Exported API

7 regular exports + 5 S3 methods:

| Function | Module |
|---|---|
| `estimate_rates()` | Inference |
| `simulate_tree()` | Simulation |
| `train_GAM()` | GAM |
| `predict_survival()` | GAM |
| `diagnose_cem()` | Diagnostics |
| `diagnose_mcem()` | Diagnostics |
| `diagnose_gam()` | Diagnostics |
| `print.emphasis_fit()` | Inference (S3) |
| `print.model_selection()` | Inference (S3) |
| `print.cem_diagnostics()` | Diagnostics (S3) |
| `print.mcem_diagnostics()` | Diagnostics (S3) |
| `print.gam_diagnostics()` | Diagnostics (S3) |

## Test Suite

58 pass, 18 skip. Test files:

| File | Tests | Coverage |
|---|---|---|
| `test-simulate.R` | Simulation helpers, model resolution | Unit + C++ integration (skip) |
| `test-augment.R` | Augmentation / IS | C++ integration (skip) |
| `test-em.R` | MCEM, estimate_rates edges, EP+exp error | Unit + C++ integration (skip) |
| `test-de.R` | CEM, IS fhat, particle init, perturbation | Unit + C++ integration (skip) |
| `test-inference.R` | Inference helpers, compare_models, IS fhat correction | Unit + C++ integration (skip) |
| `test-gam.R` | train_GAM, predict_survival, train_likelihood_GAM, find_MLE, diagnose_gam | Unit (mgcv required) |

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

**Handling**: Zero-weight trees are included in the IS denominator.
`S_completed = N_valid + rejected_zero_weights` is used in both C++
(`E_step.cpp`) and R (`.is_fhat()`). They contribute w=0 to the numerator
and 1 to the denominator.

## Recommendations -- Status

1. ~~**Track N_total alongside N_valid**~~ -- **DONE**.
2. ~~**Correct the fhat denominator**~~ -- **DONE**.
3. ~~**Report acceptance rate in diagnostics**~~ -- **DONE**. Both
   `diagnose_cem()` and `diagnose_mcem()` now compute and display
   acceptance rate (N_valid / S_completed) in their print methods.
   MCEM `final_IS` stores `n_rejected` and `rejected_zero_weights`.
4. ~~**Default to exponential link in examples**~~ -- **DONE**.
5. ~~**Remove dead C++ code (`phylodiv_tree.cpp`)**~~ -- **DONE**. Removed
   `src/phylodiv_tree.cpp` and `inst/include/phylodiv_tree.hpp`. The
   simulation pipeline uses `div_tree.hpp` (separate implementation).
6. ~~**Remove Module 6 (`R/generate.R`)**~~ -- **DONE**. Removed
   `R/generate.R` (5 exported functions), `src/generateNonHomogeneousExp.cpp`,
   `tests/testthat/test-generate.R`, and 5 Rd files. Exported API reduced
   from 12 to 7 functions.
7. ~~**Add tests for inference and GAM modules**~~ -- **DONE**. Added
   `test-inference.R` (14 tests: helpers, compare_models, IS fhat correction,
   estimate_rates validation + integration skips) and `test-gam.R` (7 tests:
   train_GAM, predict_survival, train_likelihood_GAM, find_MLE,
   diagnose_gam). Total: 58 pass, 18 skip.
8. **Exponential link for EP model** -- **OPEN**. EP + exponential is blocked
   in C++ (`model.hpp:77-79`) because the pendant-edge integral has no
   closed-form solution with exponential link. Requires numerical integration
   in C++ (e.g. Gauss-Legendre quadrature) or a different parameterization.
   EP is the only model restricted to linear link, which is prone to IS
   collapse. This is a significant C++ change and deferred.

## Recent Changes (2026-03-07)

### Session 1
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
  denominator, not `N_valid` alone.
- Updated technical report (`doc/emphasis_technical.tex`) Section 3.
- Updated README examples to use `link = "exponential"` for DD/PD/EP models.

### Session 2
- Made C++ wrappers internal: `augment_trees`, `em_cpp`, `m_cpp`, `eval_logf`.
- Made inference helpers internal: `prune_to_extant`, `estimate_rates_control`,
  `compare_models`, `select_diversification_model`.
- Added `method = "gam"` to `estimate_rates()` with `.run_gam()` dispatch.
- Added `diagnose_gam()` and `print.gam_diagnostics()` to diagnostics module.
- Made GAM pipeline helpers internal: `estimate_likelihood_surface`,
  `train_likelihood_GAM`, `find_MLE`.
- Fixed GAM basis dimension `k` capping in `train_likelihood_GAM()` to avoid
  mgcv errors when grid has fewer unique values than default `k`.
- Removed dead C++ code: `src/phylodiv_tree.cpp`, `inst/include/phylodiv_tree.hpp`.
- Removed Module 6: `R/generate.R`, `src/generateNonHomogeneousExp.cpp`,
  `tests/testthat/test-generate.R`, 5 Rd files.
- Added acceptance rate to CEM and MCEM diagnostic print methods.
- Added `n_rejected` and `rejected_zero_weights` to MCEM `final_IS`.
- Added test files: `test-inference.R` (14 tests), `test-gam.R` (7 tests).
- Test suite: 58 pass, 18 skip (up from 41 pass, 15 skip).
- Exported API: 7 functions + 5 S3 methods (down from 19 + 5).

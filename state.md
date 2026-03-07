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

## C++ Exports (src/)

| C++ function | R name | File | Purpose |
|---|---|---|---|
| `simulate_div_tree_cpp()` | `simulate_div_tree_cpp` | div_tree.cpp | Forward tree simulation |
| `generateNonHomogeneousExpCpp()` | `generateNonHomogeneousExpCpp` | generateNonHomogeneousExp.cpp | Non-homogeneous exponential |
| `eval_logf_cpp()` | `eval_logf` | loglik.cpp | Evaluate logf and logg |
| `rcpp_mce()` | `augment_trees` | rcpp_mce.cpp | E-step: augment trees |
| `rcpp_mcem()` | `em_cpp` | rcpp_mcem.cpp | Full E+M step |
| `rcpp_mcm()` | `m_cpp` | rcpp_mcm.cpp | M-step only |

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

**Current handling**: Discarded. IS estimator uses `S_valid` in denominator.

**Correctness**: Justified (see tech report Eq. 3.10). Since the rejection
probability is independent of theta, including these as zero-weight samples
would introduce bias proportional to the user-chosen limit, not the model.
These trees would also have near-zero weight under f (if the tree needed
that many lineages, it's extremely unlikely under the true model).

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

**Current handling**: Discarded from both numerator and denominator.
The fhat computation divides by N (= `sample_size`, the number of VALID
trees requested), not by N_total.

**Correctness issue**: This is subtly biased upward. In correct IS:

    fhat = log( (1/N_total) * sum_{all i} w_i )

Trees with w_i = 0 contribute 0 to the sum but 1 to the denominator. By
excluding them, we compute:

    fhat = log( (1/N_valid) * sum_{valid i} w_i )

This inflates fhat because N_valid < N_total. The bias is:

    bias = log(N_total / N_valid) = -log(acceptance_rate)

When acceptance_rate is high (most trees valid), bias is negligible.
When acceptance_rate is low (IS collapse), bias can be large.

### Recommendations

1. **Track N_total alongside N_valid** in the E-step. Currently the C++ code
   runs exactly `maxN` attempts and collects `N` valid trees. It should also
   return `N_attempted` (= maxN or the number actually tried before stopping).

2. **Correct the fhat denominator** for zero-weight trees: use
   `N_valid + rejected_zero_weights` instead of `N_valid` alone. Overruns and
   lambda rejections stay excluded (they are computational, not IS failures).

3. **Report acceptance rate** in diagnostics: `N_valid / (N_valid + rejected_zero_weights)`.
   Low acceptance rate (< 0.5) is a warning sign for IS quality.

4. **Default to exponential link** for DD/PD/EP models in examples and docs,
   since linear link is prone to IS collapse for these models.

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

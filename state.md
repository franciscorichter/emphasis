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

Status: Clean. All well-structured.

### Module 3: MCEM Engine (`R/emphasis.R`)
| Function | Exported | Purpose |
|---|---|---|
| `.mcem_dynamic_fresh()` | No | MCEM EM loop with convergence |
| `.mcem_warn_estep()` | No | Diagnose E-step failures |
| `get_required_sampling_size()` | No | Legacy; unused anywhere in package |

Status: `get_required_sampling_size()` is dead code (uses `MASS::rlm`). Candidate for removal.

### Module 4: CEM + IS helpers (`R/diagnostics.R`)
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
| `train_GAM()` | Yes | Train survival-probability GAM |
| `predict_survival()` | Yes | Predict from trained GAM |

Status: Under development. Two planned use cases:
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

### Module 7: Legacy utilities (`R/utils.R`)
| Function | Exported | Purpose |
|---|---|---|
| `get_extant()` | No | Get extant lineages at time |
| `transf()` | No | Species name -> index |
| `newick()` | No | Generate Newick string |
| `etree2phylo()` | No | etree -> phylo conversion |
| `phylo2etree()` | No | phylo -> etree conversion |
| `GPD()` | No | Pairwise phylogenetic diversity matrix |
| `n_from_time()` | No | Count lineages at time |
| `n_for_all_bt()` | No | Count lineages at all branching times |
| `extend_tree()` | No | Extend tree with events |
| `phylodiversity()` | No | Compute total PD |
| `get_current_species()` | No | Species alive at time |
| `get.time()` | No | Elapsed time |
| `L2phylo()` | No | L-table -> phylo (from DDD) |

Status: Mixed. Some used by simulate/augment (`L2phylo`, `etree2phylo`, `phylo2etree`).
Others need dead-code audit.

### Module 8: Empty file (`R/augment.R`)

Status: Remove. Contains only a comment.

## C++ Exports (src/)

| C++ function | R name | File | Purpose |
|---|---|---|---|
| `simulate_div_tree_cpp()` | `simulate_div_tree_cpp` | div_tree.cpp | Forward tree simulation |
| `generateNonHomogeneousExpCpp()` | `generateNonHomogeneousExpCpp` | generateNonHomogeneousExp.cpp | Non-homogeneous exponential |
| `eval_logf_cpp()` | `eval_logf` | loglik.cpp | Evaluate logf and logg |
| `rcpp_mce()` | `augment_trees` | rcpp_mce.cpp | E-step: augment trees |
| `rcpp_mcem()` | `em_cpp` | rcpp_mcem.cpp | Full E+M step |
| `rcpp_mcm()` | `m_cpp` | rcpp_mcm.cpp | M-step only |

## Cleanup Priorities

| Priority | Action | Impact |
|---|---|---|
| High | Remove `R/augment.R` (empty) | Cleanliness |
| High | Remove `get_required_sampling_size()` (unused) | Dead code; may drop `MASS` import |
| Medium | Consider removing `R/generate.R` (superseded) | Reduces API surface |
| Low | Audit `R/utils.R` for dead functions | Cleanliness |

## Recent Bug Fixes (2026-03-07)

- Fixed `emphasis.hpp` copy-paste bug: error message printed `rejected_overruns`
  twice instead of `rejected_zero_weights`, hiding the real failure cause.
- Fixed `maxN` default from 10 to 2000: the old value was less than `sample_size`,
  making E-step always fail.
- Added `.mcem_warn_estep()` diagnostic for zero-weight IS failures.
- Fixed all R CMD check issues: non-ASCII chars, C++14 spec, Rd braces,
  LaTeX build artifacts.
- Fixed MCEM `init_pars` for DD/PD/EP models: slope parameters now start at 0
  instead of midpoint of bounds, preventing IS collapse where `lambda = max(0,
  beta_0 + beta_X * X)` hits zero for large N or PD.
- Added proactive warning for PD/EP models with linear link and negative bounds.
- Fixed `select_diversification_model` best-model bug: `label_to_key` used
  uppercase keys (`N`, `PD`, `EP`) but `compare_models` passed lowercase names
  (`dd`, `pd`, `ep`), so `best_key` was always `NA` and fell back to `dd`
  regardless of AIC.
- Added Module 5 use-case-2 functions: `estimate_likelihood_surface()`,
  `train_likelihood_GAM()`, `find_MLE()`.

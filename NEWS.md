# emphasis 1.0.0

First release on CRAN. The package was rewritten after version 0.4
(2021--2024). The functions of that version (`emphasis()`, `emphasis_de()`,
`sim_tree_pd_cpp()`, `augmentPD()`, `loglikelihood()` and the rest) are gone;
the interface is the one described below.

## Changes since version 0.4

* The default link of `estimate_rates()`, `simulate_tree()`, `auto_bounds()`
  and `emphasis_pipeline()` is now `"exponential"`, under which the
  intercepts are log rates; it was `"linear"`. The model ladder, the
  real-tree analyses, the exact comparison on real trees and the conditioned
  estimability study use the exponential link; the earlier simulation studies
  ran on the linear link and remain valid for it. Scripts relying on the linear default
  must pass `link = "linear"`, which stays supported. The `"gaussian"` link
  is experimental: it has been used in one convergence study of the surrogate
  proposal and in no estimation experiment, and it is not available for `ED`
  or `K` models.
* A tree that is not ultrametric (one with branching times after the present)
  is refused with a message, instead of being scored with the present read as
  an observed split.
* The `D` likelihood on a bare branching-time vector counted one lineage too
  many after an observed split, so a `d` model with zero `D` coefficients did
  not reduce to `cr`; the alive lineages now follow the forward sweep by
  lineage id. Trees given as `phylo` objects were not affected.
* One simulator, `simulate_tree()`, replaces the separate forward simulators
  and augmentation functions. It simulates trees forward, in batches over a
  parameter matrix, or conditionally on an observed tree (augmentation).
* One estimator, `estimate_rates()`, replaces `emphasis()` and
  `emphasis_de()`, with the method as an argument.
* Incomplete taxon sampling: the sampling fraction `rho` enters the
  simulator, the likelihood and both augmentation samplers.
* Conditioning on survival is the `cond` argument of `estimate_rates()`: the
  likelihood is divided by the survival probability of the clade, read from
  an emulator. `train_GAM()` now fits that emulator, a binomial GAM of survival
  on the parameters; in version 0.4 it fitted a GAM to simulated
  log-likelihoods.
* Lineage-level covariates: the age-imbalance `D`, evolutionary
  distinctiveness `ED` (raw or centred) and the number of ancestral splits `K`
  join diversity `N`. Version 0.4 read diversity and the clade's phylogenetic
  diversity.
* A second augmentation proposal: the conditional birth-death-immigration
  sampler, exact for constant rates and a mean-field proposal for covariate
  models, is the default; the thinning sampler of version 0.4 remains.
* Model comparison and selection: `compare_models()`, `lineage_table()` and
  `covariate_path()`.
* Computation runs on one thread unless `num_threads` is set.

## Models

* Speciation and extinction rates depend on four lineage-level covariates read
  on the complete tree: diversity `N`, age-imbalance `D`, evolutionary
  distinctiveness `ED` (raw or centred, `EDc`) and the number of ancestral
  splits `K`. Models are named (`"cr"`, `"dd"`, `"d"`, `"nd"`, `"ed"`,
  `"ned"`, `"edc"`, `"nedc"`, `"k"`, `"nk"`) or given as a formula such as
  `~ N + ED`, on an exponential (the default), linear or gaussian
  (experimental) link.

## Estimation

* `estimate_rates()` fits a model by Monte Carlo expectation-maximization
  (`"mcem"`), the cross-entropy method (`"cem"`) or a GAM fitted to a
  likelihood surface (`"gam"`). `model = "all"` fits the nested ladder of
  models and ranks it.
* `auto_bounds()` builds a parameter box from the observed tree by forward
  simulation.
* `compare_models()` compares fits by AIC and reports the Monte Carlo error of
  each log-likelihood and the downward bias of a log of a Monte Carlo mean.
* `lineage_table()` returns the per-lineage design the likelihood is built on,
  and `covariate_path()` traces a selection path over the covariates.
* The survival emulator is a GAM from `train_GAM()`, an `nnet` object or a
  function, evaluated by `predict_survival()`.
* `diagnose_mcem()`, `diagnose_cem()` and `diagnose_gam()` report convergence
  and importance-sampling diagnostics.

## Behaviour

* Progress output of the fitting functions goes through `message()` and can be
  silenced with `suppressMessages()` or `verbose = FALSE`.
* `emphasis_pipeline()` runs on one thread unless `control$num_threads` is
  set.
* `simulate_tree()`'s `max_tries` is the number of forward draws in all, the
  first included: the default `1` is a single draw, unconditioned on
  survival, and an extinct clade is redrawn until one survives or `max_tries`
  draws have been made. `max_tries = 0`, which asked for a single draw when
  the argument counted the draws after the first, is read as `1` with a
  warning. `auto_bounds()` passes `surv_tries` (default 20) as `max_tries`,
  and the survival GAM is trained on single draws.
* `control$max_time` binds inside the E-step of an MCEM fit. Both samplers
  read the clock between draws and inside a draw, stop drawing when it runs
  out, and the fit ends with `stop_reason = "time_budget"` and a `clock` note;
  the iteration in progress is discarded. The E-step that evaluates the
  returned parameters gets `max_time` again and, when it runs out, estimates
  the log-likelihood from the draws it made, so a fit returns within about
  twice `max_time` plus one M-step. The M-step is not bounded by the clock.
* A fit carries `notes`, the named notes `estimate_rates()` states about it.
  When the surrogate proposal's mean-field iteration was relaxed below the
  full step, ran under a `damping` other than `"auto"`, or did not converge
  within its sweep budget, the `mean_field` note gives the setting, the sweeps
  taken, the relaxation and the final residual at the last E-step and whether
  it converged, and `print()` shows it. Every E-step's mean field is in the
  trace (`mf_converged`, `mf_delta`, `mf_sweeps`, `mf_omega`, `mf_budget`).
* `compare_models()` chooses the model with the fewest parameters within
  `margin` of the minimum AIC and, among those, the one with the lower AIC,
  whatever order the fits are passed in.

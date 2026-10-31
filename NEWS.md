# emphasis 1.1.0

* The trees an E-step returns under the thinning proposal carry their lineage
  ids and covariate columns again. The copy `rcpp_mcem.cpp` made of them had
  only `brts`, `n` and `t_ext`, so under `stop_rule = "mc_error"` every batch
  M-step of a `D`, `ED` or `K` model failed to build its covariate table, the
  Monte Carlo error of the step was never computed (`mc_z` was `NA` on every
  iteration), the draws never grew and no thinning fit of those models met
  the rule: it ran to `max_iter` whatever the estimate did. The surrogate
  proposal was not affected.

* `estimate_rates(cond = "surface")` conditions a fit on survival for any
  model and link: an unconditioned fit, then MCEM continued on local
  quadratic surfaces of a direct Monte Carlo estimate of log P(S) around the
  iterate, rebuilt whenever the iterate leaves the inner half of a surface's
  box, with a second run from the nested start when the unconditioned
  extinction intercept is on an edge of the box. The reported log-likelihood
  subtracts a direct estimate of log P(S) at the estimate, and the fit's
  `conditioning` records the surfaces, clades, seconds and the start kept.
  It is the recommended conditioning for the lineage-level models, whose
  survival GAM is not accurate enough (E43). The settings are `control`
  entries named `surface_*`.
* `survival_mc()` estimates log P(S) by direct simulation for any model, link
  and sampling fraction, with common random numbers across parameter values;
  `survival_surface()` fits the local surface, with `predict()` and `print()`
  methods. `predict_survival()` and `estimate_rates(cond = )` accept a
  surface.

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
  must pass `link = "linear"`, which stays supported.
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
  likelihood is divided by the survival probability of the clade, given as a
  function of the parameters, a GAM or an `nnet` object. `train_GAM()` now
  fits a binomial GAM of survival on the parameters; in version 0.4 it fitted
  a GAM to simulated log-likelihoods. That GAM is not accurate enough to
  condition on: in E43 (https://franciscorichter.github.io/emphasis/experiments/E43.html) its
  90th-percentile error of log P(S) was 0.44 to 2.56 for every lineage-level
  model at 2000, 5000 and 20 000 training clades, against 0.05 required.
  Under `cr` and `dd` condition on the exact survival probability, passed as
  a function (`?estimate_rates` gives the `cr` formula); under lineage-level
  models on a direct Monte Carlo estimate, or not at all.
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
  `~ N + ED`, on an exponential (the default) or linear link.

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
* The survival probability `cond` takes is a function, a GAM from
  `train_GAM()` or an `nnet` object, evaluated by `predict_survival()`.
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
* `control$max_time` binds inside the E-step and the M-step of an MCEM fit,
  and the fit ends with `stop_reason = "time_budget"` and a `clock` note
  naming the step it cut. Both samplers read the clock between draws and
  inside a draw and stop drawing when it runs out; the iteration in progress
  is discarded. The M-step's optimiser (NLopt's sbplx, under both samplers)
  reads the same deadline before every evaluation of its objective and stops
  once it has passed; the iteration then keeps the best point the optimiser
  had evaluated when that point improves on the current iterate, and is
  discarded otherwise. The E-step that evaluates the returned parameters gets
  `max_time` again and, when it runs out, estimates the log-likelihood from
  the draws it made, so a fit returns within about twice `max_time`.
* A fit carries `notes`, the named notes `estimate_rates()` states about it.
  When the surrogate proposal's mean-field iteration was relaxed below the
  full step, ran under a `damping` other than `"auto"`, or did not converge
  within its sweep budget, the `mean_field` note gives the setting, the sweeps
  taken, the relaxation and the final residual at the last E-step and whether
  it converged, and `print()` shows it. Every E-step's mean field is in the
  trace (`mf_converged`, `mf_delta`, `mf_sweeps`, `mf_omega`, `mf_budget`).
* `control$damping` defaults to `"default"`, which is `"auto"` at complete
  sampling (`rho = 1`) and a fixed relaxation of 1/8 when `rho < 1`. In
  experiment E10c a fixed 1/8 converged 0 of 192 units within the 20 sweeps
  at `rho = 1` against 131 for `"auto"`, at 1.4--3.3 times the E-step time
  and no gain in effective sample, and below `rho = 1` it converged more
  (142 against 131 at `rho = 0.2`). An explicit `"auto"` or number is used as
  given, and the `mean_field` note names the resolved setting and whether it
  is the default.
* `compare_models()` chooses the model with the fewest parameters within
  `margin` of the minimum AIC and, among those, the one with the lower AIC,
  whatever order the fits are passed in.

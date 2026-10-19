# Usage: one session, from a simulated tree to a chosen model

This page runs the package once, on one tree, through every route the [overview](index.html) describes: a tree is simulated under a diversity-dependent model, the ladder of nested models is fitted to its reconstructed tree, the fits are compared, the covariate path is traced on the observed tree, the constant-rate fit is conditioned on survival through its exact probability and through the survival emulator, the Monte Carlo EM run is inspected, and the earlier four-stage pipeline is run on the same tree. Every chunk was executed in the order shown, in one R session, and the block under each chunk is its console output as the session wrote it. A sentence after each block says what to read in it.

The budgets are small so that the whole session runs in minutes: 100 augmented trees per E-step and 25 EM iterations where the package defaults are 200 and 200. The estimates therefore stop at the iteration cap rather than at the stopping rule, and their Monte Carlo digits are those of a 100-draw estimate. The fits are not the measurements the [experiments](experiments/index.html) report, which use the defaults or more; they show the calls, the objects they return and how to read them.

## 1. The session

The session ran on 19 October 2026 under R 4.3.3 on x86_64 Linux with `emphasis` 0.4 at commit `4e0b9dd`, the commit that holds $D$ constant between events (the later commits on the branch change the README and a test, not the package), with `ape` 5.8.1 and `mgcv` 1.9.1. <!-- out/session.txt of the run; the commit from the reflog of the worktree: 4e0b9dd is followed by e1c2e82 (tests) and four README commits --> The package is installed with

```r
remotes::install_github("franciscorichter/emphasis")
```

Two things about reproducibility. `set.seed()` fixes the forward simulator, so the tree of Section 2 is the same on every machine. The fits use `num_threads = 4`, and augmentation on more than one thread is not reproducible: each worker draws from its own substream and the sample keeps the first draws to finish, so a re-run gives the same tree and fits that differ in their Monte Carlo digits. <!-- simulate_tree, section Reproducibility --> The elapsed and CPU time of each chunk are in Section 9.

## 2. A simulated tree

`simulate_tree()` runs the forward process from the crown. The generator is the diversity-dependent model `dd` on the exponential link, $\lambda(N) = 0.5\,e^{-0.02 N}$ and $\mu = 0.1$, over a crown age of 10. The parameter vector follows the layout of the overview's [model table](index.html#3-the-models): the $\beta$ coefficients of the active covariates, then the $\gamma$ coefficients.

```r
library(emphasis)
library(ape)

set.seed(6)
theta <- c(beta_0 = log(0.5), beta_N = -0.02, gamma_0 = log(0.1), gamma_N = 0)
sim <- simulate_tree(pars = theta, max_t = 10, model = "dd", link = "exponential")
sim$status
phy <- sim$tes
c(extant_tips = Ntip(phy), tips_with_extinct = Ntip(sim$tas), crown_age = max(branching.times(phy)))
```

```output
[1] "done"
      extant_tips tips_with_extinct         crown_age 
               59                83                10 
```

`status` is `"done"`: the clade reached the present at the first attempt. The reconstructed tree `tes` has 59 tips and the complete tree `tas` 83, so 24 lineages arose and went extinct before the present and are what the fits below have to integrate over. At 59 lineages the generating speciation rate is $0.5\,e^{-1.18} = 0.15$, against an extinction rate of 0.10. `phy` is the only object the rest of the session reads; `theta` is kept to compare the estimates with.

## 3. The ladder

`estimate_rates(model = "all")` fits `cr`, then `dd`, then `nd`, `ned` and `nk`, each started from the estimate of the model it nests and in a box that `auto_bounds()` builds for it ([overview, Section 7](index.html#7-monte-carlo-em)). The ladder returns the five fits with their `compare_models()` table.

```r
set.seed(1)
fits <- estimate_rates(phy, model = "all", link = "exponential",
                       control = list(sample_size = 100, max_iter = 25, num_threads = 4))
fits
```

```output
Warning: BDI augmentation: the mean-field iteration did not converge (delta = 9.47e-02 after 20 sweeps); the proposal uses the last sweep.
[13 further warnings of the same form, with deltas from 4.04e-03 to 1.06e+00, collapsed here]
Warning: the last E-step's effective sample is 6.7: the fit is sampler-limited and the precision of its estimate and of its log-likelihood is unquantified. The pre-registered rules exclude fits below 10 from every selection statement; more draws, the other proposal (control$sampling), or the exponential link may help.
sampling = "bdi" is not available for a K-dependent model: the sampler's mean field does not carry the ancestral-split count; the thinning proposal reads K per lineage; using the thinning sampler ("dynamic_fresh") instead.
Warning: compare_models(): the ranking is not safe against Monte Carlo precision. "ned" is behind "dd" by 8.08 AIC, but its effective sample is smaller, which inflates its AIC by up to 0.10 on its own. Raise num_trees on the lower-ESS fits and compare again.
Note: nd and ned both read the pendant edge; on 100-tip trees one path over N, D and ED chooses the pair from a generator with only one of them (0.60 from nd, 0.45 from ned), so an AIC margin between the two is a margin between two readings of the same edge.
emphasis fits (mcem, link = exponential): cr -> dd -> nd -> ned -> nk

  model n_pars    loglik      AIC    loglik_se        ESS  loglik_gap gap_heavy delta_AIC
1    dd      4 -133.9155 275.8311 1.134813e-01  21.048225 0.020330315     FALSE  0.000000
2    nk      6 -132.6659 277.3317 1.407005e-01  53.325920 0.004376303     FALSE  1.500652
3    cr      2 -139.9438 283.8875 1.733163e-14 100.000000 0.000000000     FALSE  8.056445
4   ned      6 -135.9546 283.9091 2.411951e-01   6.741967 0.069888828      TRUE  8.078026
5    nd      6 -136.2393 284.4785 1.760942e-01  14.304060 0.029955110     FALSE  8.647417
         AICw   AIC_swing within_margin
1 0.657512293  0.00000000          TRUE
2 0.310485623 -0.03190802          TRUE
3 0.011707635 -0.04066063         FALSE
4 0.011581981  0.09911703         FALSE
5 0.008712468  0.01924959         FALSE

chosen: dd (the fewest parameters within 2 AIC of the minimum)
```

The table is the `compare_models()` ranking, one row per rung, ordered by AIC. Two models are within 2 AIC of the minimum, `dd` with 4 parameters and `nk` with 6, and `chosen` is `dd`, the generating model. The `cr` row has an effective sample of 100 of 100 and a standard error of $10^{-14}$: under constant rates the surrogate proposal is the exact conditional distribution and every draw has the same weight. The `ned` row has an effective sample of 6.7 and `gap_heavy` is `TRUE`, which is what the two warnings say: a fit below 10 effective draws is excluded from every selection statement, and `ned` is the model `ranking_at_risk` names, since its Jensen gap of 0.070 exceeds the leader's 0.020 and correcting both could move it by up to 0.10 AIC, against a deficit of 8.08. The 14 warnings before the table come from the mean-field sweep of the surrogate proposal, which did not converge in 20 iterations on 14 E-steps; the proposal then used the last sweep. The message about `nk` says that the surrogate does not carry $K$, so that rung ran on the thinning proposal. Every rung stopped at the 25-iteration cap. In the `dd` estimate, $\gamma_N = 0.050847 = 3/59$ is the cap on the $N$ slope in the `auto_bounds()` box of a 59-tip tree ([overview, Section 7](index.html#7-monte-carlo-em)), so the extinction slope sits on the box's upper edge, where its generating value is 0.

## 4. Comparing the fits

A fit prints its estimate, its log-likelihood with the Monte Carlo standard error, the effective sample of the last E-step and the Jensen gap. `compare_models()` ranks any set of fits by AIC and reports `chosen`, the model with the fewest parameters within `margin` AIC of the minimum ([overview, Section 9](index.html#9-model-choice-and-the-three-routes)); it is the function the ladder calls, so calling it on the ladder's fits reproduces the table above.

```r
fits$fits$dd
cmp <- compare_models(cr = fits$fits$cr, dd = fits$fits$dd, nd = fits$fits$nd,
                      ned = fits$fits$ned, nk = fits$fits$nk, margin = 2)
cmp[, c("model", "n_pars", "loglik", "loglik_se", "ESS", "loglik_gap", "delta_AIC", "within_margin")]
attr(cmp, "chosen")
attr(cmp, "ranking_at_risk")
```

```output
emphasis fit (mcem, model = N)

Parameters:
   beta_0    beta_N   gamma_0   gamma_N 
-0.341926 -0.021438 -4.672061  0.050847 

Log-likelihood: -133.9155  (MC se: 0.1135)
  effective sample: 21; the estimate sits at least 0.020 below the truth (a log of a
  Monte Carlo mean is biased down), which moves AIC by 0.041
AIC:            275.8311 
n_pars:         4 
Iterations:     25 
Stop reason:    max_iter 
Warning: compare_models(): the ranking is not safe against Monte Carlo precision. "ned" is behind "dd" by 8.08 AIC, but its effective sample is smaller, which inflates its AIC by up to 0.10 on its own. Raise num_trees on the lower-ESS fits and compare again.
  model n_pars    loglik    loglik_se        ESS  loglik_gap delta_AIC within_margin
1    dd      4 -133.9155 1.134813e-01  21.048225 0.020330315  0.000000          TRUE
2    nk      6 -132.6659 1.407005e-01  53.325920 0.004376303  1.500652          TRUE
3    cr      2 -139.9438 1.733163e-14 100.000000 0.000000000  8.056445         FALSE
4   ned      6 -135.9546 2.411951e-01   6.741967 0.069888828  8.078026         FALSE
5    nd      6 -136.2393 1.760942e-01  14.304060 0.029955110  8.647417         FALSE
[1] "dd"
[1] "ned"
```

The printed fit carries the estimate in the compact layout $(\beta_0, \beta_N, \gamma_0, \gamma_N)$ on the log-rate scale of the exponential link, the log-likelihood $-133.92$ with its bootstrap standard error 0.11, the effective sample 21 and the lower bound 0.020 on the Jensen gap, which moves the AIC by 0.041. The table repeats the ladder's ranking on the columns selected, the first attribute is the chosen model and the second names `ned` as the one fit whose position could change once the gaps are corrected; the warning that accompanies it is the same one the ladder printed. The columns left out here are `AIC`, `AICw`, `AIC_swing` and `gap_heavy`.

## 5. The covariate path

`covariate_path()` traces one differential-geometric LARS path over the covariates a rate can read, on the lineage table of the tree, and picks a point on it by BIC. Given a single `phylo` it treats the tree as fully observed: the speciation path is exact, and the extinction path has no event to read.

```r
cp <- covariate_path(phy, covariates = c("N", "D", "ED", "K"), criterion = "BIC")
cp$speciation$entry_order
cp$speciation$entry_gamma
cp$speciation$active
cp$speciation$coef
c(rows = cp$speciation$n_rows, events = cp$speciation$n_events)
cp$extinction$n_events
```

```output
[1] "K"  "ED" "D"  "N" 
         N          D         ED          K 
0.04970953 0.36378635 0.58108040 1.04345450 
[1] "K"
(Intercept)           N           D          ED           K 
 -0.7485024   0.0000000   0.0000000   0.0000000  -0.1680652 
  rows events 
  1769     57 
[1] 0
```

The covariates enter the speciation path in the order $K$, $\mathrm{ED}$, $D$, $N$, at $\gamma$ = 1.04, 0.58, 0.36 and 0.05, and BIC keeps $K$ alone, with a coefficient of $-0.168$ per ancestral split on the log rate and an intercept of $-0.749$, a rate of $e^{-0.749} = 0.47$ at $K = 0$. $N$, the generating covariate, enters last. The lineage table has 1,769 rows, one per lineage and segment, and 57 events, the $n - 2$ splits a 59-tip tree records after the crown. The extinction path has 0 events, since no extinction is observed on a reconstructed tree. On a single observed tree the path reads the tree as complete, with 24 lineages missing from it, which is why the overview calls it a diagnostic of which covariates the data ask for first and not a route of record ([Section 9](index.html#9-model-choice-and-the-three-routes)); on the augmented trees of an E-step, with their weights, the same call gives the path the M-step sees.

## 6. Conditioning on survival

The likelihood is conditioned on the survival of both crown lineages when `estimate_rates()` is given a survival probability through `cond` ([overview, Section 8](index.html#8-conditioning-on-survival)). Under `cr` that probability is exact: with $d = \lambda - \mu$ and crown age $T$, one lineage survives with probability $p = d/(\lambda - \mu e^{-dT})$ and the clade with probability $p^2$. A function of the parameter data frame is the first of the three forms `cond` accepts, and `predict_survival()` evaluates any of the three.

```r
age <- max(branching.times(phy))
p_cr <- function(pars) {
  lambda <- exp(pars$beta_0); mu <- exp(pars$gamma_0); d <- lambda - mu
  p <- ifelse(abs(d) < 1e-12, 1 / (1 + lambda * age), d / (lambda - mu * exp(-d * age)))
  p^2
}
theta_cr <- as.data.frame(t(fits$fits$cr$pars))
theta_cr
predict_survival(p_cr, theta_cr)
```

```output
     beta_0   gamma_0
1 -1.414223 -4.313994
[1] 0.9029204
```

At the unconditioned `cr` estimate, $\lambda = e^{-1.414} = 0.243$ and $\mu = e^{-4.314} = 0.013$, one crown lineage survives 10 time units with probability 0.950 and the clade with probability $0.950^2 = 0.903$. The function takes a data frame with the parameter names of the fit, which is the form `cond` passes it.

The second form is the binomial GAM that `auto_bounds()` trains over its box, 2,000 simulated clades labelled survived or extinct. The grid compares it with the exact probability at four parameter values.

```r
set.seed(2)
ab_cr <- auto_bounds(phy, model = "cr", link = "exponential", verbose = FALSE)
grid <- data.frame(beta_0 = log(c(0.3, 0.3, 0.5, 0.5)), gamma_0 = log(c(0.05, 0.2, 0.05, 0.2)))
cbind(grid, exact = predict_survival(p_cr, grid), gam = predict_survival(ab_cr$survival_gam, grid))
```

```output
Training univariate survival GAM on 2000 simulations (survival rate: 0.478 )...
      beta_0   gamma_0     exact       gam
1 -1.2039728 -2.995732 0.7138428 0.7209476
2 -1.2039728 -1.609438 0.1950539 0.2688857
3 -0.6931472 -2.995732 0.8118027 0.8247545
4 -0.6931472 -1.609438 0.3747787 0.4011796
```

The GAM was trained on 2,000 clades of which 0.478 survived. At three of the four points it is within 0.07 of the exact probability on the log scale; at the point of lowest survival, $\lambda = 0.3$ and $\mu = 0.2$, it gives 0.269 for an exact 0.195, 0.32 above it on the log scale. The 90th percentile of that error over a whole box, measured on `cr` and `dd` at 100 and 200 tips with 2,000 clades, is 0.28 ([E36](experiments/E36.html)).

The conditioned fit of `cr` starts from the unconditioned estimate in the box `auto_bounds()` found.

```r
fit_cr_cond <- estimate_rates(phy, model = "cr", link = "exponential", cond = p_cr,
                              init_pars = fits$fits$cr$pars,
                              control = list(lower_bound = ab_cr$lower_bound,
                                             upper_bound = ab_cr$upper_bound,
                                             sample_size = 100, max_iter = 25, num_threads = 4))
rbind(unconditioned = exp(fits$fits$cr$pars), conditioned = exp(fit_cr_cond$pars))
c(unconditioned = fits$fits$cr$loglik, conditioned = fit_cr_cond$loglik)
```

```output
                 beta_0     gamma_0
unconditioned 0.2431143 0.013380009
conditioned   0.2370679 0.003973605
unconditioned   conditioned 
    -139.9438     -139.6893 
```

Conditioning moves the speciation rate from 0.243 to 0.237 and the extinction rate from 0.0134 to 0.0040, a factor of 3.4: the extinction side is the parameter the survival term acts on. The two log-likelihoods are values of different functions, the marginal likelihood (6) and the conditioned one (15) of the overview, and are not compared with each other. The box is the one `auto_bounds()` built for `cr` in the previous chunk; `init_pars` starts the run from the unconditioned estimate.

## 7. Diagnostics

`diagnose_mcem()` returns the convergence table of a fit, one row per EM iteration, and the importance-sampling quality of its final E-step; with `plot = TRUE` it also draws the traces.

```r
dg <- diagnose_mcem(fits$fits$dd, plot = FALSE)
dg
tail(dg$convergence[, c("fhat", "delta_max", "num_trees")], 5)
```

```output
MCEM diagnostics
================

Iterations:      25
Stop reason:     max_iter
Final fhat:      -133.9155
Final delta_max: 4.88e-02

IS quality (final E-step):
  num_trees:    100
  acceptance:   68%
  ESS:          21.0  (21% of num_trees)
  mean lw:      -134.0489
  sd lw:        0.9147
  fhat:         -133.9155
  loglik_se:    0.1135  (bootstrap, B=200)
        fhat  delta_max num_trees
21 -134.5362 0.06638218       100
22 -134.4519 0.02269657       100
23 -134.2903 0.01991885       100
24 -134.1171 0.05379311       100
25 -133.9953 0.04877226       100
```

The run stopped at `max_iter` after 25 iterations with a final step $\delta_{25} = 0.049$, above the $10^{-2}$ that the default rule needs on three consecutive iterations, so the run had not converged under its own rule; over the last five iterations $\hat\ell$ still rises from $-134.54$ to $-134.00$, by about one standard error (0.11) per iteration. The last E-step reports an acceptance of 68% and an effective sample of 21 of 100, with a standard deviation of 0.91 in the log weights. `Final fhat` ($-133.9155$) differs from row 25 of the trace ($-133.9953$): the trace row is the E-step at the iterate the M-step started from, and the final value is the extra E-step at the returned estimate, which is what the fit reports as its log-likelihood.

## 8. The pipeline

`emphasis_pipeline()` is the earlier workflow kept for compatibility ([overview, Section 11](index.html#11-the-simulator-the-diagnostics-the-pipeline)): `auto_bounds()` with its survival GAM, a likelihood surface fitted on a grid and maximised, a cross-entropy search over the box, and MCEM from the CEM result, conditioned on the GAM. The grid and the search are cut to a fraction of their defaults here (30 grid points of 30 draws where the defaults are 150 of 200; 16 particles over 3 iterations where the defaults are 50 over 20).

```r
set.seed(3)
pl <- emphasis_pipeline(phy, model = "dd", link = "exponential",
                        control = list(num_threads = 4,
                                       gam = list(n_grid = 30, sample_size = 30),
                                       cem = list(num_particles = 16, max_iter = 3),
                                       mcem = list(sample_size = 100, max_iter = 25)),
                        verbose = FALSE)
pl
rbind(generating = exp(theta), pipeline = exp(pl$pars), ladder_dd = exp(fits$fits$dd$pars))
```

```output
Training univariate survival GAM on 2000 simulations (survival rate: 0.427 )...
Training univariate likelihood GAM on 24 grid points...
emphasis pipeline (model = N, conditioned)
================================================== 

Stage log:
    bounds      ok  loglik=       ---  AIC=       ---   34.0s
    gam         ok  loglik=   1663.75  AIC=  -3319.50    2.4s
    cem         ok  loglik=   -136.57  AIC=    281.15    4.8s
   *mcem        ok  loglik=   -135.53  AIC=    279.06   11.3s

Best stage: mcem
Parameters:
   beta_0    beta_N   gamma_0   gamma_N 
-0.599915 -0.022359 -5.522949  0.050847 

Log-likelihood: -135.5295
AIC:            279.0590
Total time:     52.5s
              beta_0    beta_N     gamma_0  gamma_N
generating 0.5000000 0.9801987 0.100000000 1.000000
pipeline   0.5488583 0.9778891 0.003994050 1.052162
ladder_dd  0.7104009 0.9787901 0.009352975 1.052162
```

The stage log gives each stage's time: 34.0 s for the bounds and the survival GAM, 2.4 s for the surface on the 24 of the 30 grid points at which the likelihood estimate was finite, 4.8 s for the cross-entropy search and 11.3 s for the 25 MCEM iterations. The GAM stage's value of 1663.75 is the maximum of a smoothed surface extrapolated from 24 points of 30 draws and is not a log-likelihood of the tree; the pipeline starts MCEM from the CEM result when one exists and reports the MCEM stage, marked with the asterisk, as its result. Its log-likelihood $-135.53$ is conditioned on the survival GAM and is not compared with the ladder's unconditioned $-133.92$. The last table is on the rate scale: the intercepts are rates and the slopes are factors per lineage. The pipeline's speciation intercept is 0.549 against a generating 0.5 and the ladder's 0.710; the two $N$ factors, 0.978 and 0.979, sit just below the generating 0.980; the extinction intercept is 0.0040 against a generating 0.1, and the extinction slope sits at the box edge in both fits. Extinction is the side a reconstructed tree carries least, and at this size and turnover the measurements say it is not read ([overview, Section 10](index.html#10-what-is-established)).

## 9. What the session cost

Elapsed and CPU seconds of each chunk on four threads of an AMD Ryzen 9 9955HX, as `proc.time()` measured them around each chunk. <!-- out/times.txt of the run -->

| Chunk | Elapsed s | CPU s |
|---|---|---|
| `01-simulate` | 0.1 | 0.1 |
| `02-ladder` | 445.7 | 663.8 |
| `03-compare` | 0.0 | 0.0 |
| `04-path` | 0.0 | 0.0 |
| `05-survival` | 39.7 | 47.6 |
| `06-diagnose` | 0.0 | 0.0 |
| `07-pipeline` | 52.5 | 150.6 |

The session took 538 s elapsed and 862 s of CPU, 83% of the elapsed time in the ladder, whose five rungs each ran 25 E-steps of 100 draws and a last E-step at the estimate. The ratio of CPU to elapsed time is 1.5 on the ladder and 2.9 on the pipeline against the 4 threads asked for: the parts of the run that are R code use one thread: the surrogate proposal's draws and their mean-field sweep (`R/bdi.R`), and the `auto_bounds()` calls of the ladder, which does not pass `num_threads` on to them; the thinning proposal, the likelihood and the M-step use four, as do the forward simulations of `auto_bounds()` inside the pipeline. At the package defaults, 200 draws and up to 200 iterations, the same ladder costs up to 16 times the E-step work shown here; the ladder on 170 real trees of 50 to 200 tips ran with 300 draws and 300 iterations under a four-hour clock per fit ([Real phylogenies](real-phylogenies.html)).

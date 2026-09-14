# The convergence criterion, and why the stopping rule was not replaced

**Status: measured, not shipped.** Branch `worktree-wf_1b9d56bb-780-5` carries a working
implementation; it is not merged, and the reason is below.

## What prompted it

Convergence is the most-failed criterion in the validation study: of 124 cells only 51 converge in
at least 95 % of their fits, and `i:convergence` is the failing criterion in 68 cells, more than
any other. The rule in force compares a parameter step to the parameter's own scale
(`.rel_change` / `.rel_floor`, R/inference.R) and never to the Monte Carlo noise of the iterate it
is testing. At `num_trees = 200` a quiet draw near the fixed point has probability about 0.5, so
three quiet iterations is a geometric wait rather than evidence.

The replacement follows the MCEM literature: measure the per-iteration ascent of the EM surrogate
`Q` on the trees the M-step already consumed, bootstrap its Monte Carlo spread by resampling those
trees, declare convergence when the ascent can no longer be distinguished from that spread, and
grow the sample when it cannot.

## The measurement

Two independent measurements were made on disjoint samples of the cells that fail the convergence
criterion, each against exact MLEs (`ll_cr_nee` for cr, `DDD::dd_loglik` for dd). The verifying
measurement, 8 cells x 3 trees x 3 arms, 72 fits, deficits against MLEs recomputed from scratch
(they reproduce the study's reference to 2e-11 nats):

| arm | converged | median | p10 | worst | within 0.01 | within 0.1 | median iter |
|---|---|---|---|---|---|---|---|
| step (in force) | 0.708 | -0.0021 | -0.127 | -0.239 | 0.667 | 0.875 | 180 |
| ascent, mc_tol 0.002 | 0.500 | -0.0108 | -0.089 | -0.129 | 0.417 | 0.958 | 74 |
| ascent, mc_tol 0.01 | 1.000 | -0.0163 | -0.174 | -0.353 | 0.458 | 0.750 | 52 |

Cost, in E-step draws: +26 % over the sample, +62 % at the median fit, 26x at the worst pair.

## What it says

The two distributions cross. The ascent rule trims the bad tail and degrades the good middle: the
median fit lands five times further from the maximum, and the fraction within 0.01 nats falls from
0.67 to 0.42. `mc_tol` trades one against the other -- at 0.01 every fit reports "converged" and
the tail is the worst of the three arms.

The useful conclusion is about the study's own result rather than about either rule. **The cells
that failed the convergence criterion were not landing far from the maximum.** Under the rule in
force the median deficit over those same failing cells is -0.002 nats, which is well inside the
study's own 0.1-nat bar. What fails is the *declaration*, not the estimate: the rule is
conservative, runs to `max_iter`, and reports that it did not converge while sitting on the answer.
A rule that declares convergence more often is easy to obtain and is not by itself an improvement.

So the 51-of-124 convergence rate should be read as a property of the stopping rule's reporting,
not as evidence that those fits are unconverged. The report's recommendation to make convergence a
property of the Monte Carlo noise stands as a direction, but the gain it buys is in the tail, and
on this evidence it is paid for in the centre and in draws.

## What is worth keeping from the attempt

Two negative results, established by measurement during the design phase:

- Testing a *drift* of `Q` over a window against its noise looks like the natural statistic and is
  wrong: the per-iteration ascent is where the surrogate's bias is second order. The drift variant
  is what a reimplementation tries first and it fails in the way that most resembles success.
- Rebuilding the observed-data likelihood by importance-reweighting the same draws,
  `log mean exp(logf - logg)`, is biased upward near the optimum (an exact change of -1.7e-4 nats
  read as +1.2e-3) and its bootstrap spread can come out negative. The surrogate `Q` is the scale
  that behaves.

Also measured, and independently useful: `m_cpp` replayed on `em_cpp`'s returned trees and weights
reproduces `em_cpp`'s estimates bit-identically, so the M-step can be bootstrapped without new
draws on either sampler.

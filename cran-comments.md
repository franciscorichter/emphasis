## Resubmission

This replaces emphasis 1.0.0, submitted on 2026-10-22, which has not yet been
reviewed. It is submitted as version 1.1.0 because it adds functionality:

* `estimate_rates(cond = "surface")` conditions a fit on survival for any
  model and link, through local quadratic surfaces of a direct Monte Carlo
  estimate of the survival probability. Two new exported functions support
  it: `survival_mc()` and `survival_surface()`.
* A fix: under `stop_rule = "mc_error"` the thinning proposal's fits of the
  lineage-level models could never meet the rule and ran to `max_iter`.

## Test environments

* Local: macOS 26.6.2 (arm64), R 4.4.0, Apple clang 17.0.0:
  `R CMD check --as-cran --no-manual` on the built tarball.
* GitHub Actions (`r-lib/actions`), passing at the submitted commit:
  macOS (latest) R release; Windows (latest) R release; Ubuntu (latest)
  R devel, R release and R oldrel-1.

## R CMD check results

0 errors | 0 warnings | 4 notes locally (one from the package, three from
the test machine). Tests: 2274 passed, 0 failed, 92 skipped on CRAN.

* New submission.

* Possibly misspelled words in DESCRIPTION (flagged by CRAN's pre-tests of
  1.0.0): Etienne, Haegeman, phylogenies, speciation. These are spelled
  correctly: Etienne and Haegeman are the authors of the cited reference,
  and phylogenies and speciation are standard terms in evolutionary biology.

* Installed size: on Linux with R's default `-g` flag the shared library is
  about 11.7 Mb, which is debugging information; stripped of its debug
  sections it is 0.65 Mb. The package is C++ code built on the header
  libraries of Rcpp, RcppParallel (TBB), nloptr and BH.

## Notes from the local environment only

* "Package suggested but not available for checking: 'dglars'": not
  installed on the test machine; the tests that use it are skipped when it
  is absent.
* "unable to verify current time": the machine had no access to a time
  server during the check.
* "Examples with CPU (user + system) or elapsed time > 5s: compare_models
  (user 1.35 s, elapsed 6.75 s)": the CPU time is under 5 s; the elapsed
  time reflects other load on the test machine.

## Reverse dependencies

There are no reverse dependencies (new package).

## Threads, time and files

* Every function runs on one thread unless `num_threads` is set; tests
  that use more than two threads are skipped on CRAN.
* Examples run in under 5 s each; the longer ones are in `\donttest{}`.
* Nothing is written outside `tempdir()`.

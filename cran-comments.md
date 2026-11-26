## Submission

This is a patch release of emphasis 1.1.0. It fixes the test failure on
r-devel-linux-x86_64-fedora-clang reported on the CRAN check page (the CRAN
team's message of 2026-11-26, deadline 2026-12-17).

The failing test (`tests/testthat/test-reproducibility.R`) asserted that two
runs with `num_threads = 2` under one seed share at least one draw. That
depends on how the threads are scheduled, and on that platform the two runs
shared none. The test now asserts only what holds whatever the scheduling: a
full sample of finite draws. No code under `R/` or `src/` changed.

## Test environments

* Local: macOS (arm64), R 4.4.0: `R CMD check --as-cran` on the built tarball.

## R CMD check results

0 errors | 0 warnings | 1 note

* "checking for future file timestamps ... NOTE: unable to verify current time"
  comes from the local machine's network and does not concern the package.

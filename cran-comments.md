## Submission

First submission of emphasis (version 1.0.0).

## Test environments

* Pop!_OS 24.04 LTS (Ubuntu 24.04), x86_64, R 4.3.3, g++ 13.3.0:
  `R CMD check --as-cran --no-manual`.

To be run before submission: win-builder (R-devel and R-release) and
macOS (mac.r-project.org/macbuilder), and R-hub.

## R CMD check results

0 errors | 0 warnings | 5 notes (two from the package, three from the test
machine). Tests: 2270 passed, 80 skipped on CRAN, 181 s elapsed.

* New submission.

* Installed size is 12.4 Mb, of which the shared library is 11.6 Mb.
  The size is debugging information: R's default flags on the test machine
  include `-g`, and the library stripped of its debug sections
  (`strip --strip-debug`) is 0.67 Mb. The package is C++ code built on
  the header libraries of Rcpp, RcppParallel (TBB), nloptr and BH.

## Notes from the local environment only

These come from the test machine, not from the package:

* "Compilation used the following non-portable flag(s):
  '-mno-omit-leaf-frame-pointer'". The flag is in the Ubuntu build of R's
  `CXXFLAGS` (`R CMD config CXXFLAGS`); the package's `src/Makevars` and
  `src/Makevars.win` set only include paths, `-DRCPP_PARALLEL_USE_TBB=1`
  (Windows) and the RcppParallel link flags.
* "unable to verify current time": the machine had no access to a time
  server during the check.
* "Files 'README.md' or 'NEWS.md' cannot be checked without 'pandoc' being
  installed": pandoc is not installed on the test machine.

## Threads, time and files

* Every function runs on one thread unless `num_threads` is set; tests
  that use more than two threads are skipped on CRAN.
* Examples run in under 5 s each; the longer ones are in `\donttest{}`.
* Nothing is written outside `tempdir()`.

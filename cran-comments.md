## Resubmission

This replaces the emphasis 1.0.0 submitted on 2026-10-21, which has not yet
been reviewed. The only change is the removal of the experimental gaussian
link function (`link = "gaussian"`): it made the coefficients poorly
identified and duplicated the other links where it worked. The package now
offers the linear and exponential links. The checks below were rerun on the
resubmitted version.

## Submission

First submission of emphasis (version 1.0.0).

## Test environments

* Local: Pop!_OS 24.04 LTS (Ubuntu 24.04), x86_64, R 4.3.3, g++ 13.3.0:
  `R CMD check --as-cran --no-manual`.
* GitHub Actions (`r-lib/actions`), all passing at commit 63fadf6 and later:
  macOS (latest) R release; Windows (latest) R release; Ubuntu (latest)
  R devel, R release and R oldrel-1.
* win-builder: R-devel and R-release.
* macOS builder (mac.r-project.org/macbuilder, R 4.6.1 Patched, macOS 26.6,
  arm64): the check stopped at "checking package dependencies" because the
  builder had none of the CRAN packages emphasis imports (ape, DDD, nloptr,
  progress, Rcpp, RcppParallel, BH) available; the package itself was not
  checked there. macOS is covered by the GitHub Actions run above.

## R CMD check results

0 errors | 0 warnings | 5 notes (two from the package, three from the test
machine). CRAN's incoming pre-tests (Windows and Debian, r-devel) gave 1 note,
the spelling note below. Tests: 2198 passed, 0 failed, 87 skipped on CRAN.

* New submission.

* Possibly misspelled words in DESCRIPTION: Etienne, Haegeman, phylogenies,
  speciation. These are spelled correctly: Etienne and Haegeman are the
  authors of the cited reference, and phylogenies and speciation are standard
  terms in evolutionary biology.

* Installed size is 12.5 Mb, of which the shared library is 11.7 Mb.
  The size is debugging information: R's default flags on the test machine
  include `-g`, and the library stripped of its debug sections
  (`strip --strip-debug`) is 0.65 Mb. The package is C++ code built on
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

## Reverse dependencies

There are no reverse dependencies (new package).

## Threads, time and files

* Every function runs on one thread unless `num_threads` is set; tests
  that use more than two threads are skipped on CRAN.
* Examples run in under 5 s each; the longer ones are in `\donttest{}`.
* Nothing is written outside `tempdir()`.

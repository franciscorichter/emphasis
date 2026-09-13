## H72 — abort/assert entry points in compiled code (R CMD check WARNING).
## Read-only evidence collection: no compilation. Attributes the symbols to
## their headers and shows they are NDEBUG-gated, i.e. an artefact of the
## pkgbuild "debug build" (-UNDEBUG -O0) objects that were left in src/ and
## reused by both the March R CMD check and the scratch install.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
pkg <- "/Users/pancho/Code/emphasis"
so  <- file.path(pkg, "src", "emphasis.so")
sym <- system2("nm", c("-g", so), stdout = TRUE)
cat("Symbols in src/emphasis.so:\n")
print(grep("assert_rtn| _abort$", sym, value = TRUE))

## 1. Where the assert lives: Boost multiprecision (BH), not package code.
src_assert <- system2("grep", c("-rn", "'assert('", file.path(pkg, "src"),
                                file.path(pkg, "inst/include")), stdout = TRUE)
cat("\nassert( calls in package sources (excluding static_assert):\n")
print(grep("static_assert", src_assert, value = TRUE, invert = TRUE))  # expect character(0)

bh <- system.file("include/boost/assert.hpp", package = "BH")
cat("\nBoost assert.hpp default mapping:\n")
print(grep("define BOOST_ASSERT\\(expr\\) assert|defined\\(NDEBUG\\)", readLines(bh), value = TRUE))

## 2. Where the abort lives: Rcpp r_cast.h, inside '#ifndef NDEBUG'.
rc <- readLines(system.file("include/Rcpp/r_cast.h", package = "Rcpp"))
i  <- grep("abort\\(\\)", rc)
cat("\nRcpp r_cast.h abort() sites and their guards:\n")
for (k in i) cat(sprintf("  line %d: %s   [guard: %s]\n", k, trimws(rc[k]),
                         trimws(rc[max(which(grepl("#ifndef NDEBUG", rc[1:k])))])))

## 3. R's own compile flags always define NDEBUG; the objects in src/ were
##    produced by pkgbuild's debug build, which appends -UNDEBUG -O0.
mk <- readLines(file.path(R.home(), "etc", "Makeconf"))
cat("\nMakeconf ALL_CPPFLAGS:", grep("^ALL_CPPFLAGS", mk, value = TRUE), "\n")
log <- "/Users/pancho/.claude/jobs/867af780/tmp/emphasis-doc-test.log"
if (file.exists(log)) {
  l <- readLines(log)
  cat("pkgbuild line:", grep("Re-compiling", l, value = TRUE), "\n")
  fl <- strsplit(grep("^clang\\+\\+", l, value = TRUE)[1], " ")[[1]]
  cat("flags of interest:", fl[grepl("NDEBUG|^-O", fl)], "\n")
}
opt <- system2("dwarfdump", c("--debug-info", file.path(pkg, "src/E_step.o")), stdout = TRUE)
cat("DW_AT_APPLE_optimized units in E_step.o:", sum(grepl("DW_AT_APPLE_optimized", opt)),
    "(0 => built at -O0, i.e. the pkgbuild debug build)\n")

## 4. Reachability of the one assert from package code: calc_sum_w receives
##    only weights that passed std::isfinite (E_step.cpp:94), so
##    cpp_dec_float(double) never sees NaN/Inf from the E-step.
cat("\nE_step.cpp guard:", grep("isfinite", readLines(file.path(pkg, "src/E_step.cpp")), value = TRUE), "\n")

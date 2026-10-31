# Evaluate `expr`, muffling only the warnings whose message matches `pattern`;
# every other warning propagates.  Used where a test sets up the condition a
# warning reports (a starved fit, a stress configuration, tied counts) and the
# warning is the package doing its job, not a finding of the test.
muffle_expected <- function(expr, pattern) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl(pattern, conditionMessage(w))) invokeRestart("muffleWarning")
  })
}

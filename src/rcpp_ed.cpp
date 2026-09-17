#include <Rcpp.h>
#include <vector>
#include "ed_covariate.hpp"

//' Evolutionary distinctiveness of a lineage forest at a time
//'
//' The fair-proportion routine the likelihood and the simulator share
//' (\code{inst/include/ed_covariate.hpp}), exposed so that it can be checked
//' against an independent computation on an \code{ape} tree.
//'
//' @param parent Integer vector: for each lineage, the (1-based) index of the
//'   lineage it split from, or \code{0} for a crown lineage.
//' @param birth Numeric vector: each lineage's birth time (forward).
//' @param alive Logical vector: whether each lineage is alive on the segment
//'   being evaluated.
//' @param t The evaluation time.
//' @return A data frame with one row per lineage: \code{c} and \code{ts} such
//'   that \code{ED(u) = c + (u - ts)} for \code{u} in the segment (\code{NA}
//'   for a lineage not alive), \code{ed} the value at \code{t}, and
//'   \code{n_desc} the alive lineages in its subtree.
//' @keywords internal
// [[Rcpp::export(name = "ed_fair_proportion")]]
Rcpp::DataFrame rcpp_ed_fair_proportion(Rcpp::IntegerVector parent,
                                        Rcpp::NumericVector birth,
                                        Rcpp::LogicalVector alive,
                                        double t)
{
  const R_xlen_t n = parent.size();
  if (birth.size() != n || alive.size() != n) {
    Rcpp::stop("ed_fair_proportion: parent, birth and alive must have the same length");
  }
  std::vector<int> par(static_cast<std::size_t>(n));
  std::vector<double> b(static_cast<std::size_t>(n));
  std::vector<char> a(static_cast<std::size_t>(n));
  for (R_xlen_t i = 0; i < n; ++i) {
    par[static_cast<std::size_t>(i)] = parent[i] - 1;      // 0 -> -1 (no parent)
    b[static_cast<std::size_t>(i)] = birth[i];
    a[static_cast<std::size_t>(i)] = alive[i] ? 1 : 0;
  }
  emphasis::ed::result_t r;
  emphasis::ed::fair_proportion(par, b, a, t, r);
  Rcpp::NumericVector c(n), ts(n), ed(n);
  Rcpp::IntegerVector nd(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    const std::size_t k = static_cast<std::size_t>(i);
    c[i]  = r.c[k];
    ts[i] = r.ts[k];
    ed[i] = a[k] ? r.c[k] + (t - r.ts[k]) : NA_REAL;
    nd[i] = r.n_desc[k];
  }
  return Rcpp::DataFrame::create(Rcpp::Named("c") = c, Rcpp::Named("ts") = ts,
                                 Rcpp::Named("ed") = ed, Rcpp::Named("n_desc") = nd);
}

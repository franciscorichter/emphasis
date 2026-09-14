
#include <Rcpp.h>
#include <array>
#include <stdexcept>
#include <string>
#include <vector>
#include <tbb/tbb.h>
#include "emphasis.hpp"
#include "augment_tree.hpp"
#include "model.hpp"
#include "rinit.h"
#include "unpack.h"
using namespace Rcpp;

// src/loglik.cpp — an augmented-tree data frame back into a tree_t.
namespace loglik { emphasis::tree_t pack(const Rcpp::DataFrame& r_tree); }


//' Draw augmented trees via importance sampling
//'
//' Augments an observed extant tree (given by its branching times) with
//' stochastically drawn extinct lineages, and returns the simulated trees
//' together with their log proposal probabilities.  Scoring (\code{logf})
//' and IS aggregation (\code{fhat}) are intentionally left to the R layer
//' via \code{\link{eval_logf}} and \code{.is_fhat}.
//'
//' @param brts Numeric vector of branching times (crown age first,
//'   sorted decreasing).
//' @param pars Numeric vector of 8 model parameters
//'   \code{c(beta_0, beta_N, beta_M, beta_D, gamma_0, gamma_N, gamma_M, gamma_D)}.
//' @param sample_size Number of valid augmented trees to collect.
//' @param maxN Maximum total augmentation attempts (including failures).
//' @param max_missing Maximum extinct lineages per augmented tree.
//' @param max_lambda Upper bound on the speciation rate (thinning bound).
//' @param num_threads Threads for parallel augmentation.
//' @param model Integer vector \code{c(use_N, use_P, use_E)}.
//' @param link Link function: \code{0} = linear, \code{1} = exponential.
//' @param rho Sampling fraction in \code{(0, 1]}.
//' @param parent_tip_start Tip start of the lineage that splits at each
//'   observed branching event, in the same forward-time order as \code{brts}
//'   (that is, \code{rev} of the decreasing \code{brts}); one entry per event,
//'   or one per node with the last ignored. Empty (the default) means the
//'   topology is not available: every observed lineage is then recorded as
//'   dating from the crown and every observed event gets \code{D = 0}, which
//'   is what a bare branching-time vector has always produced.
//' @return A named list:
//' \describe{
//'   \item{trees}{List of augmented-tree data frames.}
//'   \item{logf}{Per-tree log p(obs, z | theta).}
//'   \item{logg}{Per-tree log q(z | obs, theta).}
//'   \item{rejected}{Unhandled rejections.}
//'   \item{rejected_overruns}{Rejected: too many extinct lineages.}
//'   \item{rejected_lambda}{Rejected: lambda bound exceeded.}
//'   \item{rejected_zero_weights}{Rejected: zero IS weight (log weight -Inf).}
//'   \item{rejected_nonfinite}{Rejected: log weight +Inf or NaN.}
//'   \item{num_trees}{Number of trees returned (\code{length(trees)}).}
//'   \item{envelope_violations}{Thinning candidates drawn during this call
//'     whose acceptance probability exceeded 1, i.e. the envelope did not
//'     dominate the rate. Non-zero means the draws are not from the density
//'     \code{logg} charges them; a warning is issued.}
//'   \item{time}{Elapsed time (ms).}
//' }
//' @keywords internal
// [[Rcpp::export(name = "augment_trees")]]
List rcpp_mce(const std::vector<double>& brts,
              const std::vector<double>& pars,
              int sample_size,
              int maxN,
              int max_missing,
              double max_lambda,
              int num_threads,
              Rcpp::IntegerVector model = Rcpp::IntegerVector::create(0, 0, 0),
              int link = 0,
              double rho = 1.0,
              Rcpp::NumericVector parent_tip_start = Rcpp::NumericVector::create())
{
  const std::vector<double> pts(parent_tip_start.begin(), parent_tip_start.end());
  if (pars.size() != 8) {
    throw std::invalid_argument("augment_trees: pars must have length 8 (got " +
      std::to_string(pars.size()) + ")");
  }
  std::vector<int> model_bin = {model[0], model[1], model[2]};
  // The proposal reads M = P/N at the segment start. Without the observed
  // topology P is the legacy quantity (every observed lineage dated from the
  // crown), which the scorer cannot reproduce from the finished tree, so
  // log q would not be the density the sampler drew from. cr, dd, d and nd
  // never reach this; an M-active model passed as a binary vector does.
  if (model_bin[1] == 1 && pts.empty()) {
    throw std::invalid_argument(
      "augment_trees: an M-dependent model (model[2] == 1) needs the observed topology. "
      "Pass parent_tip_start; without it log q is not the density the sampler "
      "draws from.");
  }
  // Bounds are only needed for M-step (nlopt); E-step does not use them.
  std::vector<double> lb8(8, -1e6), ub8(8, 1e6);
  auto mdl = emphasis::Model(lb8, ub8, model_bin, link, rho);

  // The envelope counter is process-wide and accumulates across calls and
  // threads; the difference over the call is this call's count, and leaves
  // the counter readable through thinning_envelope_violations().
  const long long viol_before = emphasis::thinning_envelope_violations(false);

  auto E = emphasis::E_step(sample_size,
                            maxN,
                            pars,
                            brts,
                            mdl,
                            max_missing,
                            max_lambda,
                            num_threads,
                            0.0,
                            pts);
  List ret;
  List trees;
  for (const emphasis::tree_t& tree : E.trees) {
    trees.push_back(unpack(tree));
  }
  ret["trees"]                 = trees;
  ret["logf"]                  = E.logf_;
  ret["logg"]                  = E.logg_;
  ret["rejected"]              = E.info.rejected;
  ret["rejected_overruns"]     = E.info.rejected_overruns;
  ret["rejected_lambda"]       = E.info.rejected_lambda;
  ret["rejected_zero_weights"] = E.info.rejected_zero_weights;
  ret["rejected_nonfinite"]    = E.info.rejected_nonfinite;
  ret["num_trees"]             = E.info.num_trees;
  const long long viol = emphasis::thinning_envelope_violations(false) - viol_before;
  ret["envelope_violations"]   = static_cast<double>(viol);
  ret["time"]                  = E.info.elapsed;
  if (viol > 0) {
    Rcpp::warning("augment_trees: thinning envelope did not dominate the rate at " +
      std::to_string(viol) + " candidate(s); those trees are drawn from a "
      "different density than logg charges");
  }
  return ret;
}


//' What the proposal could attach an augmented lineage to
//'
//' Replays an augmented tree through the same forward sweep the thinning
//' sampler carries, and reports at every augmented birth the lineages the
//' proposal had to choose from, the labelled attachments they carry, and the
//' parent that was recorded.  \code{candidates} must be the alive count
//' \code{n} and \code{attachments} the \code{2 * tips + Ne} that
//' \code{Model::sampling_prob} charges \code{-log} of; the two crown lineages
//' carry no node, so they appear only here and in their reserved ids.
//'
//' @param tree One augmented-tree data frame from \code{\link{augment_trees}}.
//' @return A data frame with one row per augmented birth: \code{brts},
//'   \code{n}, \code{candidates}, \code{attachments}, \code{parent_id} and
//'   \code{parent_alive}.
//' @keywords internal
// [[Rcpp::export(name = "eval_attachments")]]
DataFrame rcpp_attachments(const Rcpp::DataFrame& tree)
{
  auto local_tree = loglik::pack(tree);
  const auto rep = emphasis::attachment_report(local_tree);
  NumericVector brts(rep.size()), n(rep.size()), cand(rep.size()), att(rep.size());
  IntegerVector pid(rep.size());
  LogicalVector alive(rep.size());
  for (size_t i = 0; i < rep.size(); ++i) {
    brts[i] = rep[i].brts;
    n[i] = rep[i].n;
    cand[i] = static_cast<double>(rep[i].candidates);
    att[i] = static_cast<double>(rep[i].attachments);
    pid[i] = rep[i].parent_id;
    alive[i] = rep[i].parent_alive;
  }
  return DataFrame::create(Named("brts") = brts, Named("n") = n,
                           Named("candidates") = cand,
                           Named("attachments") = att,
                           Named("parent_id") = pid,
                           Named("parent_alive") = alive);
}


//' Count thinning candidates with acceptance probability above 1
//'
//' The thinning sampler accepts a candidate speciation time with probability
//' \code{nh(t) / lambda_max}. A value above 1 means the envelope
//' \code{lambda_max} did not dominate the rate on that segment. The counter
//' accumulates over every augmentation call in the session.
//'
//' @param reset Logical; zero the counter after reading it.
//' @return Number of candidates with acceptance probability above 1 since the
//'   last reset.
//' @keywords internal
// [[Rcpp::export(name = "thinning_envelope_violations")]]
double rcpp_thinning_envelope_violations(bool reset = false)
{
  return static_cast<double>(emphasis::thinning_envelope_violations(reset));
}

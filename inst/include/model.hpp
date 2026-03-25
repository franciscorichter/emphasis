/* Copyright (c) 2007-2014 Massachusetts Institute of Technology
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

//
// C++-API for plug-in diversification model
// Hanno Hildenbrandt 2020
//

#ifndef EMPHASIS_PLUGIN_HPP_INCLUDED
#define EMPHASIS_PLUGIN_HPP_INCLUDED

#include <vector>
#include <stdexcept>
#include <functional>
#include "model_helpers.hpp"

using namespace emphasis::detail;

namespace emphasis {

  using param_t = std::vector<double>;                  // unspecific parameters
  using tree_t = std::vector<node_t>;                   // tree, sorted by note_t::brts

  using reng_t = std::mt19937_64;   // we need doubles

  // Link functions for rate computation
  // 0 = linear:      rate = max(0, eta)
  // 1 = exponential: rate = exp(eta)
  enum class LinkType : int { linear = 0, exponential = 1 };

  // General diversification model
  //
  // Parameters layout (always 8 elements):
  //   {beta_0, beta_N, beta_P, beta_E, gamma_0, gamma_N, gamma_P, gamma_E}
  //
  // model_bin_ = {use_N, use_P, use_E} selects active covariates.
  // Inactive parameters should be 0 (and their bounds pinned to 0).
  //
  // Linear link (default):
  //   lambda(N,P) = max(0, beta_0  + beta_N*N  + beta_P*P)
  //   mu(N,P)     = max(0, gamma_0 + gamma_N*N + gamma_P*P)
  //
  // Exponential link:
  //   lambda(N,P) = exp(beta_0  + beta_N*N  + beta_P*P)
  //   mu(N,P)     = exp(gamma_0 + gamma_N*N + gamma_P*P)
  //
  class Model
  {
  public:
    // Full constructor: 8-element bounds + model binary vector + link type
    Model(const param_t& lb8, const param_t& ub8,
          const std::vector<int>& model_bin,
          int link = 0)
      : lower_bound_(lb8), upper_bound_(ub8), model_bin_(model_bin),
        link_(static_cast<LinkType>(link))
    {
      if (model_bin_.size() != 3) model_bin_ = {0, 0, 0};
    }

    // Legacy constructor: maps old 4-param rpd5c layout to 8-param
    Model(const param_t& lb, const param_t& ub)
      : model_bin_({1, 1, 0}), link_(LinkType::linear)
    {
      if (lb.size() == 8 && ub.size() == 8) {
        lower_bound_ = lb;
        upper_bound_ = ub;
      } else {
        lower_bound_ = {lb.size() > 1 ? lb[1] : 0, lb.size() > 2 ? lb[2] : 0,
                        lb.size() > 3 ? lb[3] : 0, 0,
                        lb.size() > 0 ? lb[0] : 0, 0, 0, 0};
        upper_bound_ = {ub.size() > 1 ? ub[1] : 0, ub.size() > 2 ? ub[2] : 0,
                        ub.size() > 3 ? ub[3] : 0, 0,
                        ub.size() > 0 ? ub[0] : 0, 0, 0, 0};
      }
    }

   ~Model() = default;

    const char* description() const { return "general diversification model"; }
    bool is_threadsafe() const { return true; }
    bool numerical_max_lambda() const { return false; }
    int nparams() const { return 8; }

    // Apply link function to linear predictor
    double apply_link(double eta) const {
      if (link_ == LinkType::exponential) {
        return std::exp(eta);
      }
      return std::max(0.0, eta);
    }

    // Speciation rate
    double speciation_rate(const param_t& pars, const node_t& node) const {
      const double eta = pars[0] + pars[1] * node.n + pars[2] * node.pd;
      return apply_link(eta);
    }

    // Extinction rate
    double extinction_rate(const param_t& pars, const node_t& node) const {
      const double eta = pars[4] + pars[5] * node.n + pars[6] * node.pd;
      return apply_link(eta);
    }

    // EP covariate E_s: exact pendant age of the focal lineage at event time.
    //
    // For extinction nodes: E = brts - tip_start (exact, stored t_spec).
    // For speciation nodes with known parent (parent_id >= 0):
    //   E = brts - focal_tip_start (exact parent pendant age).
    // For initial tree nodes (parent_id == -1, unknown parent):
    //   E = P/N (mean-field = correct marginalization over possible parents,
    //   since Σ_s λ(s)/N = λ(E=P/N) for linear link).
    double e_s(const node_t& node) const {
      if (detail::is_extinction(node)) {
        return node.brts - node.tip_start;
      }
      if (node.parent_id >= 0) {
        return node.brts - node.focal_tip_start;
      }
      return (node.n > 0.0) ? (node.pd / node.n) : 0.0;
    }

    // EP-aware speciation rate (mean-field E for non-extinction nodes)
    double speciation_rate_ep(const param_t& pars, const node_t& node) const {
      const double E = e_s(node);
      return apply_link(pars[0] + pars[1] * node.n + pars[2] * node.pd + pars[3] * E);
    }

    // EP-aware extinction rate (exact E for extinction nodes, mean-field otherwise)
    double extinction_rate_ep(const param_t& pars, const node_t& node) const {
      const double E = e_s(node);
      return apply_link(pars[4] + pars[5] * node.n + pars[6] * node.pd + pars[7] * E);
    }

    // Draw extinction time from truncated exponential with rate = mu at speciation time
    double extinction_time(double t_speciation, const param_t& pars, const tree_t& tree) const {
      static thread_local reng_t reng_ = make_random_engine<reng_t>();
      auto it = lower_bound_node(t_speciation, tree.size(),
                                 reinterpret_cast<const node_t*>(tree.data()));
      double mu = extinction_rate(pars, *it);
      if (mu <= 0.0) mu = 1e-10;
      return t_speciation + emphasis::detail::trunc_exp(tree.back().brts - t_speciation, mu, reng_);
    }

    // Non-homogeneous rate for thinning: N * lambda * (1 - exp(-mu*(T-t)))
    double nh_rate(double t, const param_t& pars, const tree_t& tree) const {
      auto it = lower_bound_node(t,
                                 tree.size(),
                                 reinterpret_cast<const node_t*>(tree.data()));
      const double pd = calculate_pendant_pd(t, tree);
      node_t cur = *it;
      cur.pd = pd;
      const double lambda = model_bin_[2] ? speciation_rate_ep(pars, cur)
                                          : speciation_rate(pars, cur);
      const double mu     = model_bin_[2] ? extinction_rate_ep(pars, cur)
                                          : extinction_rate(pars, cur);
      const double T = tree.back().brts;
      const double mu_eff = std::max(mu, 1e-10);
      return lambda * it->n * (1.0 - std::exp(-mu_eff * (T - t)));
    }


    double sampling_prob(const param_t& pars, const tree_t& tree) const {
      const bool ep_exp = model_bin_[2] && (link_ == LinkType::exponential);

      double inte = 0;
      double logg = 0;
      double prev_brts = 0;
      double tips = tree[0].n;
      double Ne = 0.0;
      const double T = tree.back().brts;

      // Running sums for EP+exp
      double sum_exp_bE = ep_exp ? tree[0].n : 0.0;
      double sum_exp_gE = ep_exp ? tree[0].n : 0.0;

      for (unsigned i = 0; i < tree.size(); ++i) {
        const auto& node = tree[i];
        const double lambda = model_bin_[2] ? speciation_rate_ep(pars, node)
                                            : speciation_rate(pars, node);
        const double mu = std::max(model_bin_[2] ? extinction_rate_ep(pars, node)
                                                 : extinction_rate(pars, node), 1e-10);
        {
          double dt = node.brts - prev_brts;
          if (ep_exp && dt > 0.0) {
            // For EP+exp, the thinning integral is:
            // ∫ Σ_s λ_s(t) * (1 - exp(-μ_s(t)*(T-t))) dt
            // This doesn't have a simple closed form when both λ and μ
            // depend on pendant age. Use the same mean-field approximation
            // for the survival term but exact integration for the rate part.
            // Approximation: use node-level mu for the survival factor
            double int_segment = dt - (1.0/mu) * (std::exp(-mu*(T - node.brts)) - std::exp(-mu*(T - prev_brts)));
            // Use exact EP+exp rate sum for lambda * n part
            double N_seg = node.n;
            double PD_seg = node.pd;
            double A_lam = pars[0] + pars[1]*N_seg + pars[2]*PD_seg;
            // Mean lambda over the segment using exact integral
            double lam_integral = sum_exp_bE * exp_integral(A_lam, pars[3], prev_brts, node.brts);
            double mean_lam_n = (dt > 0.0) ? lam_integral / dt : N_seg * lambda;
            inte += mean_lam_n * int_segment;
          } else {
            double int_segment = dt - (1.0/mu) * (std::exp(-mu*(T - node.brts)) - std::exp(-mu*(T - prev_brts)));
            inte += node.n * lambda * int_segment;
          }
        }
        tips += is_tip(node);
        Ne -= is_extinction(node);
        if (is_missing(node)) {
          const double lifespan = node.t_ext - node.brts;
          logg += std::log(node.n * mu * lambda) - mu * lifespan - std::log(2.0 * tips + Ne++);
        }

        // Update running sums
        if (ep_exp) {
          if (is_missing(node)) {
            sum_exp_bE += std::exp(-pars[3] * node.brts);
            sum_exp_gE += std::exp(-pars[7] * node.brts);
          } else if (is_extinction(node)) {
            sum_exp_bE -= std::exp(-pars[3] * node.tip_start);
            sum_exp_gE -= std::exp(-pars[7] * node.tip_start);
          } else if (i != tree.size() - 1) {
            sum_exp_bE += std::exp(-pars[3] * node.brts);
            sum_exp_gE += std::exp(-pars[7] * node.brts);
          }
        }

        prev_brts = node.brts;
      }
      return logg - inte;
    }

    // Helper: integral of exp(a + b*t) from t1 to t2
    // = exp(a) * [exp(b*t2) - exp(b*t1)] / b
    // Handles b→0 limit: exp(a) * (t2 - t1)
    static double exp_integral(double a, double b, double t1, double t2) {
      if (std::abs(b) < 1e-12) return std::exp(a) * (t2 - t1);
      return std::exp(a) * (std::exp(b * t2) - std::exp(b * t1)) / b;
    }

    double loglik(const param_t& pars, const tree_t& tree) const {
      const bool ep_exp = model_bin_[2] && (link_ == LinkType::exponential);

      log_sum log_lambda{};
      double log_mu_sum = 0.0;
      double inte = 0.0;
      double prev_brts = 0.0;

      // For EP+exp: running sums of exp(-β_E * ts_s) over active lineages
      double sum_exp_bE = 0.0;  // Σ_s exp(-pars[3] * ts_s)
      double sum_exp_gE = 0.0;  // Σ_s exp(-pars[7] * ts_s)

      if (ep_exp) {
        // Initialize running sums from initial lineages
        // At t=0, all lineages have tip_start = 0 for the observed tree
        // The initial n lineages contribute n * exp(0) = n
        sum_exp_bE = tree[0].n;
        sum_exp_gE = tree[0].n;
      }

      for (unsigned i = 0; i < tree.size(); ++i) {
        const auto& node = tree[i];
        const double dt = node.brts - prev_brts;

        if (ep_exp && dt > 0.0) {
          // Exact integral: ∫_{prev}^{brts} Σ_s [λ_s(t) + μ_s(t)] dt
          // λ_s(t) = exp(β₀ + β_N·N + β_P·PD + β_E·t) · exp(-β_E·ts_s)
          // ∫ Σ_s λ_s(t) dt = sum_exp_bE · ∫ exp(A + β_E·t) dt
          //   where A = β₀ + β_N·N + β_P·PD (constant in segment)
          //
          // N, PD are from previous node (constant in segment)
          double N_seg = node.n;
          double PD_seg = node.pd;  // PD at segment start
          double A_lam = pars[0] + pars[1]*N_seg + pars[2]*PD_seg;
          double A_mu  = pars[4] + pars[5]*N_seg + pars[6]*PD_seg;
          inte += sum_exp_bE * exp_integral(A_lam, pars[3], prev_brts, node.brts)
                + sum_exp_gE * exp_integral(A_mu,  pars[7], prev_brts, node.brts);
        } else {
          // Standard: piecewise constant rates
          const double lambda = model_bin_[2] ? speciation_rate_ep(pars, node)
                                              : speciation_rate(pars, node);
          const double mu     = model_bin_[2] ? extinction_rate_ep(pars, node)
                                              : extinction_rate(pars, node);
          inte += dt * node.n * (lambda + mu);
        }

        // Event contributions (speciation / extinction log-probabilities)
        const double lambda = model_bin_[2] ? speciation_rate_ep(pars, node)
                                            : speciation_rate(pars, node);
        const double mu     = model_bin_[2] ? extinction_rate_ep(pars, node)
                                            : extinction_rate(pars, node);

        if (is_extinction(node)) {
          log_mu_sum += std::log(std::max(mu, 1e-300));
        }
        else if (i != tree.size() - 1) {
          log_lambda += lambda;
        }

        // Update running sums for new/removed lineages
        if (ep_exp) {
          if (is_missing(node)) {
            // New augmented lineage born at node.brts
            sum_exp_bE += std::exp(-pars[3] * node.brts);
            sum_exp_gE += std::exp(-pars[7] * node.brts);
          } else if (is_extinction(node)) {
            // Lineage with tip_start removed
            sum_exp_bE -= std::exp(-pars[3] * node.tip_start);
            sum_exp_gE -= std::exp(-pars[7] * node.tip_start);
          } else if (!is_extinction(node) && i != tree.size() - 1) {
            // Observed speciation: daughter born at node.brts
            sum_exp_bE += std::exp(-pars[3] * node.brts);
            sum_exp_gE += std::exp(-pars[7] * node.brts);
          }
        }

        prev_brts = node.brts;
      }
      const double loglik = log_mu_sum + log_lambda.result() - inte;
      return loglik;
    }

    param_t lower_bound() const { return lower_bound_; }
    param_t upper_bound() const { return upper_bound_; }
  private:
    param_t lower_bound_;
    param_t upper_bound_;
    std::vector<int> model_bin_ = {0, 0, 0};
    LinkType link_ = LinkType::linear;
  };
}

#endif

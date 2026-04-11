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
  // 2 = gaussian:    rate = beta_0 * exp(-(eta_cov - 1)^2 / 2)
  enum class LinkType : int { linear = 0, exponential = 1, gaussian = 2 };

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
  // Gaussian link:
  //   lambda(N,P) = beta_0  * exp(-( beta_N*N +  beta_P*P +  beta_E*E - 1)^2 / 2)
  //   mu(N,P)     = gamma_0 * exp(-(gamma_N*N + gamma_P*P + gamma_E*E - 1)^2 / 2)
  //
  class Model
  {
  public:
    // Full constructor: 8-element bounds + model binary vector + link type + rho
    Model(const param_t& lb8, const param_t& ub8,
          const std::vector<int>& model_bin,
          int link = 0,
          double rho = 1.0)
      : lower_bound_(lb8), upper_bound_(ub8), model_bin_(model_bin),
        link_(static_cast<LinkType>(link)), rho_(rho)
    {
      if (model_bin_.size() != 3) model_bin_ = {0, 0, 0};
      if (rho_ <= 0.0 || rho_ > 1.0) rho_ = 1.0;
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
      // For gaussian, apply_link should not be called directly — use gaussian_rate
      return std::max(0.0, eta);
    }

    // Gaussian link: rate = intercept * exp(-(eta_cov - 1)^2 / 2)
    // where eta_cov = covariate terms only (no intercept)
    double gaussian_rate(double intercept, double eta_cov) const {
      const double d = eta_cov - 1.0;
      return intercept * std::exp(-0.5 * d * d);
    }

    // Speciation rate
    double speciation_rate(const param_t& pars, const node_t& node) const {
      if (link_ == LinkType::gaussian) {
        const double eta_cov = pars[1] * node.n + pars[2] * node.pd;
        return gaussian_rate(pars[0], eta_cov);
      }
      const double eta = pars[0] + pars[1] * node.n + pars[2] * node.pd;
      return apply_link(eta);
    }

    // Extinction rate
    double extinction_rate(const param_t& pars, const node_t& node) const {
      if (link_ == LinkType::gaussian) {
        const double eta_cov = pars[5] * node.n + pars[6] * node.pd;
        return gaussian_rate(pars[4], eta_cov);
      }
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
      if (link_ == LinkType::gaussian) {
        const double eta_cov = pars[1] * node.n + pars[2] * node.pd + pars[3] * E;
        return gaussian_rate(pars[0], eta_cov);
      }
      return apply_link(pars[0] + pars[1] * node.n + pars[2] * node.pd + pars[3] * E);
    }

    // EP-aware extinction rate (exact E for extinction nodes, mean-field otherwise)
    double extinction_rate_ep(const param_t& pars, const node_t& node) const {
      const double E = e_s(node);
      if (link_ == LinkType::gaussian) {
        const double eta_cov = pars[5] * node.n + pars[6] * node.pd + pars[7] * E;
        return gaussian_rate(pars[4], eta_cov);
      }
      return apply_link(pars[4] + pars[5] * node.n + pars[6] * node.pd + pars[7] * E);
    }

    // Draw extinction time from truncated exponential with rate = mu at speciation time.
    // With rho < 1, the lineage may be an unsampled extant species (returns t > T).
    // Caller should check: if result >= T, treat as unsampled tip (t_ext = t_ext_tip).
    double extinction_time(double t_speciation, const param_t& pars, const tree_t& tree) const {
      static thread_local reng_t reng_ = make_random_engine<reng_t>();
      auto it = lower_bound_node(t_speciation, tree.size(),
                                 reinterpret_cast<const node_t*>(tree.data()));
      double mu = extinction_rate(pars, *it);
      if (mu <= 0.0) mu = 1e-10;
      const double T = tree.back().brts;
      const double remaining = T - t_speciation;
      if (rho_ < 1.0) {
        // P(unsampled extant) = (1-rho)*exp(-mu*remaining) / (1 - rho*exp(-mu*remaining))
        double surv = std::exp(-mu * remaining);
        double p_unsampled = (1.0 - rho_) * surv / (1.0 - rho_ * surv);
        double u = std::uniform_real_distribution<>()(reng_);
        if (u < p_unsampled) {
          return T + 1.0;  // sentinel: unsampled extant species
        }
      }
      return t_speciation + emphasis::detail::trunc_exp(remaining, mu, reng_);
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
      return lambda * it->n * (1.0 - rho_ * std::exp(-mu_eff * (T - t)));
    }


    double sampling_prob(const param_t& pars, const tree_t& tree) const {
      const bool ep_exp = model_bin_[2] && (link_ == LinkType::exponential);
      const bool ep_gauss = model_bin_[2] && (link_ == LinkType::gaussian);

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
            // ∫ Σ_s λ_s(t) * (1 - ρ*exp(-μ_s(t)*(T-t))) dt
            // Approximation: use node-level mu for the survival factor
            double int_segment = dt - rho_ * (1.0/mu) * (std::exp(-mu*(T - node.brts)) - std::exp(-mu*(T - prev_brts)));
            // Use exact EP+exp rate sum for lambda * n part
            double N_seg = node.n;
            double PD_seg = node.pd;
            double A_lam = pars[0] + pars[1]*N_seg + pars[2]*PD_seg;
            // Mean lambda over the segment using exact integral
            double lam_integral = sum_exp_bE * exp_integral(A_lam, pars[3], prev_brts, node.brts);
            double mean_lam_n = (dt > 0.0) ? lam_integral / dt : N_seg * lambda;
            inte += mean_lam_n * int_segment;
          } else {
            double int_segment = dt - rho_ * (1.0/mu) * (std::exp(-mu*(T - node.brts)) - std::exp(-mu*(T - prev_brts)));
            inte += node.n * lambda * int_segment;
          }
        }
        tips += is_tip(node);
        Ne -= is_extinction(node);
        if (is_missing(node)) {
          // Augmented extinct species: logg = log(N*λ*μ) - μ*lifespan - log(K)
          const double lifespan = node.t_ext - node.brts;
          logg += std::log(node.n * mu * lambda) - mu * lifespan - std::log(2.0 * tips + Ne++);
        }
        else if (is_unsampled(node)) {
          // Unsampled extant species (rho < 1): logg = log(N*λ*(1-ρ)) - μ*(T-t) - log(K)
          logg += std::log(node.n * lambda * (1.0 - rho_)) - mu * (T - node.brts) - std::log(2.0 * tips + Ne++);
        }

        // Update running sums
        if (ep_exp) {
          if (is_missing(node) || is_unsampled(node)) {
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

    // Helper: integral of intercept * exp(-(A + b*t - 1)^2 / 2) from t1 to t2
    // where A = covariate sum excluding E, b = beta_E coefficient.
    // Uses error function: = intercept * sqrt(pi/2) / |b| * [erf(z2) - erf(z1)]
    // where z_i = (A + b*t_i - 1) / sqrt(2)
    // Handles b→0 limit: intercept * exp(-(A-1)^2/2) * (t2 - t1)
    static double gauss_integral(double intercept, double A, double b, double t1, double t2) {
      if (std::abs(b) < 1e-12) {
        double d = A - 1.0;
        return intercept * std::exp(-0.5 * d * d) * (t2 - t1);
      }
      const double inv_sqrt2 = 1.0 / std::sqrt(2.0);
      double z1 = (A + b * t1 - 1.0) * inv_sqrt2;
      double z2 = (A + b * t2 - 1.0) * inv_sqrt2;
      return intercept * std::sqrt(M_PI / 2.0) / std::abs(b) * (std::erf(z2) - std::erf(z1));
    }

    double loglik(const param_t& pars, const tree_t& tree) const {
      const bool ep_exp = model_bin_[2] && (link_ == LinkType::exponential);
      const bool ep_gauss = model_bin_[2] && (link_ == LinkType::gaussian);

      log_sum log_lambda{};
      double log_mu_sum = 0.0;
      double inte = 0.0;
      double prev_brts = 0.0;

      // For EP+exp: running sums of exp(-β_E * ts_s) over active lineages
      double sum_exp_bE = 0.0;  // Σ_s exp(-pars[3] * ts_s)
      double sum_exp_gE = 0.0;  // Σ_s exp(-pars[7] * ts_s)

      if (ep_exp) {
        sum_exp_bE = tree[0].n;
        sum_exp_gE = tree[0].n;
      }

      for (unsigned i = 0; i < tree.size(); ++i) {
        const auto& node = tree[i];
        const double dt = node.brts - prev_brts;

        if (ep_exp && dt > 0.0) {
          double N_seg = node.n;
          double PD_seg = node.pd;
          double A_lam = pars[0] + pars[1]*N_seg + pars[2]*PD_seg;
          double A_mu  = pars[4] + pars[5]*N_seg + pars[6]*PD_seg;
          inte += sum_exp_bE * exp_integral(A_lam, pars[3], prev_brts, node.brts)
                + sum_exp_gE * exp_integral(A_mu,  pars[7], prev_brts, node.brts);
        } else if (ep_gauss && dt > 0.0) {
          // EP + Gaussian: per-lineage exact integrals via erf.
          // For each alive lineage s with tip_start ts_s:
          //   λ_s(t) = β₀ · exp(-(β_N·N + β_P·P + β_E·(t - ts_s) - 1)² / 2)
          //          = β₀ · exp(-((A_lam + β_E·t) - β_E·ts_s - 1)² / 2)
          //          = β₀ · exp(-((A_lam - β_E·ts_s) + β_E·t - 1)² / 2)
          // So for each lineage, A_eff = A_lam_cov - β_E·ts_s
          // We need to iterate over alive lineages (no factorization possible)
          // But we can still use gauss_integral per lineage.
          // NOTE: we don't have a direct list of alive lineages with tip_starts here,
          // so we fall through to piecewise-constant approximation using the EP rate.
          const double lambda = speciation_rate_ep(pars, node);
          const double mu     = extinction_rate_ep(pars, node);
          inte += dt * node.n * (lambda + mu);
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
          if (is_missing(node) || is_unsampled(node)) {
            sum_exp_bE += std::exp(-pars[3] * node.brts);
            sum_exp_gE += std::exp(-pars[7] * node.brts);
          } else if (is_extinction(node)) {
            sum_exp_bE -= std::exp(-pars[3] * node.tip_start);
            sum_exp_gE -= std::exp(-pars[7] * node.tip_start);
          } else if (!is_extinction(node) && i != tree.size() - 1) {
            sum_exp_bE += std::exp(-pars[3] * node.brts);
            sum_exp_gE += std::exp(-pars[7] * node.brts);
          }
        }

        prev_brts = node.brts;
      }
      // Sampling terms for incomplete taxon sampling (rho < 1)
      // Each observed tip contributes log(rho), each unsampled extant contributes log(1-rho)
      double log_sampling = 0.0;
      if (rho_ < 1.0) {
        int n_obs = 1;  // start at 1: crown-age tree has n tips but only n-1 is_tip events
        int n_unsamp = 0;
        for (const auto& node : tree) {
          if (is_tip(node)) ++n_obs;
          else if (is_unsampled(node)) ++n_unsamp;
        }
        log_sampling = n_obs * std::log(rho_) + n_unsamp * std::log(1.0 - rho_);
      }

      const double loglik = log_mu_sum + log_lambda.result() - inte + log_sampling;
      return loglik;
    }

    param_t lower_bound() const { return lower_bound_; }
    param_t upper_bound() const { return upper_bound_; }
    double rho() const { return rho_; }
  private:
    param_t lower_bound_;
    param_t upper_bound_;
    std::vector<int> model_bin_ = {0, 0, 0};
    LinkType link_ = LinkType::linear;
    double rho_ = 1.0;
  };
}

#endif

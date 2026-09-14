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
#include <atomic>
#include <cstdint>
#include "model_helpers.hpp"

using namespace emphasis::detail;

namespace emphasis {

  using param_t = std::vector<double>;                  // unspecific parameters
  using tree_t = std::vector<node_t>;                   // tree, sorted by note_t::brts

  using reng_t = std::mt19937_64;   // the one engine type; we need doubles


  // -------------------------------------------------------------------------
  // Seeding of the C++ samplers.
  //
  // The thinning augmenter (src/augment_tree.cpp), the extinction-time draw in
  // Model::extinction_time below and the forward simulator
  // (inst/include/general_tree.hpp) each held an engine seeded from the wall
  // clock XOR the thread id.  set.seed() did not reach any of them -- 0 of 4
  // repeats identical against 4 of 4 for the pure-R BDI path -- and a
  // thread_local engine is copied byte for byte into every forked child, so
  // the children of a primed parent drew the same numbers rather than
  // independent ones (H47).
  //
  // An Rcpp entry point now opens a stream epoch with rng::set_seed(s), where
  // s is the integer the R layer passed: sample.int(.Machine$integer.max, 1L)
  // by default, so set.seed() propagates and a forked child, whose R stream
  // differs from its parent's, draws its own s.  Within an epoch each thread
  // draws from one substream, seeded from a std::seed_seq over (s, substream);
  // the substream is the TBB worker index, or an item index where the caller
  // has one, which is what keeps parallel workers independent of each other
  // while every one of them stays a function of s alone.
  //
  // Re-seeding happens when the epoch or the substream changes, not once per
  // call: re-seeding a thread to its own substream on every call would hand
  // every call on that thread the same draws.
  namespace rng {

    namespace detail_rng {

      inline std::atomic<uint64_t>& seed_slot()
      {
        static std::atomic<uint64_t> s{ 0 };
        return s;
      }

      // Bumped by every set_seed, so that a thread whose engine was seeded in
      // an earlier epoch re-seeds on its next draw.  Not part of the seed: two
      // calls carrying the same s must produce the same stream.
      inline std::atomic<uint64_t>& epoch_slot()
      {
        static std::atomic<uint64_t> e{ 0 };
        return e;
      }

      struct stream_t
      {
        reng_t eng{};
        uint64_t epoch = ~uint64_t(0);   // never seeded
        uint64_t sub = 0;
      };

      inline stream_t& local()
      {
        static thread_local stream_t s;
        return s;
      }

      // std::seed_seq consumes 32-bit words, so each 64-bit input enters as
      // its low and its high half.
      inline reng_t make(uint64_t seed, uint64_t sub)
      {
        const uint32_t w[5] = {
          static_cast<uint32_t>(seed), static_cast<uint32_t>(seed >> 32),
          static_cast<uint32_t>(sub),  static_cast<uint32_t>(sub >> 32),
          0x9e3779b9u
        };
        std::seed_seq sseq(w, w + 5);
        return reng_t(sseq);
      }

    }


    // Open a stream epoch: every thread re-seeds from `seed` on its next draw.
    inline void set_seed(uint64_t seed)
    {
      detail_rng::seed_slot().store(seed);
      detail_rng::epoch_slot().fetch_add(1);
    }


    // Put the calling thread on substream `sub` of the current epoch and hand
    // back its engine.  Callers that know which worker or which item they are
    // select the substream here; everything below them draws from engine().
    inline reng_t& stream(uint64_t sub)
    {
      auto& s = detail_rng::local();
      const uint64_t ep = detail_rng::epoch_slot().load();
      if (s.epoch != ep || s.sub != sub) {
        s.eng = detail_rng::make(detail_rng::seed_slot().load(), sub);
        s.epoch = ep;
        s.sub = sub;
      }
      return s.eng;
    }


    // The calling thread's engine, on the substream it is already on.
    inline reng_t& engine()
    {
      return stream(detail_rng::local().sub);
    }

  }


  // Link functions for rate computation
  // 0 = linear:      rate = max(0, eta)
  // 1 = exponential: rate = exp(eta)
  // 2 = gaussian:    rate = beta_0 * exp(-(eta_cov - 1)^2 / 2)
  enum class LinkType : int { linear = 0, exponential = 1, gaussian = 2 };

  // General diversification model
  //
  // Orthogonal covariate basis (per lineage s at event time t):
  //   N       = number of alive lineages          (diversity level)
  //   M = P/N = mean pendant age of alive lineages (maturity; P = sum of pendant ages)
  //   D = E - M                                    (focal isolation relative to the
  //                                                 clade mean; Sum_s D = 0)
  //
  // A lineage's pendant age at t is t - tip_start, where tip_start is the time
  // it last became a pendant tip: its birth time, reset every time it
  // speciates.  P = N*t - Sum_s tip_start_s is accumulated by the forward
  // sweep in compute_pendant_pd() (src/augment_tree.cpp) and read off
  // node.pd; E comes from node.focal_tip_start via e_s().  Both are the
  // definitions the simulator uses (inst/include/general_tree.hpp).
  //
  // Parameters layout (always 8 elements); slots 2/6 are coefficients on M and
  // slots 3/7 on D:
  //   {beta_0, beta_N, beta_M, beta_D, gamma_0, gamma_N, gamma_M, gamma_D}
  //
  // model_bin_ = {use_N, use_M, use_D} selects active covariates.
  // Inactive parameters should be 0 (and their bounds pinned to 0).
  //
  // Linear link (default):
  //   lambda = max(0, beta_0  + beta_N*N  + beta_M*M  + beta_D*D)
  //   mu     = max(0, gamma_0 + gamma_N*N + gamma_M*M + gamma_D*D)
  //
  // Exponential link:
  //   lambda = exp(beta_0  + beta_N*N  + beta_M*M  + beta_D*D)
  //   mu     = exp(gamma_0 + gamma_N*N + gamma_M*M + gamma_D*D)
  //
  // Gaussian link (eta_cov excludes the intercept):
  //   lambda = beta_0  * exp(-( beta_N*N +  beta_M*M +  beta_D*D - 1)^2 / 2)
  //   mu     = gamma_0 * exp(-(gamma_N*N + gamma_M*M + gamma_D*D - 1)^2 / 2)
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

    // ---------------------------------------------------------------------
    // Orthogonal covariate basis {N, M, D}:
    //   N       = node.n                     (lineage count / diversity level)
    //   M = P/N = node.pd / node.n           (mean pendant age / "maturity",
    //                                         over the N lineages alive on the
    //                                         segment that ends at the node)
    //   D = E - M                            (focal lineage's isolation
    //                                         relative to the clade mean; Σ_s D = 0)
    // The parameter slots are reused: pars[2]/pars[6] are now coefficients on M
    // (beta_M / gamma_M) and pars[3]/pars[7] on D (beta_D / gamma_D).
    // ---------------------------------------------------------------------

    // Mean pendant age M = P/N (global covariate).
    double m_cov(const node_t& node) const {
      return (node.n > 0.0) ? (node.pd / node.n) : 0.0;
    }

    // Speciation rate (no-D path: eta = beta_0 + beta_N*N + beta_M*M)
    double speciation_rate(const param_t& pars, const node_t& node) const {
      const double M = m_cov(node);
      if (link_ == LinkType::gaussian) {
        const double eta_cov = pars[1] * node.n + pars[2] * M;
        return gaussian_rate(pars[0], eta_cov);
      }
      const double eta = pars[0] + pars[1] * node.n + pars[2] * M;
      return apply_link(eta);
    }

    // Extinction rate (no-D path)
    double extinction_rate(const param_t& pars, const node_t& node) const {
      const double M = m_cov(node);
      if (link_ == LinkType::gaussian) {
        const double eta_cov = pars[5] * node.n + pars[6] * M;
        return gaussian_rate(pars[4], eta_cov);
      }
      const double eta = pars[4] + pars[5] * node.n + pars[6] * M;
      return apply_link(eta);
    }

    // Raw focal isolation E_s = pendant age of the focal lineage at event time.
    //
    // For extinction nodes: E = brts - tip_start (the dying lineage's own
    //   tip_start, reset by any split it made since its birth).
    // For an event whose splitting lineage is on record: E = brts -
    //   focal_tip_start.  That covers every observed branching event of a tree
    //   passed with its topology and every augmented lineage drawn from a
    //   recorded parent.
    // For an event whose splitting lineage is not on record (focal_tip_start
    //   is ts_unknown: a branching-time vector with no topology, or a lineage
    //   drawn when nothing was on record): E = P/N = M, a mean-field
    //   marginalization over possible parents, which sets D = E - M = 0.
    double e_s(const node_t& node) const {
      if (detail::is_extinction(node)) {
        return node.brts - node.tip_start;
      }
      if (node.focal_tip_start >= 0.0) {
        return node.brts - node.focal_tip_start;
      }
      return (node.n > 0.0) ? (node.pd / node.n) : 0.0;
    }

    // Centered focal deviation D = E - M.
    double d_cov(const node_t& node) const {
      return e_s(node) - m_cov(node);
    }

    // D-aware speciation rate: eta = beta_0 + beta_N*N + beta_M*M + beta_D*D
    double speciation_rate_ep(const param_t& pars, const node_t& node) const {
      const double M = m_cov(node);
      const double D = e_s(node) - M;
      if (link_ == LinkType::gaussian) {
        const double eta_cov = pars[1] * node.n + pars[2] * M + pars[3] * D;
        return gaussian_rate(pars[0], eta_cov);
      }
      return apply_link(pars[0] + pars[1] * node.n + pars[2] * M + pars[3] * D);
    }

    // D-aware extinction rate (exact D for extinction nodes, mean-field D=0 otherwise)
    double extinction_rate_ep(const param_t& pars, const node_t& node) const {
      const double M = m_cov(node);
      const double D = e_s(node) - M;
      if (link_ == LinkType::gaussian) {
        const double eta_cov = pars[5] * node.n + pars[6] * M + pars[7] * D;
        return gaussian_rate(pars[4], eta_cov);
      }
      return apply_link(pars[4] + pars[5] * node.n + pars[6] * M + pars[7] * D);
    }

    // Draw extinction time from truncated exponential with rate = mu at speciation time.
    // With rho < 1, the lineage may be an unsampled extant species (returns t > T).
    // Caller should check: if result >= T, treat as unsampled tip (t_ext = t_ext_tip).
    //
    // The rate is the segment's mu (proposal_rates), the same one the survival
    // factor of nh_rate used to decide that the lineage is missing at all and
    // the same one sampling_prob charges the lifetime with.
    double extinction_time(double t_speciation, const param_t& pars, const tree_t& tree) const {
      // The stream the augmentation that called this is already drawing from,
      // so the lifetime and the birth times it goes with come off one seeded
      // engine rather than off a second, clock-seeded one.
      reng_t& reng_ = rng::engine();
      const node_t* first = reinterpret_cast<const node_t*>(tree.data());
      auto it = lower_bound_node(t_speciation, tree.size(), first);
      double lambda = 0.0, mu = 0.0;
      proposal_rates(pars, *it, segment_start(it, first), lambda, mu);
      (void)lambda;
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

    // Pendant PD at an arbitrary time inside the segment `node` governs.
    //
    // node.pd is P at the node's own time over the N = node.n lineages alive on
    // the segment that ends there.  No lineage is born and none dies strictly
    // inside a segment, so on it both N and S = Sum_s tip_start_s are constant
    // and P is the straight line through (node.brts, node.pd) with slope N:
    //
    //   P(t) = N*t - S = node.pd + node.n * (t - node.brts)
    //
    // This is the one definition of P the model has.  The scorer reads it at
    // the nodes (node.pd itself); the sampler, which needs it at candidate
    // times between them, reads it here, so the two cannot disagree.  Exact,
    // and O(1) rather than a scan over the event list.
    static double pendant_pd(const node_t& node, double t) {
      return node.pd + node.n * (t - node.brts);
    }

    // The pendant PD at t, off the node that governs the segment t falls in.
    //
    // lower_bound_node returns the first node with brts >= t, which is the node
    // whose segment (brts_{i-1}, brts_i] contains t: its (pd, n, brts) are the
    // ones in force there.  Verified against the topology in
    // tests/testthat/test-d-compensator.R.
    double pendant_pd_at(double t, const tree_t& tree) const {
      auto it = lower_bound_node(t,
                                 tree.size(),
                                 reinterpret_cast<const node_t*>(tree.data()));
      return pendant_pd(*it, t);
    }

    // The time at which the segment ending at `it` begins: the previous node's
    // time, or 0 for the first segment.
    static double segment_start(const node_t* it, const node_t* first) {
      return (it == first) ? 0.0 : (it - 1)->brts;
    }

    // The two rates the proposal holds on the segment (t0, node.brts].
    //
    // The proposal reads only the alive set: the diversity N = node.n and the
    // mean pendant age M = P/N, with P extrapolated to the time the segment
    // begins.  It never reads D.  That is what makes the sampler and the
    // scorer read the same numbers: a birth splits the segment it falls in, so
    // the sampler (which sees the segment whole) and the scorer (which sees
    // the two halves) disagree about which node ends it, but they agree about
    // the alive set on it and about where it begins.
    //
    // Under the linear link Sum_s D_s = 0 over the alive lineages, so
    // N * lambda(N, M) is the model's own total speciation intensity whenever
    // no lineage's rate is clipped at zero.
    void proposal_rates(const param_t& pars, const node_t& node, double t0,
                        double& lambda, double& mu) const {
      node_t cur = node;
      cur.pd = pendant_pd(node, t0);
      lambda = speciation_rate(pars, cur);
      mu = std::max(extinction_rate(pars, cur), 1e-10);
    }

    // The intensity of the thinning sampler's birth process:
    //
    //   nh(t) = N * lambda_seg * (1 - rho * exp(-mu_seg * (T - t)))
    //
    // lambda_seg and mu_seg are the segment's rates (proposal_rates) and the
    // survival factor is the probability that a lineage born at t leaves no
    // sampled descendant, under a lifetime that is exponential with rate
    // mu_seg.  Only that factor carries t, and it falls with t, so nh is
    // decreasing on every segment and its value where the segment begins
    // dominates it there — for every model and every link.
    double nh_rate(double t, const param_t& pars, const tree_t& tree) const {
      const node_t* first = reinterpret_cast<const node_t*>(tree.data());
      auto it = lower_bound_node(t, tree.size(), first);
      double lambda = 0.0, mu = 0.0;
      proposal_rates(pars, *it, segment_start(it, first), lambda, mu);
      const double T = tree.back().brts;
      return lambda * it->n * (1.0 - rho_ * std::exp(-mu * (T - t)));
    }


    // log q(z | obs, theta): the log density of the augmentation the thinning
    // sampler drew (src/augment_tree.cpp, do_augment_tree_cont).
    //
    // The sampler draws, in forward time,
    //
    //   birth times  an inhomogeneous Poisson process of intensity nh(t),
    //                realised by thinning a homogeneous process of rate
    //                nh(segment start) — a dominating envelope, since nh falls
    //                within a segment — so the times have density
    //                prod_k nh(t_k) * exp(-int_0^T nh(t) dt);
    //   a lifetime   for the lineage born at t_k: unsampled extant with
    //                probability (1-rho)e^{-mu r}/(1-rho e^{-mu r}), otherwise
    //                exponential with rate mu_seg truncated to (0, r),
    //                r = T - t_k;
    //   a parent     uniformly over the labelled attachments, 2 per observed
    //                lineage alive (including the two crown lineages) and 1 per
    //                missing one: 2*tips + Ne of them (H5).
    //
    // The survival factor cancels against the lifetime's normalisation, which
    // is why each event contributes log(N * lambda * mu) - mu * lifespan (extinct)
    // or log(N * lambda * (1-rho)) - mu * (T - t) (unsampled extant).
    //
    // N, lambda_seg and mu_seg are the same numbers the sampler used, so the
    // compensator is int nh dt exactly: on a segment nh is N*lambda_seg times a
    // single exponential in t, whose integral is the closed form below for
    // every link.
    double sampling_prob(const param_t& pars, const tree_t& tree) const {
      double inte = 0;
      double logg = 0;
      double prev_brts = 0;
      double tips = tree[0].n;
      double Ne = 0.0;
      const double T = tree.back().brts;

      for (unsigned i = 0; i < tree.size(); ++i) {
        const auto& node = tree[i];
        double lambda = 0.0, mu = 0.0;
        proposal_rates(pars, node, prev_brts, lambda, mu);
        {
          const double dt = node.brts - prev_brts;
          const double int_segment = dt - rho_ * (1.0/mu) * (std::exp(-mu*(T - node.brts)) - std::exp(-mu*(T - prev_brts)));
          inte += node.n * lambda * int_segment;
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
    // Substituting u = (A + b*t - 1)/sqrt(2), dt = sqrt(2)/b du, gives
    //   = intercept * sqrt(pi/2) / b * [erf(z2) - erf(z1)]
    // where z_i = (A + b*t_i - 1) / sqrt(2).
    // NOTE: divide by the SIGNED b, not |b|.  For b < 0 we have z2 < z1, so
    // erf(z2) - erf(z1) < 0, and dividing by the negative b restores the
    // correct positive integral (the integrand is strictly positive).  Using
    // |b| would flip the sign of the hazard integral for negative beta_E.
    // Handles b→0 limit: intercept * exp(-(A-1)^2/2) * (t2 - t1)
    // Helper: integral of max(0, c + b*t) from t1 to t2 — one lineage's
    // contribution to the compensator under the linear link, where the rate is
    // a line clipped at zero.
    //
    // b == 0 leaves a constant rate max(0, c).  Otherwise f(t) = c + b*t has a
    // single root at t* = -c/b and is positive on one side of it: above t* when
    // b > 0, below it when b < 0.  The integral is the area of f over the part
    // of [t1, t2] on that side, which is empty when t* lies past the far
    // endpoint and the whole interval when it lies past the near one.  Over the
    // surviving [lo, hi] the integrand is a line, so the exact area is its
    // width times its value at the midpoint.
    static double relu_integral(double c, double b, double t1, double t2) {
      if (t2 <= t1) return 0.0;
      if (b == 0.0) return std::max(0.0, c) * (t2 - t1);
      const double root = -c / b;
      const double lo = (b > 0.0) ? std::max(t1, root) : t1;
      const double hi = (b > 0.0) ? t2 : std::min(t2, root);
      if (hi <= lo) return 0.0;
      return (hi - lo) * (c + b * (0.5 * (lo + hi)));
    }

    static double gauss_integral(double intercept, double A, double b, double t1, double t2) {
      if (std::abs(b) < 1e-12) {
        double d = A - 1.0;
        return intercept * std::exp(-0.5 * d * d) * (t2 - t1);
      }
      const double inv_sqrt2 = 1.0 / std::sqrt(2.0);
      double z1 = (A + b * t1 - 1.0) * inv_sqrt2;
      double z2 = (A + b * t2 - 1.0) * inv_sqrt2;
      return intercept * std::sqrt(M_PI / 2.0) / b * (std::erf(z2) - std::erf(z1));
    }

    double loglik(const param_t& pars, const tree_t& tree) const {
      const bool ep_exp = model_bin_[2] && (link_ == LinkType::exponential);
      const bool ep_gauss = model_bin_[2] && (link_ == LinkType::gaussian);
      const bool ep_linear = model_bin_[2] && (link_ == LinkType::linear);

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
          // Orthogonal basis {N, M=P/N, D=E-M}. D is constant within the
          // segment, so eta = beta_0 + beta_N*N + beta_M*M + beta_D*(E - M)
          //  = [beta_0 + beta_N*N + (beta_M - beta_D)*M] + beta_D*(t - ts_s).
          // The bracketed constant is A_lam; the running sums over exp(-beta_D*ts)
          // are unchanged (beta_D = pars[3], same slot).
          double N_seg = node.n;
          double M_seg = (N_seg > 0.0) ? node.pd / N_seg : 0.0;
          double A_lam = pars[0] + pars[1]*N_seg + (pars[2] - pars[3])*M_seg;
          double A_mu  = pars[4] + pars[5]*N_seg + (pars[6] - pars[7])*M_seg;
          inte += sum_exp_bE * exp_integral(A_lam, pars[3], prev_brts, node.brts)
                + sum_exp_gE * exp_integral(A_mu,  pars[7], prev_brts, node.brts);
        } else if (ep_gauss && dt > 0.0) {
          // EP + Gaussian: exact per-lineage integral over the segment via erf.
          // Orthogonal basis {N, M=P/N, D=E-M}; N and M are held constant per
          // segment (M at the node value), and only D = (t - ts_s) - M carries an
          // explicit time dependence.  For each lineage s alive throughout the
          // segment (born at tip_start ts_s), the covariate argument is
          //   eta_cov = beta_N*N + beta_M*M + beta_D*((t - ts_s) - M)
          //           = [beta_N*N + (beta_M - beta_D)*M - beta_D*ts_s] + beta_D*t
          // so A_lam,s = beta_N*N + (beta_M - beta_D)*M - beta_D*ts_s and b = beta_D.
          // The Gaussian rate does not factor into a running sum (unlike the
          // exponential link), so this is an O(alive) inner loop per segment.
          const double N_seg  = node.n;
          const double M_seg  = (N_seg > 0.0) ? node.pd / N_seg : 0.0;
          const double base_lam = pars[1] * N_seg + (pars[2] - pars[3]) * M_seg;
          const double base_mu  = pars[5] * N_seg + (pars[6] - pars[7]) * M_seg;
          for (const auto& s : tree) {
            if (is_extinction(s)) continue;
            // lineage s alive throughout (prev_brts, node.brts)?
            if (s.brts <= prev_brts && s.t_ext >= node.brts) {
              const double A_lam = base_lam - pars[3] * s.tip_start;
              const double A_mu  = base_mu  - pars[7] * s.tip_start;
              inte += gauss_integral(pars[0], A_lam, pars[3], prev_brts, node.brts);
              inte += gauss_integral(pars[4], A_mu,  pars[7], prev_brts, node.brts);
            }
          }
        } else if (ep_linear && dt > 0.0) {
          // EP + linear: exact per-lineage integral over the segment.
          //
          // Same decomposition as the exponential and gaussian branches.  N and
          // M are the segment's values; only D = (t - ts_s) - M carries time
          // within the segment, so for a lineage born at ts_s
          //   eta_s(t) = beta_0 + beta_N*N + beta_M*M + beta_D*((t - ts_s) - M)
          //            = [beta_0 + beta_N*N + (beta_M - beta_D)*M - beta_D*ts_s]
          //              + beta_D*t
          // and the linear link makes the rate max(0, eta_s(t)): a line clipped
          // at zero, whose kink falls inside the segment whenever the rate
          // crosses zero there.  relu_integral does that case split exactly.
          //
          // The alive set is the exponential branch's: the two crown lineages
          // (tip_start 0, alive throughout), plus every lineage whose birth
          // node lies at or before the segment start and which has not died
          // before the segment ends.  That is node.n lineages, which is what
          // makes the branch collapse onto dt*n*(lambda+mu) when beta_D and
          // gamma_D are zero.
          const double N_seg = node.n;
          const double M_seg = (N_seg > 0.0) ? node.pd / N_seg : 0.0;
          const double base_lam = pars[0] + pars[1] * N_seg + (pars[2] - pars[3]) * M_seg;
          const double base_mu  = pars[4] + pars[5] * N_seg + (pars[6] - pars[7]) * M_seg;
          double seg = 0.0;
          auto add_lineage = [&](double ts) {
            seg += relu_integral(base_lam - pars[3] * ts, pars[3], prev_brts, node.brts);
            seg += relu_integral(base_mu  - pars[7] * ts, pars[7], prev_brts, node.brts);
          };
          add_lineage(0.0);           // the two crown lineages
          add_lineage(0.0);
          const unsigned last = static_cast<unsigned>(tree.size()) - 1u;
          for (unsigned j = 0; j < tree.size(); ++j) {
            const auto& s = tree[j];
            if (is_extinction(s) || j == last) continue;   // the last node marks the present
            if (s.brts <= prev_brts && s.t_ext >= node.brts) add_lineage(s.tip_start);
          }
          inte += seg;
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

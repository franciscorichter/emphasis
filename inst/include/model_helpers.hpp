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

// some model utilities
// Hanno 2020

#ifndef EMPHASIS_MODEL_HELPRES_HPP_INCLUDED
#define EMPHASIS_MODEL_HELPRES_HPP_INCLUDED

#include <limits>
#include <cmath>
#include <random>
#include <array>
#include <chrono>
#include <thread>
#include <numeric>
#include <algorithm>
#include "model.hpp"


#ifndef EMPHASIS_LOGSUM_LOWER_TRESHOLD
#define EMPHASIS_LOGSUM_LOWER_TRESHOLD 10e-20
#endif

#ifndef EMPHASIS_LOGSUM_UPPER_TRESHOLD
#define EMPHASIS_LOGSUM_UPPER_TRESHOLD 10e+20
#endif

#define t_ext_tip 10e10         /* t_ext for present nodes (observed tips) */
#define t_ext_extinct 0.0       /* t_ext for extinction nodes */
#define t_ext_unsampled 5e10    /* t_ext for unsampled extant species (rho < 1) */


namespace emphasis {

  // A tip_start is a forward time and is therefore never negative; this value
  // in focal_tip_start marks an event whose splitting lineage is not known
  // (a branching-time vector carrying no topology, or an augmented lineage
  // drawn when no parent was on record).  Model::e_s falls back to the
  // mean pendant age M there, which sets D = 0 for that event.
  constexpr double ts_unknown = -1.0;

  // node_t::clade on an observed node: 1 when the node carries the tip_start
  // of the lineage that splits at it (the caller supplied the topology),
  // 0 otherwise.
  constexpr int clade_topology = 1;

  /* tree node */
  struct node_t
  {
    double brts;
    double n;             /* n[i] = number of species in [time_i-1, time_i) */
    double t_ext;         /* emp_t_ext_tip for present-day species;  emp_t_ext_extinct for extinction nodes */
    double pd;            // pendant PD P(brts) = sum over the lineages alive on
                          // the segment ending here of (brts - their tip_start)
    double tip_start;     // forward time at which the lineage BORN at this node
                          // last became a pendant tip (= brts for a speciation
                          // node, the dying lineage's tip_start for an
                          // extinction node)
    double focal_tip_start; // tip_start of the lineage whose event this is (the
                          // lineage that splits, or the one that dies);
                          // ts_unknown when no parent is on record
    int clade;            // clade_topology on an observed node whose
                          // focal_tip_start came from the observed topology
    int id;               // unique stable lineage ID (assigned at creation); -1 = unset
    int parent_id;        // id of parent lineage; -1 = root / initial tree
  };
  
  
  
  namespace detail {

    static constexpr double huge = std::numeric_limits<double>::max();


    // returns low-entropy 512 bit array for seed sequence
    // based on std::chrono::high_resolution_clock.
    // ripped from rndutils
    inline auto make_low_entropy_seed_array() noexcept->std::array<uint64_t, 8>
    {
      // the classic: time, advertised with nano-second resolution.
      const auto e1 = static_cast<uint64_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
      // different between invocations from different threads within one app: thread-id
      const auto tid = std::this_thread::get_id();
      const uint64_t e2{ std::hash<typename std::remove_const<decltype(tid)>::type>()(tid) };
      return std::array<uint64_t, 8>{ {
          e1, e2,
          0x000000003c10b019, 0x2bf820b4dd7c1a8a,
          0x9901cf90a40883da, 0x5a3686b2e1de6e51,
          0x000000cc0494d228, 0x000000cc04b66740
      }};
    }


    // random number generator from low-entropy seed sequence
    // ripped from rndutils
    template <typename URNG>
    inline auto make_random_engine() -> URNG
    {
      auto seed_array = make_low_entropy_seed_array();
      std::seed_seq sseq(seed_array.cbegin(), seed_array.cend());
      return URNG(sseq);
    }


    inline bool is_extinction(const node_t& node) { return node.t_ext == t_ext_extinct; }
    inline bool is_tip(const node_t& node) { return node.t_ext == t_ext_tip; }
    inline bool is_unsampled(const node_t& node) { return node.t_ext == t_ext_unsampled; }
    inline bool is_missing(const node_t& node) { return !(is_extinction(node) || is_tip(node) || is_unsampled(node)); }


    struct node_less
    {
      bool operator()(const node_t& a, const node_t& b) const noexcept { return a.brts < b.brts; };
      bool operator()(const node_t& a, double val) const noexcept { return a.brts < val; };
      bool operator()(double val, const node_t& a) const noexcept { return val < a.brts; };
    };


    inline const node_t* lower_bound_node(double t, unsigned n, const node_t* tree)
    {
      auto it = std::lower_bound(tree, tree + n, t, node_less{});
      return std::min(it, tree + n - 1);
    }


    // Running sum of log(val) over positive rates, accumulated as a product
    // that is folded into sum_ whenever it leaves (LOWER, UPPER).
    // A rate <= 0 has log-density -Inf regardless of the other terms, so it is
    // recorded as a flag and never enters prod_ or sum_: result() then returns
    // -Inf and the sign of the result is never read off sum_.
    class log_sum
    {
    public:
      double result() const 
      { 
        if (zero_) return -std::numeric_limits<double>::infinity();
        return std::log(prod_) + sum_;
      }

      void operator+=(double val)
      {
        if (val <= 0.0) {
          zero_ = true;
          return;
        }
        if ((prod_ > EMPHASIS_LOGSUM_LOWER_TRESHOLD) && (prod_ < EMPHASIS_LOGSUM_UPPER_TRESHOLD)) {
          prod_ *= val;
        }
        else {
          sum_ += std::log(prod_) + std::log(val);
          prod_ = 1;
        }
      }

    private:
      double prod_ = 1;
      double sum_ = 0;
      bool zero_ = false;
    };


    template <typename RENG>
    inline double trunc_exp(double upper, double rate, RENG& reng)
    {
      std::exponential_distribution<double> exp_dist(rate);
      double result = exp_dist(reng);
      while (result > upper) {
        result = exp_dist(reng);
      }
      return result;
    }


    // Int_0_t1 (1-exp(-mu*(tm-t)))
    class mu_integral
    {
    public:
      mu_integral(double mu, double tm)
      : mu_(mu),
        s_(1.0 / (mu * std::exp(mu * tm)))
      {}

      double operator()(double t0, double t1)
      {
        const double expt0 = (pt1_ == t0) ? expt1_ : std::exp(mu_ * t0);
        expt1_ = std::exp(mu_ * t1);
        return (t1 - t0) - s_ * (expt1_ - expt0);
      }

    private:
      double pt1_ = -1.0;
      double expt1_ = 0;
      const double mu_ = 0;
      const double s_ = 0;
    };


    inline double calculate_pd(double tm, unsigned n, const node_t* tree)
    {
      double brts = 0.0;
      double prev_brts = 0;
      double ni = tree[0].n;
      double pd = 0.0;
      for (unsigned i = 0; (i < n) && (tree[i].brts <= tm); ++i) {
        const auto& node = tree[i];
        brts = node.brts;
        if ((node.t_ext > tm)) {
          pd += (brts - prev_brts) * ni++;
          prev_brts = brts;
        }
      }
      return pd + (tm - prev_brts) * ni;   // remainder
    }
    
    inline double calculate_pd2(double tm, const std::vector<node_t>& tree)
    {
      return calculate_pd(tm, static_cast<unsigned>(tree.size()), tree.data());
    }

    // Pendant PD: sum of pendant edge lengths of alive lineages at time tm.
    // Each alive lineage i contributes (tm - tip_start_i).
    //
    // One lineage per node, so the two crown lineages are not counted and a
    // tip_start reset by a later split is not seen.  This is the value the
    // thinning envelope reads at an arbitrary candidate time inside
    // Model::nh_rate, where no event-list sweep is available; the pd stored on
    // the nodes comes from compute_pendant_pd(), which keeps the running
    // (N, sum tip_start) state and is exact.
    inline double calculate_pendant_pd(double tm, const std::vector<node_t>& tree)
    {
      double ppd = 0.0;
      for (const auto& node : tree) {
        if (node.brts > tm) break;  // tree is sorted by brts
        // A lineage is "alive" at tm if it is not an extinction node and
        // its t_ext > tm (either a tip surviving to present or an extinction
        // that happens after tm).
        if (!is_extinction(node) && node.t_ext > tm) {
          ppd += (tm - node.tip_start);
        }
      }
      return ppd;
    }

  }

}

#endif

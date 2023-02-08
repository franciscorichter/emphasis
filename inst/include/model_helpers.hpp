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

#define t_ext_tip 10e10     /* t_ext for present nodes */
#define t_ext_extinct 0.0   /* t_ext for extinction nodes */


namespace emphasis {
  /* tree node */
  struct node_t
  {
    double brts;
    double n;         /* n[i] = number of species in [time_i-1, time_i) */
  double t_ext;     /* emp_t_ext_tip for present-day species;  emp_t_ext_extinct for extinction nodes */
  double pd;
  int clade; // required for sim_tree
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
    inline bool is_missing(const node_t& node) { return !(is_extinction(node) || is_tip(node)); }


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


    class log_sum
    {
    public:
      double result() const 
      { 
        double r = std::log(prod_) + sum_;
        if (!std::isfinite(r)) {
          const double s = std::signbit(sum_) ? -1.0 : 1.0;
          r = s * std::numeric_limits<double>::infinity();
        }
        return r; 
      }

      void operator+=(double val)
      {
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

  }

}

#endif

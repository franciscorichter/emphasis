// Part of the emphasis R package (MIT licence, see LICENSE). This file builds on
// the earlier C++ implementation by Hanno Hildenbrandt (2020).

// public C++ API
// Hanno 2020

#ifndef EMPHASIS_EMPHASIS_HPP_INCLUDED
#define EMPHASIS_EMPHASIS_HPP_INCLUDED

#include <chrono>
#include <stdexcept>
#include <string>
#include <limits>
#include <memory>
#include <vector>
#include <functional>
#include "model.hpp"

namespace emphasis {

  static constexpr int default_max_missing_branches = 10000;
  static constexpr double default_max_aug_lambda = 500.0;

  using brts_t = std::vector<double>;     // input tree


  // The wall-clock deadline of a fit, passed down from the R driver as the
  // seconds it has left.  Read-only once set, so every TBB worker may compare
  // the clock against it without synchronisation.
  using steady_clock_t = std::chrono::steady_clock;
  struct deadline_t
  {
    bool active = false;
    steady_clock_t::time_point at{};

    // `seconds_left` < 0 or not finite: no deadline.  0: already passed.
    static deadline_t in_seconds(double seconds_left)
    {
      deadline_t d;
      if (seconds_left >= 0.0 && seconds_left < 1e9) {
        d.active = true;
        d.at = steady_clock_t::now() +
          std::chrono::duration_cast<steady_clock_t::duration>(std::chrono::duration<double>(seconds_left));
      }
      return d;
    }
    bool passed() const { return active && steady_clock_t::now() >= at; }
  };


  struct E_step_info_t
  {
    double fhat = 0;                    // mean, unscaled, weight
    int num_trees = 0;                  // number of trees augmented
    int rejected_overruns = 0;          // trees rejected because overrun of missing branches
    int rejected_lambda = 0;            // trees rejected because of lambda overrun
    int rejected_zero_weights = 0;      // trees rejected because of zero-weight (log_w = -Inf)
    int rejected_nonfinite = 0;         // trees rejected because log_w is +Inf or NaN
    int rejected = 0;                   // trees rejected because of unhandled exception
    double elapsed = 0;                 // elapsed runtime [ms]
    double logf = 0;                    // likelihood of simulated tree
    double logg = 0;                    // sampling probability
    bool timed_out = false;             // the deadline passed before N trees were drawn
  };


  // results from e
  struct E_step_t
  {
    E_step_t(E_step_t&&) = default;
    E_step_t& operator=(E_step_t&&) = default;
    E_step_t(const E_step_t&) = delete;
    E_step_t& operator=(const E_step_t&) = delete;
    E_step_t() {};
    ~E_step_t() {};

    std::vector<tree_t> trees;          // augmented trees
    std::vector<double> weights;
    std::vector<double> logg_;
    std::vector<double> logf_;

    E_step_info_t info;
  };


  class emphasis_error : public std::runtime_error
  {
  public:
    explicit emphasis_error(const std::string& what) : std::runtime_error(what) {}
    explicit emphasis_error(const char* what) : std::runtime_error(what) {}
  };
  

  class emphasis_error_E : public std::runtime_error
  {
  private:
    static std::string make_msg(const E_step_t& E) {
      std::string msg = "maxN exceeded with rejection reasons: ";
      msg += std::to_string(E.info.rejected_lambda) + " lambda; ";
      msg += std::to_string(E.info.rejected_overruns) + " overruns; ";
      msg += std::to_string(E.info.rejected_zero_weights) + " zero weights; ";
      msg += std::to_string(E.info.rejected_nonfinite) + " non-finite weights; ";
      msg += std::to_string(E.info.rejected) + " unhandled exception; ";
      msg += " Trees so far: " + std::to_string(E.info.num_trees);
      return msg;
    }
    
  public:
    explicit emphasis_error_E(const E_step_t& E) : std::runtime_error(make_msg(E)), info(E.info) {}
    E_step_info_t info;
  };


  // parent_tip_start: tip_start of the lineage that splits at each observed
  // branching event, in the forward-time order of `brts`.  Empty when the
  // caller has only branching times and no topology.
  // parent_id: the id of that splitting lineage per observed event, in the
  // same order -- a crown id (crown_id_a / crown_id_b) or the index of the
  // earlier event whose daughter it is -- which is what the ED covariate
  // needs to know which lineage is which.  Empty when not supplied.
  E_step_t E_step(int N,      // sample size
                  int maxN,   // max number of augmented trees (incl. invalid)
                  const param_t& pars,
                  const brts_t& brts,
                  const Model& model,
                  int max_missing = default_max_missing_branches,
                  double max_lambda = default_max_aug_lambda,
                  int num_threads = 0,
                  double max_time_seconds = 0,  // 0 = no time limit
                  const std::vector<double>& parent_tip_start = {},
                  const std::vector<int>& parent_id = {},
                  const deadline_t& deadline = {});  // inactive: no deadline


  // results from m
  struct M_step_t
  {
    param_t estimates;
    int opt = -1;                       // nlopt result
    double minf = 0.0;
    double elapsed = 0.0;               // elapsed runtime [ms]
  };


  using conditional_fun_t = std::function<double(const param_t&)>;
  
  
  M_step_t M_step(const param_t& pars,
                  const std::vector<tree_t>& trees,          // augmented trees
                  const std::vector<double>& weights,
                  const Model& model,
                  const param_t& lower_bound = {}, // overrides model.lower_bound
                  const param_t& upper_bound = {}, // overrides model.upper.bound
                  double xtol_rel = 0.001,
                  int num_threads = 0,
                  conditional_fun_t* conditional = nullptr);


  // results from mcem
  struct mcem_t
  {
    mcem_t(mcem_t&&) = default;
    mcem_t& operator=(mcem_t&&) = default;
    mcem_t(const mcem_t&) = delete;
    mcem_t& operator=(const mcem_t&) = delete;
    mcem_t() {};
    ~mcem_t() {};

    E_step_t e;
    M_step_t m;
  };

  
  mcem_t mcem(int N,      // sample size
              int maxN,   // max. number of augmented trees (incl. invalid)
              const param_t& pars,
              const brts_t& brts,
              const Model& model,
              int max_missing = default_max_missing_branches,
              double max_lambda = default_max_aug_lambda,
              const param_t& lower_bound = {}, // overrides model.lower_bound
              const param_t& upper_bound = {}, // overrides model.upper.bound
              double xtol_rel = 0.001,
              int num_threads = 0,
              conditional_fun_t* conditional = nullptr,
              const std::vector<double>& parent_tip_start = {},
              const std::vector<int>& parent_id = {},
              const deadline_t& deadline = {});

}

#endif

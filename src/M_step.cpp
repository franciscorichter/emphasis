#include <thread>
#include <mutex>
#include <chrono>
#include <numeric>
#include <cmath>
#include <limits>
#include <algorithm>
#include <functional>
#include <memory>
#include <tbb/tbb.h>
#include "model.hpp"
#include "emphasis.hpp"
#include "sbplx.hpp"


namespace emphasis {

  namespace {

    struct nlopt_f_data
    {
      nlopt_f_data(const Model& M,
                   const std::vector<tree_t>& Trees,
                   const std::vector<double>& W,
                   conditional_fun_t* Conditional)
        : model(M), trees(Trees), w(W), conditional(Conditional),
          sum_w(std::accumulate(W.cbegin(), W.cend(), 0.0))
      {
      }

      ~nlopt_f_data()
      {
        auto empty = tree_t{};
      }

      const Model& model;
      const std::vector<tree_t>& trees;
      const std::vector<double>& w;
      conditional_fun_t* conditional;
      const double sum_w;   // sum of the weights over the trees passed
    };


    // Q(pars) = sum_i w_i * loglik(pars, tree_i), a weighted sum with
    // unnormalised weights (BDI: mean 1; thinning: max-scaled).
    // Terms with w_i == 0 are skipped, so a tree of zero density under the
    // current parameters (loglik = -Inf, weight 0) does not turn Q into NaN.
    // A non-finite loglik on a tree with w_i != 0 makes the objective +Inf:
    // nlopt treats the point as infeasible instead of stalling on NaN.
    double objective(unsigned int n, const double* x, double*, void* func_data)
    {
      auto psd = reinterpret_cast<nlopt_f_data*>(func_data);
      param_t pars(x, x + n);
      constexpr double inf = std::numeric_limits<double>::infinity();
      const double Q = tbb::parallel_reduce(tbb::blocked_range<size_t>(0, psd->trees.size()), 0.0,
        [&](const tbb::blocked_range<size_t>& r, double q) -> double {
          for (size_t i = r.begin(); i < r.end(); ++i) {
            if (psd->w[i] == 0.0) continue;
            const double loglik = psd->model.loglik(pars, psd->trees[i]);
            if (!std::isfinite(loglik)) return -inf;   // partial sums are finite or -Inf
            q += loglik * psd->w[i];
          }
          return q;
        },
        std::plus<double>{}
      );
      if (!std::isfinite(Q)) {
        return inf;
      }
      if (nullptr == psd->conditional) {
        return -Q;
      }
      // conditional(pars) returns log(P_tree(pars)).
      // Conditioned MLE maximizes Q/sum_w - log(P_tree); multiplied by sum_w
      // that is Q - sum_w * log(P_tree), so we minimize
      // -Q + sum_w * log(P_tree). The argmin is invariant to a positive
      // rescaling of the weights.
      return -Q + psd->sum_w * psd->conditional->operator()(pars);
    }
  }


  M_step_t M_step(const param_t& pars,
                  const std::vector<tree_t>& trees,          // augmented trees
                  const std::vector<double>& weights,
                  const Model& model,
                  const param_t& lower_bound, // overrides model.lower_bound
                  const param_t& upper_bound, // overrides model.upper.bound
                  double xtol_rel,
                  int num_threads,
                  conditional_fun_t* conditional)
  {
 //   if (!model.is_threadsafe()) num_threads = 1;
    tbb::task_scheduler_init _tbb((num_threads > 0) ? num_threads : tbb::task_scheduler_init::automatic);
    auto T0 = std::chrono::high_resolution_clock::now();
    nlopt_f_data sd{ model, trees, weights, conditional };
    auto M = M_step_t{};
    sbplx nlopt(pars.size());
    M.estimates = pars;
    nlopt.set_xtol_rel(xtol_rel);
    auto lower = lower_bound.empty() ? model.lower_bound() : lower_bound;
    if (!lower.empty()) nlopt.set_lower_bounds(lower);
    auto upper = upper_bound.empty() ? model.upper_bound() : upper_bound;
    if (!upper.empty()) nlopt.set_upper_bounds(upper);
    nlopt.set_min_objective(objective, &sd);
    M.minf = nlopt.optimize(M.estimates);
    M.opt = static_cast<int>(nlopt.result());
    auto T1 = std::chrono::high_resolution_clock::now();
    M.elapsed = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(T1 - T0).count());
    return M;
  }

}

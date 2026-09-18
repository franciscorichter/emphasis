#include <thread>
#include <mutex>
#include <chrono>
#include <numeric>
#include <cmath>
#include <limits>
#include <algorithm>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
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
    // nlopt treats the point as infeasible instead of stalling on NaN.  The
    // conditional's return value is guarded the same way: it is an arbitrary R
    // function (mcem() takes one from the caller, and predict.gam returns NA on
    // a newdata row it cannot evaluate), and a NaN there leaves the objective
    // NaN at every point, on which sbplx returns the initial parameters with
    // status XTOL_REACHED and the driver reads the iteration as converged.
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
            // Marks the whole subrange infeasible, discarding the q it had
            // accumulated.  This is not a claim that every partial sum is
            // finite or -Inf: a subrange of finite terms can still overflow to
            // +Inf, and -Inf + (+Inf) is NaN.  The outer !isfinite(Q) test
            // below is what turns any of those into +Inf; do not remove it.
            if (!std::isfinite(loglik)) return -inf;
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
      const double logp = psd->conditional->operator()(pars);
      if (!std::isfinite(logp)) {
        return inf;
      }
      return -Q + psd->sum_w * logp;
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
    // A weight multiplies a log-density inside Q and scales the conditional
    // through sum_w.  A NaN weight makes the objective the constant +Inf, on
    // which sbplx returns a moved point with status XTOL_REACHED; a negative
    // weight silently minimizes the density of that tree.  Both are rejected
    // here, next to the length checks m_cpp already performs.  Zero is allowed:
    // it is how the drivers drop a zero-density draw.
    for (size_t i = 0; i < weights.size(); ++i) {
      if (!std::isfinite(weights[i]) || weights[i] < 0.0) {
        throw std::invalid_argument("M_step: weight " + std::to_string(i + 1) + " is " +
          std::to_string(weights[i]) + "; weights must be finite and non-negative");
      }
    }
    // tbb::task_scheduler_init was removed in oneTBB; a task_arena bounds the
    // parallelism in both the legacy and the current library.  objective()
    // runs a parallel_reduce, so the optimiser call is what must execute
    // inside the arena: a parallel algorithm started outside one runs on the
    // default arena, which is sized by hardware concurrency and ignores
    // num_threads.
    const int nth = (num_threads > 0)
      ? std::min(num_threads, static_cast<int>(std::thread::hardware_concurrency()))
      : static_cast<int>(std::thread::hardware_concurrency());
    tbb::task_arena arena(std::max(1, nth));
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
    // The first simplex's step, per coordinate: a tenth of the box where the
    // box is finite and not the +-1e6 placeholder, NLopt's own default
    // otherwise (a quarter of |x|, or 1 at x = 0).  NLopt's default under
    // finite bounds is 0.75 of the distance to the nearer bound, which for a
    // start a hair inside a bound -- the warm start of a nested fit clipped
    // to the box lands there -- is that hair: sbplx then satisfies its
    // x-tolerance on the spot and returns the start unchanged, at every EM
    // iteration, as XTOL_REACHED.  Measured on the ED simulation arm
    // (2026-09-18): a start 1.7e-18 above the lower bound of beta_N froze
    // every M-step of the fit; the same start 1e-6 inside moved.
    if (!lower.empty() && !upper.empty() && lower.size() == pars.size() && upper.size() == pars.size()) {
      std::vector<double> step(pars.size());
      for (size_t i = 0; i < pars.size(); ++i) {
        const double width = upper[i] - lower[i];
        const bool boxed = std::isfinite(width) && width > 0.0 && width <= 1e4;
        step[i] = boxed ? 0.1 * width : (pars[i] != 0.0 ? 0.25 * std::abs(pars[i]) : 1.0);
      }
      nlopt.set_initial_step(step);
    }
    nlopt.set_min_objective(objective, &sd);
    arena.execute([&] {
      M.minf = nlopt.optimize(M.estimates);
      M.opt = static_cast<int>(nlopt.result());
    });
    auto T1 = std::chrono::high_resolution_clock::now();
    M.elapsed = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(T1 - T0).count());
    return M;
  }

}

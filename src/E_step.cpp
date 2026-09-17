#include <mutex>
#include <atomic>
#include <algorithm>
#include <numeric>
#include <limits>
#include <string>
#include <tbb/tbb.h>
#include "model.hpp"
#include "emphasis.hpp"
#include "augment_tree.hpp"
#include "model_helpers.hpp"
#include "precision_weights.hpp"


namespace emphasis {
  namespace detail {
    // this little addition reduces the load to memory allocator massively.
    tree_t thread_local pooled_tree;

    void inplace_cumsum_of_diff(brts_t& input)
    {
      double sum = 0.0;
      for (size_t i = 1; i < input.size(); ++i) {
        input[i - 1] = sum += (input[i - 1] - input[i]);
      }
      input.back() += sum;
    }


    // Observed tree from its branching times, and — when the caller has the
    // topology — from the tip_start of the lineage that splits at each
    // observed branching event.
    //
    // After inplace_cumsum_of_diff the entries of `brts` are forward times:
    // entry i < size-1 is the observed branching event i, and the last entry
    // is the present, a terminal marker rather than an event.  The two crown
    // lineages have no node of their own.
    //
    // `parent_tip_start` is in that same order and may be
    //   empty                : no topology; every observed event records an
    //                          unknown splitting lineage (tip_start 0,
    //                          focal_tip_start ts_unknown), which is what a
    //                          bare branching-time vector can support.
    //   brts.size() - 1      : one entry per observed branching event.
    //   brts.size()          : one entry per node; the last is ignored.
    tree_t create_tree(brts_t brts, const std::vector<double>& parent_tip_start,
                       const std::vector<int>& parent_id)
    {
      inplace_cumsum_of_diff(brts);
      const size_t m = brts.size();
      const size_t np = parent_tip_start.size();
      if (np != 0 && np != m && np + 1 != m) {
        throw emphasis_error("create_tree: parent_tip_start must be empty, or hold one "
                             "entry per observed branching event");
      }
      const size_t nid = parent_id.size();
      if (nid != 0 && nid != m && nid + 1 != m) {
        throw emphasis_error("create_tree: parent_id must be empty, or hold one "
                             "entry per observed branching event");
      }
      const bool topology = (np != 0);
      tree_t tree;
      tree.reserve(m);
      for (size_t i = 0; i < m; ++i) {
        node_t node{};
        node.brts = brts[i];
        node.n = 2.0 + i;
        node.t_ext = t_ext_tip;
        node.pd = 0.0;
        // The lineage born at this node became a pendant tip here; without a
        // topology the old convention (every observed lineage dates from the
        // crown) is kept, so that a bare branching-time vector reproduces the
        // values it produced before.
        node.tip_start = topology ? brts[i] : 0.0;
        node.focal_tip_start = (topology && (i + 1 < m)) ? parent_tip_start[i] : ts_unknown;
        node.clade = topology ? clade_topology : 0;
        node.id = -1;
        // The lineage this event splits, when the caller named it; the
        // augmentation assigns node ids in this same forward-time order, so
        // the caller's indices and the ids agree.
        node.parent_id = (nid != 0 && i < nid) ? parent_id[i] : -1;
        tree.push_back(node);
      }
      std::sort(tree.begin(), tree.end(), detail::node_less{});
      return(tree);
    }

  }


  E_step_t E_step(int N,
                  int maxN,
                  const param_t& pars,
                  const brts_t& brts,
                  const Model& model,
                  int max_missing,
                  double max_lambda,
                  int num_threads,
                  double max_time_seconds,  // 0 = default 120s safety limit
                  const std::vector<double>& parent_tip_start,
                  const std::vector<int>& parent_id)
  {
    // N < 1 leaves the sample empty on every path: nothing is stored, the
    // `num_trees < N` test below is false, and max_element/fhat would then
    // read an empty weight vector and divide by zero.
    if (N < 1) {
      throw emphasis_error("E_step: sample size must be at least 1");
    }
    if (!model.is_threadsafe()) num_threads = 1;
    num_threads = std::max(1, std::min(num_threads, static_cast<int>(std::thread::hardware_concurrency())));
    std::mutex mutex;
    std::atomic<bool> stop{ false };

    tree_t init_tree = detail::create_tree(brts, parent_tip_start, parent_id);
    auto E = E_step_t{};
    auto T0 = std::chrono::high_resolution_clock::now();

    // Integer division gives 0 when maxN < num_threads, which is not a valid
    // grainsize.
    const int grainsize = std::max(1, maxN / std::max(1, num_threads));
    tbb::task_arena arena(std::max(1, num_threads));

    // `stop` is written only under `mutex` and only ever set to true, so an
    // attempt that observes it under the lock is either the one that filled
    // the sample or one that finished after the sample was full; the latter
    // is neither stored nor counted.  Every counter therefore refers to an
    // attempt completed before tree N was pushed, and
    // trees.size() + rejected_zero_weights is the number of completed
    // augmentations that fhat averages over.  Attempts still in flight when
    // the N-th tree is pushed are dropped from numerator and denominator
    // alike; the drop does not depend on their outcome.
    // A parallel algorithm started outside the arena runs on the default one,
    // which is sized by hardware concurrency: num_threads would have no effect.
    arena.execute([&] {
    tbb::parallel_for(tbb::blocked_range<unsigned>(0, maxN, grainsize), [&](const tbb::blocked_range<unsigned>& r) {
      for (unsigned i = r.begin(); i < r.end(); ++i) {
        if (stop) break;
        try {
          // Time budget: prevent infinite spinning when IS collapse
          // causes nearly all augmentations to fail.
          // Default 120s safety limit; set max_time_seconds > 0 to override.
          {
            double limit = (max_time_seconds > 0) ? max_time_seconds : 120.0;
            auto now = std::chrono::high_resolution_clock::now();
            double elapsed = std::chrono::duration<double>(now - T0).count();
            if (elapsed > limit) {
              std::lock_guard<std::mutex> _(mutex);
              stop = true;
              break;
            }
          }
          // reuse tree from pool
          auto& pool_tree = detail::pooled_tree;

          emphasis::augment_tree(pars, init_tree, model, max_missing, max_lambda, pool_tree);

          double log_w = 0.0;
          double logf = 0.0;
          double logg = 0.0;
          {
            logf = model.loglik(pars, pool_tree);
            logg = model.sampling_prob(pars, pool_tree);
            log_w = logf - logg;
          }
          std::lock_guard<std::mutex> _(mutex);
          if (stop) break;
          if (std::isfinite(log_w)) {
            // any finite log-weight is accepted; the magnitude is handled by
            // the max-shift in calc_sum_w and downstream.
            if (static_cast<int>(E.trees.size()) < N) {
              E.trees.emplace_back(pool_tree.cbegin(), pool_tree.cend());
              E.weights.push_back(log_w);
              E.logf_.push_back(logf);
              E.logg_.push_back(logg);
              if (static_cast<int>(E.trees.size()) == N) {
                stop = true;
              }
            }
          }
          else if (logf == -std::numeric_limits<double>::infinity()) {
            // completed augmentation with f = 0 (e.g. lambda = 0 at a
            // speciation node): zero weight, enters the fhat denominator.
            // Classification is on logf, not on log_w: when logg has also
            // underflowed to -Inf the difference is NaN, and the draw is
            // still a completed augmentation of probability zero.
            ++E.info.rejected_zero_weights;
          }
          else {
            // logf = +Inf, or a NaN log_w whose logf is finite (g = 0 or
            // g = Inf): not a valid draw.
            ++E.info.rejected_nonfinite;
          }
        }
        catch (const augmentation_overrun&) {
          std::lock_guard<std::mutex> _(mutex);
          if (stop) break;
          ++E.info.rejected_overruns;
        }
        catch (const augmentation_lambda&) {
          std::lock_guard<std::mutex> _(mutex);
          if (stop) break;
          ++E.info.rejected_lambda;
        }
	      catch (...) {
          std::lock_guard<std::mutex> _(mutex);
          if (stop) break;
          ++E.info.rejected;
	      }
      }
    });
    });
    E.info.num_trees = static_cast<int>(E.trees.size());
    auto T1 = std::chrono::high_resolution_clock::now();
    E.info.elapsed = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(T1 - T0).count());

    // An empty sample is thrown on unconditionally: max_element below has no
    // defined result on an empty range, and S_completed would be zero.
    if (E.info.num_trees < N || E.weights.empty()) {
      throw emphasis_error_E(E);
    }
  
    const double max_log_w = *std::max_element(E.weights.cbegin(), E.weights.cend());
    // NOTE: calc_sum_w has a deliberate side-effect — it OVERWRITES each
    // E.weights[i] in place, converting it from a LOG weight (log w_i, as
    // stored during augmentation) to a LINEAR weight exp(log w_i - max_log_w).
    // After this call E.weights holds linear (max-scaled) weights, which is
    // exactly what the M-step objective expects. The M-step maximizer is
    // invariant to the omitted exp(-max_log_w) and normalization factors.
    double sum_w = calc_sum_w(E.weights.begin(), E.weights.end(), max_log_w);
    // Denominator: trees summed plus zero-weight trees (completed
    // augmentations with w=0 due to lambda=0).  Overrun/lambda/non-finite
    // rejections are excluded: those are computational failures independent
    // of theta (see tech report).
    int S_completed = static_cast<int>(E.trees.size()) + E.info.rejected_zero_weights;
    E.info.fhat = std::log(sum_w / S_completed) + max_log_w;
    return E;
  }


}

#include <cassert>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <atomic>
#include <tuple>
#include <memory>
#include <thread>
#include <unordered_map>
#include <set>
#include <tbb/tbb.h>
#include "model.hpp"
#include "augment_tree.hpp"
#include "model_helpers.hpp"


namespace emphasis {

  namespace {


    template <typename IT>
    inline IT make_node(IT it, double t, double n, double t_ext, int id = -1, int parent_id = -1)
    {
      it->brts = t; it->n = n; it->t_ext = t_ext; it->pd = 0.0;
      it->tip_start = t; it->focal_tip_start = 0.0;
      it->clade = 0; it->id = id; it->parent_id = parent_id;
      return it;
    }


    template <typename IT>
    inline IT make_extinct_node(IT it, double t, double n, double t_spec, int id = -1, int parent_id = -1)
    {
      it->brts = t; it->n = n; it->t_ext = t_ext_extinct; it->pd = 0.0;
      it->tip_start = t_spec; it->focal_tip_start = t_spec;
      it->clade = 0; it->id = id; it->parent_id = parent_id;
      return it;
    }


    auto thread_local reng = detail::make_random_engine<std::default_random_engine>();


    // Pendant PD and per-lineage tip_start by one forward sweep over the event
    // list — the bookkeeping the simulator itself runs
    // (inst/include/general_tree.hpp):
    //
    //   P(t) = N*t - S,   S = sum of tip_start over the lineages alive at t
    //   start          : two crown lineages, tip_start 0, so S = 0 and N = 2
    //   p speciates    : S += 2*t - ts[p];  ts[p] = ts[c] = t;  ++N
    //   s goes extinct : S -= ts[s];  --N
    //
    // node.pd is P at the node's own time over the N lineages alive on the
    // segment that ends there, which is the count node.n carries.  One pass
    // with a log-time lookup per event; the value it replaces summed
    // (t - tip_start) over the whole tree at every node, quadratically.
    //
    // The tip_start of the splitting lineage is known for an observed event
    // when the caller supplied the topology (create_tree stored it in
    // focal_tip_start), for an augmented lineage from the parent it was drawn
    // from, and for an augmented lineage drawn before the first observed event
    // — its parent is one of the two crown lineages, which both still carry
    // tip_start 0 at that time.
    void compute_pendant_pd_topology(tree_t& tree)
    {
      // The alive lineages, as the multiset of their tip_start values: the
      // authority for both N (its size) and S (its sum).  A lineage is fully
      // described by its tip_start here, so the two lineages an event leaves
      // behind are interchangeable and the multiset needs no names.
      //
      // `by_id` names the lineages the augmentation can pick as a parent.  The
      // value it holds is a key into the multiset, not a second copy of the
      // state: `take` removes the alive lineage nearest to the requested
      // tip_start and returns the value it actually removed, so S stays the
      // sum of a real alive set whatever the key says.  Asking for a lineage
      // that is no longer pendant at the age claimed then costs the nearest
      // real one, instead of subtracting a tip_start no lineage has and
      // driving P above N*t.
      std::multiset<double> alive{ 0.0, 0.0 };   // the two crown lineages
      std::unordered_map<int, double> by_id;     // lineage id -> its tip_start
      double S = 0.0;

      auto take = [&alive, &S](double want) {
        if (alive.empty()) return want;
        auto it = alive.lower_bound(want);
        if (it == alive.end()) --it;
        else if (it != alive.begin()) {
          auto prev = std::prev(it);
          if ((want - *prev) < (*it - want)) it = prev;
        }
        const double got = *it;
        alive.erase(it);
        S -= got;
        return got;
      };
      auto put = [&alive, &S](double v) { alive.insert(v); S += v; };

      const size_t last = tree.size() - 1;
      for (size_t i = 0; i < tree.size(); ++i) {
        auto& node = tree[i];
        const double t = node.brts;
        node.pd = static_cast<double>(alive.size()) * t - S;
        if (detail::is_extinction(node)) {
          double want = node.tip_start;    // birth time recorded at insertion
          if (node.id >= 0) {
            auto it = by_id.find(node.id);
            if (it != by_id.end()) {
              want = it->second;           // reset by a split since that birth
              by_id.erase(it);
            }
          }
          const double ts_s = take(want);
          node.tip_start = ts_s;
          node.focal_tip_start = ts_s;
        }
        else if (i == last) {
          // The last node marks the present, not an event: no lineage is born
          // and none splits there.
          node.tip_start = t;
          node.focal_tip_start = ts_unknown;
        }
        else {
          double want = 0.0;
          bool known = true;
          if (detail::is_tip(node)) {
            // Observed event: create_tree parked the splitting lineage's
            // tip_start here.  Without a topology it is ts_unknown.
            want = node.focal_tip_start;
            if (want < 0.0) { want = 0.0; known = false; }
          }
          else if (node.parent_id >= 0) {
            auto it = by_id.find(node.parent_id);
            if (it != by_id.end()) want = it->second;
          }
          // else: drawn with no lineage on record, so its parent is a crown
          // lineage, still at tip_start 0.
          const double ts_p = take(want);      // the splitting lineage leaves
          put(t); put(t);                      // parent and daughter, both tips
          node.focal_tip_start = known ? ts_p : ts_unknown;
          node.tip_start = t;
          if (node.id >= 0) by_id[node.id] = t;
          if (node.parent_id >= 0) {
            auto it = by_id.find(node.parent_id);
            if (it != by_id.end()) it->second = t;
          }
        }
      }
    }


    // The pendant PD of a tree whose observed lineages carry no topology:
    // every observed lineage is recorded as dating from the crown, one lineage
    // per node, and P is summed over the tree at each node.  Kept so that a
    // bare branching-time vector returns the values it always has.
    void compute_pendant_pd_no_topology(tree_t& tree)
    {
      // Pass 1: set tip_start for speciation nodes, initialize focal_tip_start
      for (auto& node : tree) {
        if (!detail::is_extinction(node)) {
          node.tip_start = (node.parent_id == -1) ? 0.0 : node.brts;
        }
        // An event whose parent is not on record leaves D mean-field (E = M).
        node.focal_tip_start = (node.parent_id >= 0) ? 0.0 : ts_unknown;
      }

      // Pass 2: compute pendant PD per node
      for (auto& node : tree) {
        node.pd = detail::calculate_pendant_pd(node.brts, tree);
      }

      // Pass 3: compute focal_tip_start via lineage tracking (O(N))
      std::unordered_map<int, double> alive_ts;
      for (auto& node : tree) {
        if (detail::is_extinction(node)) {
          node.focal_tip_start = node.tip_start;
          if (node.id >= 0) alive_ts.erase(node.id);
        } else {
          if (node.parent_id >= 0) {
            auto pit = alive_ts.find(node.parent_id);
            if (pit != alive_ts.end()) {
              node.focal_tip_start = pit->second;
              pit->second = node.brts;
            }
          }
          if (node.id >= 0) alive_ts[node.id] = node.brts;
        }
      }
    }


    // After augmentation, assign tip_start, focal_tip_start, and pendant PD.
    // The last node is the observed terminal marker at the present and carries
    // the flag create_tree set: whether the observed topology was supplied.
    void compute_pendant_pd(tree_t& tree)
    {
      if (tree.empty()) return;
      if (tree.back().clade == clade_topology) compute_pendant_pd_topology(tree);
      else compute_pendant_pd_no_topology(tree);
    }


    double get_next_bt(const tree_t& tree, double cbt)
    {
      auto it = std::upper_bound(tree.cbegin(), tree.cend(), cbt, detail::node_less{});
      return (it != tree.cend()) ? it->brts : tree.back().brts;
    }


    // insert speciation node_t before t_spec,
    // inserts extinction node_t before t_ext
    // and tracks n
    void insert_species(double t_spec, double t_ext, tree_t& tree, int id = -1, int parent_id = -1)
    {
      auto n_after = [](tree_t::iterator it) {
        const auto to = detail::is_extinction(*it) ? -1.0 : 1.0;
        return it->n + to;
      };
      tree.reserve(tree.size() + 2);   // keep iterators valid
      auto first = std::lower_bound(tree.begin(), tree.end(), t_spec, detail::node_less{});
      auto n = (first != tree.begin()) ? n_after(first - 1) : tree.front().n;
      first = make_node(tree.emplace(first), t_spec, n, t_ext, id, parent_id);
      // recalculate dirty range
      for (++first; first->brts < t_ext; ++first) {
        first->n = n_after(first - 1);
      }
      make_extinct_node(tree.emplace(first), t_ext, n_after(first - 1), t_spec, id, parent_id);
    }


    // Insert an unsampled extant species: speciation node only (no extinction),
    // marked with t_ext_unsampled so sampling_prob can identify it.
    // Lineage stays alive from t_spec to T (increases n for all subsequent nodes).
    void insert_unsampled_species(double t_spec, tree_t& tree, int id = -1, int parent_id = -1)
    {
      auto n_after = [](tree_t::iterator it) {
        const auto to = detail::is_extinction(*it) ? -1.0 : 1.0;
        return it->n + to;
      };
      tree.reserve(tree.size() + 1);
      auto first = std::lower_bound(tree.begin(), tree.end(), t_spec, detail::node_less{});
      auto n = (first != tree.begin()) ? n_after(first - 1) : tree.front().n;
      first = make_node(tree.emplace(first), t_spec, n, t_ext_unsampled, id, parent_id);
      // recalculate n for all subsequent nodes (lineage stays alive to T)
      for (++first; first != tree.end(); ++first) {
        first->n = n_after(first - 1);
      }
    }


    // Number of thinning candidates whose acceptance probability exceeded 1,
    // i.e. the envelope did not dominate nh(t) at the candidate time. Read
    // and reset through thinning_envelope_violations(); augment_trees()
    // reports the count accumulated over one call.
    std::atomic<long long> envelope_violations{ 0 };

    // Factor applied to the larger endpoint rate when lambda or mu vary
    // within a segment (an M or D coefficient is active): nh(t) is then not
    // monotone in t and its maximum need not lie at an endpoint.
    //
    // The factor is not a bound. In
    // nh(t) = n * lambda(t) * (1 - rho * exp(-mu(t) * (T - t)))
    // the survival factor falls with t, so when lambda rises with t the
    // product peaks between the endpoints, and the peak exceeds twice the
    // larger endpoint value at ordinary parameters: model = "nd", linear
    // link, beta_D = -0.5 gives 8 to 15 violations per 2000 draws on the
    // 7-tip tree of tests/testthat/test-thinning-envelope.R, -2.0 gives
    // 24 to 46. A dominating envelope there is max(lambda) over the segment
    // (attained at an endpoint, since eta is affine in t within a segment)
    // times max(survival) over the segment (at most
    // 1 - rho * exp(-max(mu) * (T - cbt))), which needs lambda and mu
    // separately; Model exposes neither, nor its link and model_bin.
    // Deferred to wave 3 with H6, the larger error on the same path.
    constexpr double envelope_safety = 2.0;

    // Rounding allowance on pt <= 1 (the survival factor at the candidate
    // and at the segment start are two separate exp() evaluations).
    constexpr double pt_tolerance = 1e-12;


    // nh(t) at the start of the segment (cbt, next_bt].
    // Model::nh_rate selects the node by lower_bound(t); at t == cbt that is
    // the node at cbt itself, whose n is the count before its event. The
    // count on the segment is that of the node found by upper_bound(cbt).
    // Evaluating at the first representable time after cbt selects that
    // node and leaves the survival factor unchanged to machine precision.
    double segment_start_rate(double cbt, double next_bt, const param_t& pars, const tree_t& tree, const Model& model)
    {
      const double t = std::nextafter(cbt, next_bt);
      return std::max(0.0, model.nh_rate(t, pars, tree));
    }


    void do_augment_tree_cont(const param_t& pars, tree_t& tree, const Model& model, int max_missing, double max_lambda, int& next_id)
    {
      double cbt = 0;
      tree.reserve(5 * tree.size());    // just a guess, should cover most 'normal' cases
      int num_missing_branches = 0;
      const double b = tree.back().brts;
      // With no M or D coefficient active, lambda and mu depend on n only and
      // are constant within a segment; nh(t) = n * lambda * (1 - rho * exp(-mu * (T - t)))
      // is then decreasing in t and its value at the segment start dominates
      // it on the whole segment.
      const bool constant_rates = (pars[2] == 0.0) && (pars[3] == 0.0) &&
                                  (pars[6] == 0.0) && (pars[7] == 0.0);
      double lambda_max = 0.0;
      bool new_interval = true;   // (re)compute the envelope: start, tree changed, or next_bt reached
      while (cbt < b) {
        const double next_bt = get_next_bt(tree, cbt);
        if (new_interval) {
          const double lambda_start = segment_start_rate(cbt, next_bt, pars, tree, model);
          const double endpoint_max = constant_rates
            ? lambda_start
            : std::max(lambda_start, std::max(0.0, model.nh_rate(next_bt, pars, tree)));
          // max_lambda bounds the rate, not the envelope, so it is tested
          // against the endpoint maximum; the safety factor is applied after.
          if (endpoint_max > max_lambda) throw augmentation_lambda{};
          lambda_max = constant_rates ? endpoint_max : envelope_safety * endpoint_max;
          new_interval = false;
        }
        double next_speciation_time = next_bt;
        if (0.0 != lambda_max) {
          const double u1 = std::uniform_real_distribution<>()(reng);
          next_speciation_time = cbt - std::log(u1) / lambda_max;
        }
        if (next_speciation_time < next_bt) {
          const double u2 = std::uniform_real_distribution<>()(reng);
          double pt = std::max(0.0, model.nh_rate(next_speciation_time, pars, tree)) / lambda_max;
          if (pt > 1.0 + pt_tolerance) {
            ++envelope_violations;
          }
          // no effect on the draw: u2 is in [0, 1), so u2 < pt already held
          // for every pt >= 1. The candidate is accepted either way.
          pt = std::min(pt, 1.0);
          if (u2 < pt) {
            double ext_time = model.extinction_time(next_speciation_time, pars, tree);
            // find lineages alive at next_speciation_time and pick one as parent
            std::vector<int> alive_ids;
            for (const auto& node : tree) {
              if (!detail::is_extinction(node) &&
                  node.brts < next_speciation_time &&
                  node.t_ext > next_speciation_time) {
                alive_ids.push_back(node.id);
              }
            }
            int chosen_parent_id = -1;
            if (!alive_ids.empty()) {
              std::uniform_int_distribution<size_t> uid(0, alive_ids.size() - 1);
              chosen_parent_id = alive_ids[uid(reng)];
            }
            if (ext_time >= b) {
              // Unsampled extant species (rho < 1): insert into tree
              // so that N(t) is correct for diversity-dependent models.
              int new_id = next_id++;
              insert_unsampled_species(next_speciation_time, tree, new_id, chosen_parent_id);
              num_missing_branches++;
              if (num_missing_branches > max_missing) {
                throw augmentation_overrun{};
              }
              new_interval = true;   // tree changed
            } else {
              int new_id = next_id++;
              insert_species(next_speciation_time, ext_time, tree, new_id, chosen_parent_id);
              num_missing_branches++;
              if (num_missing_branches > max_missing) {
                throw augmentation_overrun{};
              }
              new_interval = true;   // tree changed
            }
          }
          // a rejected candidate keeps lambda_max for the rest of the segment
        }
        else {
          new_interval = true;   // next_bt reached
        }
        cbt = std::min(next_speciation_time, next_bt);
      }
      compute_pendant_pd(tree);
    }

  } // namespace augment


  long long thinning_envelope_violations(bool reset)
  {
    return reset ? envelope_violations.exchange(0) : envelope_violations.load();
  }


  void augment_tree(const param_t& pars, const tree_t& input_tree, const Model& model, int max_missing, double max_lambda, tree_t& pooled)
  {
    pooled.resize(input_tree.size());
    std::copy(input_tree.cbegin(), input_tree.cend(), pooled.begin());
    // assign sequential IDs to initial tree nodes; augmented nodes get IDs starting after
    int next_id = 0;
    for (auto& node : pooled) {
      if (!detail::is_extinction(node)) {
        node.id = next_id++;
        node.parent_id = -1;
      } else {
        node.id = -1;
        node.parent_id = -1;
      }
    }
    do_augment_tree_cont(pars, pooled, model, max_missing, max_lambda, next_id);
  }


  // returns one augmented tree per vpars
  // failures results in empty tree
  std::vector<tree_t> augment_trees(const std::vector<param_t>& vpars, const tree_t& input_tree, const Model& model, int max_missing, double max_lambda, int num_threads)
  {
    if (!model.is_threadsafe()) num_threads = 1;
    num_threads = std::max(1, std::min(num_threads, static_cast<int>(std::thread::hardware_concurrency())));
    tbb::task_arena arena(std::max(1, num_threads));
    std::vector<tree_t> trees(vpars.size(), input_tree);
    // Integer division gives 0 when there are fewer parameter sets than
    // threads, which is not a valid grainsize.
    const size_t grainsize = std::max<size_t>(1, vpars.size() / static_cast<size_t>(std::max(1, num_threads)));
    // Started outside the arena, a parallel algorithm runs on the default one
    // and ignores num_threads.
    arena.execute([&] {
    tbb::parallel_for(tbb::blocked_range<size_t>(0ull, vpars.size(), grainsize), [&](const tbb::blocked_range<size_t>& r) {
      for (size_t i = r.begin(); i < r.end(); ++i) {
        try {
          int next_id = 0;
          for (auto& node : trees[i]) {
            if (!detail::is_extinction(node)) {
              node.id = next_id++;
              node.parent_id = -1;
            } else {
              node.id = -1;
              node.parent_id = -1;
            }
          }
          do_augment_tree_cont(vpars[i], trees[i], model, max_missing, max_lambda, next_id);
        } catch (...) {
          trees[i].clear();
        }
      }
    });
    });
    return trees;
  }


}

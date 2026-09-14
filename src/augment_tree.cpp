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
    // segment that ends there, which is the count node.n carries.
    //
    // The sweep consumes one event at a time and reports P as it goes, so the
    // augmentation loop can carry it forward across the lineages it inserts
    // rather than re-deriving P from the whole event list at every candidate
    // time.  That is what keeps node.pd — the only thing Model reads P off —
    // current while a tree is still being built.
    //
    // Two conventions, chosen once per tree by create_tree and recorded in the
    // `clade` flag of the last node:
    //
    //   topology  the lineages the observed tree really has: two crown lineages
    //             at tip_start 0 and a reset at every split, of the splitting
    //             lineage as well as of its daughter.  The splitting lineage is
    //             on record at every event: an observed event names it by the
    //             tip_start the observed topology gives it (create_tree parked
    //             that in focal_tip_start), an augmented one by the id of the
    //             lineage it was drawn from — a node id, or one of the two
    //             crown ids.
    //
    //   legacy    what a bare branching-time vector can support: one lineage per
    //             node, every observed lineage dating from the crown, the two
    //             crown lineages absent from the sum, and the node counting
    //             itself.  Kept so that an input with no topology returns the
    //             values it always has; the sweep is the O(N) form of that same
    //             sum.
    //
    // Both conventions carry the same alive set of named lineages, because that
    // is what the augmentation draws a parent from; only what the sum P is taken
    // over differs.
    class pendant_sweep
    {
    public:
      // A lineage alive at the sweep's current time.
      //
      //   ts      its tip_start now: the last time it became a pendant tip,
      //           reset by every split it makes, observed or augmented.
      //   obs_ts  its tip_start as the observed topology records it, which an
      //           augmented birth does not touch.  parent_tip_start is measured
      //           on the observed tree alone, so it is obs_ts, not ts, that
      //           identifies the lineage splitting at an observed event.
      //   observed  true for a lineage of the observed tree (the two crown
      //           lineages and every lineage born at an observed branching),
      //           false for an augmented one.  The proposal's parent draw and
      //           the -log(2*tips + Ne) the density charges count an observed
      //           lineage twice and an augmented one once.
      struct lineage { double ts; double obs_ts; int id; bool observed; };

      explicit pendant_sweep(bool topology) : topology_(topology)
      {
        alive_.push_back(lineage{ 0.0, 0.0, crown_id_a, true });
        alive_.push_back(lineage{ 0.0, 0.0, crown_id_b, true });
      }

      // The lineages alive now, in the order they were born.  The augmentation
      // draws its parent from this.
      const std::vector<lineage>& alive() const { return alive_; }

      // The number of labelled attachments a birth now has: 2 per observed
      // lineage alive, 1 per augmented one.  This is the 2*tips + Ne of
      // Model::sampling_prob, counted from the alive set rather than from the
      // node list, and the sampler and the density must agree on it.
      long long attachments() const
      {
        long long k = 0;
        for (const auto& l : alive_) k += l.observed ? 2 : 1;
        return k;
      }

      // The pendant PD this sweep will store on `node`.
      double pd_of(const node_t& node) const
      {
        return topology_ ? P(node.brts)                      // before its event
                         : P(node.brts) + legacy_delta(node); // after it
      }

      // Consume the node's event: assign its tip_start and focal_tip_start, and
      // move the alive set past it.  Each event is applied exactly once.
      void advance(node_t& node, bool last)
      {
        if (topology_) advance_topology(node, last);
        else           advance_legacy(node);
      }

    private:
      using alive_it = std::vector<lineage>::iterator;

      // P over the lineages the sweep currently holds alive.
      double P(double t) const { return count() * t - S_; }
      double count() const
      {
        return topology_ ? static_cast<double>(alive_.size()) : count_;
      }

      // The alive lineage with this id.  Every lineage the augmentation can
      // name has one: a node id for a lineage born at a node, a crown id for
      // one of the two the crown starts with.
      alive_it by_id(int id)
      {
        return std::find_if(alive_.begin(), alive_.end(),
                            [id](const lineage& l) { return l.id == id; });
      }

      // The observed lineage the topology means when it reports that the
      // lineage splitting here was last a tip at `want`.
      //
      // The match is on obs_ts — what the observed tree records — because that
      // is what parent_tip_start is measured from, and it is nearest rather
      // than exact because the caller's tip starts and the caller's branching
      // times are two separate floating-point computations over the same tree.
      // Two observed lineages share an obs_ts exactly when they are the two
      // halves of one observed split, which the topology does not tell apart;
      // the first-born of them is taken, in both sweeps alike, so the identity
      // the sweep settles on is the same one every time it replays the tree.
      alive_it nearest_observed(double want)
      {
        auto best = alive_.end();
        double best_d = 0.0;
        for (auto it = alive_.begin(); it != alive_.end(); ++it) {
          if (!it->observed) continue;
          const double d = std::abs(it->obs_ts - want);
          if (best == alive_.end() || d < best_d) { best = it; best_d = d; }
        }
        return best;
      }

      void advance_topology(node_t& node, bool last)
      {
        const double t = node.brts;
        if (detail::is_extinction(node)) {
          auto it = by_id(node.id);
          const double ts_s = (it != alive_.end()) ? it->ts : node.tip_start;
          if (it != alive_.end()) { S_ -= it->ts; alive_.erase(it); }
          node.tip_start = ts_s;
          node.focal_tip_start = ts_s;
        }
        else if (last) {
          // The last node marks the present, not an event: no lineage is born
          // and none splits there.
          node.tip_start = t;
          node.focal_tip_start = ts_unknown;
        }
        else {
          const bool observed = detail::is_tip(node);
          auto it = alive_.end();
          bool known = true;
          if (observed) {
            double want = node.focal_tip_start;
            if (want < 0.0) { want = 0.0; known = false; }
            it = nearest_observed(want);
          }
          else {
            it = by_id(node.parent_id);
            if (it == alive_.end()) known = false;   // no lineage on record
          }
          double ts_p = 0.0;
          if (it != alive_.end()) {
            ts_p = it->ts;
            S_ += t - ts_p;                // the splitting lineage is a tip again
            it->ts = t;
            // An observed split is the only kind the observed tree records.
            if (observed) it->obs_ts = t;
          }
          alive_.push_back(lineage{ t, t, node.id, observed });   // the daughter
          S_ += t;
          node.focal_tip_start = known ? ts_p : ts_unknown;
          node.tip_start = t;
        }
      }

      // The tip_start the legacy convention gives a node, and the change its
      // own arrival or departure makes to P at its own time.
      //
      // An observed node has no parent on record — every observed lineage
      // dates from the crown under this convention — and an augmented one
      // names the lineage it was drawn from, so it dates from its own birth.
      static double legacy_ts(const node_t& node)
      {
        return has_parent(node.parent_id) ? node.brts : 0.0;
      }

      double legacy_delta(const node_t& node) const
      {
        if (detail::is_extinction(node)) {
          auto it = born_.find(node.id);
          const double ts = (node.id >= 0 && it != born_.end()) ? it->second
                                                                : node.tip_start;
          return -(node.brts - ts);         // the dying lineage no longer counts
        }
        return node.brts - legacy_ts(node); // the node counts itself
      }

      void advance_legacy(node_t& node)
      {
        const double t = node.brts;
        if (detail::is_extinction(node)) {
          double ts = node.tip_start;
          auto it = born_.find(node.id);
          if (node.id >= 0 && it != born_.end()) { ts = it->second; born_.erase(it); }
          count_ -= 1.0;
          S_ -= ts;
          node.focal_tip_start = node.tip_start;
          if (node.id >= 0) reset_.erase(node.id);
          auto a = by_id(node.id);
          if (a != alive_.end()) alive_.erase(a);
        }
        else {
          node.tip_start = legacy_ts(node);
          // An event whose parent is not on record leaves D mean-field (E = M).
          node.focal_tip_start = has_parent(node.parent_id) ? 0.0 : ts_unknown;
          if (has_parent(node.parent_id)) {
            auto pit = reset_.find(node.parent_id);
            if (pit != reset_.end()) {
              node.focal_tip_start = pit->second;
              pit->second = t;
            }
          }
          count_ += 1.0;
          S_ += node.tip_start;
          if (node.id >= 0) { born_[node.id] = node.tip_start; reset_[node.id] = t; }
          alive_.push_back(lineage{ t, t, node.id, detail::is_tip(node) });
        }
      }

      const bool topology_;
      // The alive lineages, in birth order, crown lineages first.  Under the
      // topology convention it is also the authority for N (its size) and S
      // (the sum of its tip starts); under the legacy convention P is the sum
      // over nodes that count_ and S_ carry, and alive_ only names the
      // lineages a birth can be drawn from.
      std::vector<lineage> alive_;
      std::unordered_map<int, double> born_;      // legacy: id -> the ts it counts with
      // legacy: id -> the tip_start a split of that lineage reports.  The two
      // crown lineages are seeded at 0, the tip_start the legacy convention
      // gives every observed lineage, so that a birth drawn from one of them
      // reads the same E as a birth drawn from any other observed lineage.
      std::unordered_map<int, double> reset_{ { crown_id_a, 0.0 },
                                              { crown_id_b, 0.0 } };
      double count_ = 0.0;                        // legacy: lineages alive
      double S_ = 0.0;
    };


    // Assign tip_start, focal_tip_start and pendant PD over a whole tree.  The
    // last node is the observed terminal marker at the present and carries the
    // flag create_tree set: whether the observed topology was supplied.
    void compute_pendant_pd(tree_t& tree)
    {
      if (tree.empty()) return;
      pendant_sweep sweep(tree.back().clade == clade_topology);
      const size_t last = tree.size() - 1;
      for (size_t i = 0; i < tree.size(); ++i) {
        tree[i].pd = sweep.pd_of(tree[i]);
        sweep.advance(tree[i], i == last);
      }
    }


    // The candidate set and the parent at each augmented birth, by replaying
    // the finished tree through the sweep the sampler carried forward.
    std::vector<attachment_t> report(const tree_t& tree)
    {
      std::vector<attachment_t> out;
      if (tree.empty()) return out;
      tree_t work(tree.cbegin(), tree.cend());   // advance() writes the tip starts
      pendant_sweep sweep(work.back().clade == clade_topology);
      const size_t last = work.size() - 1;
      for (size_t i = 0; i < work.size(); ++i) {
        node_t& node = work[i];
        if (i != last && !detail::is_extinction(node) && !detail::is_tip(node)) {
          const auto& alive = sweep.alive();
          bool found = false;
          for (const auto& l : alive) if (l.id == node.parent_id) { found = true; break; }
          out.push_back(attachment_t{ node.brts, node.n,
                                      static_cast<long long>(alive.size()),
                                      sweep.attachments(), node.parent_id, found });
        }
        sweep.advance(node, i == last);
      }
      return out;
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
    //
    // The envelope is nh at the start of the segment, which dominates nh on
    // the whole segment: Model::nh_rate holds lambda and mu at the segment's
    // own values, so the only thing that carries t is the survival factor
    // 1 - rho * exp(-mu * (T - t)), and it falls with t. The counter is kept
    // as a running assertion on that, not as a record of a known gap.
    std::atomic<long long> envelope_violations{ 0 };

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
      double lambda_max = 0.0;
      bool new_interval = true;   // (re)compute the envelope: start, tree changed, or next_bt reached
      // The pendant-age state of the segment the sampler is in.  Model reads P
      // off node.pd, which create_tree leaves at zero and an insertion cannot
      // know, so the sweep is carried forward with cbt and the node that governs
      // the current segment is stamped with its pd before any rate on that
      // segment is evaluated.  Every event is consumed exactly once, in forward
      // time, so the values the sampler saw are the ones compute_pendant_pd
      // writes at the end.
      pendant_sweep sweep(tree.back().clade == clade_topology);
      while (cbt < b) {
        auto next_it = std::upper_bound(tree.begin(), tree.end(), cbt, detail::node_less{});
        if (next_it == tree.end()) next_it = tree.end() - 1;
        const double next_bt = next_it->brts;
        if (new_interval) {
          // P on this segment, for Model::pendant_pd to extrapolate from.
          next_it->pd = sweep.pd_of(*next_it);
          // nh falls within a segment, so the rate at its start is the
          // dominating envelope — there is nothing to inflate, and the
          // envelope and the bound max_lambda test the same number.
          const double lambda_start = segment_start_rate(cbt, next_bt, pars, tree, model);
          if (lambda_start > max_lambda) throw augmentation_lambda{};
          lambda_max = lambda_start;
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
            // The parent: uniform over the labelled attachments, which is the
            // K = 2*tips + Ne that Model::sampling_prob charges -log K for.
            //
            // An attachment is a complete labelled history, and a birth off an
            // observed lineage has two of them — the observed branch and the
            // missing one are the two halves of one split, and either half may
            // carry the observed lineage's own label — while a birth off an
            // augmented lineage has one, because there the two halves differ by
            // their death times and the swap is a different augmentation, not a
            // second name for this one.  So an observed lineage alive weighs 2
            // and an augmented one weighs 1, and the sampler and the density
            // count the same K.
            //
            // The sweep is the authority for which lineages are alive: it holds
            // every event up to cbt and none after it, and the two crown
            // lineages, which carry no node, are in it from the start.
            const auto& alive = sweep.alive();
            const long long attachments = sweep.attachments();
            int chosen_parent_id = no_parent;
            {
              std::uniform_int_distribution<long long> uid(0, attachments - 1);
              long long k = uid(reng);
              for (const auto& l : alive) {
                const long long w = l.observed ? 2 : 1;
                if (k < w) { chosen_parent_id = l.id; break; }
                k -= w;
              }
            }
            int new_id = next_id++;
            if (ext_time >= b) {
              // Unsampled extant species (rho < 1): insert into tree
              // so that N(t) is correct for diversity-dependent models.
              insert_unsampled_species(next_speciation_time, tree, new_id, chosen_parent_id);
            } else {
              insert_species(next_speciation_time, ext_time, tree, new_id, chosen_parent_id);
            }
            // The birth is an event at cbt: consume it, so the state stays that
            // of the segment the sampler moves into.  The extinction node the
            // insertion parked in the future is consumed when cbt reaches it.
            auto born = std::lower_bound(tree.begin(), tree.end(), next_speciation_time,
                                         detail::node_less{});
            sweep.advance(*born, false);
            num_missing_branches++;
            if (num_missing_branches > max_missing) {
              throw augmentation_overrun{};
            }
            new_interval = true;   // tree changed
          }
          // a rejected candidate keeps lambda_max for the rest of the segment
        }
        else {
          // next_bt reached: consume its event before moving onto the next segment.
          sweep.advance(*next_it, next_it == tree.end() - 1);
          new_interval = true;
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


  std::vector<attachment_t> attachment_report(const tree_t& tree)
  {
    return report(tree);
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

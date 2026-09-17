// Evolutionary distinctiveness (ED) of every lineage alive at a time t, on the
// complete tree: the fair-proportion score
//
//   ED(s, t) = sum over the branches b on the path from the crown to s of
//              length(b) / (number of alive lineages descending from b)
//
// (Isaac et al. 2007).  Two identities tie it to the covariates the model
// already has.  The last branch on the path -- from s's most recent split to
// t -- has one alive descendant, s itself, so it contributes the pendant age
// E(s, t) = t - tip_start(s); ED truncated to that term is the EP model.  And
// the sum of ED over the alive lineages is the total length of the branches
// with an alive descendant, Faith's phylogenetic diversity P_tree(t).
//
// The tree is given as a lineage forest, the representation both the
// likelihood (a node list with lineage ids) and the forward simulator (an
// L-table) reduce to: for lineage i, parent[i] is the index of the lineage it
// split from (no_parent_idx for a crown lineage), birth[i] the time of that
// split, and alive[i] whether it is alive on the segment the caller is
// evaluating.  A lineage dead by t still carries the branches its living
// descendants inherit, so it is kept; a lineage born after t does not exist
// yet and is skipped.
//
// Within a segment between events no lineage is born and none dies, so every
// ancestral share is constant and only the pendant branch grows, at rate one.
// The routine therefore returns, per alive lineage, the constant c[i] and the
// tip start ts[i] with
//
//   ED_i(u) = c[i] + (u - ts[i])        for u in the segment,
//
// which is what the compensator integrates in closed form and what the event
// term reads at u = the event time.  ts[i] is the time of i's most recent
// split at or before t (its birth if it has not split), the same tip_start
// the pendant sweep assigns.
//
// Cost.  The likelihood evaluates every segment of a tree at every parameter
// value the M-step tries, so the per-segment pass must not allocate or sort.
// The structure that depends on the tree alone -- lineages in birth order and
// each lineage's children in birth order -- is built once per tree
// (forest_index); a pass then reads a prefix of that order (the lineages born
// by t: the order is sorted by birth) with two O(N) sweeps over reused
// workspace.  fair_proportion() without an index is the convenience form for
// a single evaluation (the simulator, the R export).
//
// Shared by inst/include/model.hpp (the likelihood, via the node list) and
// inst/include/general_tree.hpp (the simulator, via its L-table), so the two
// cannot disagree about what ED is.

#ifndef EMPHASIS_ED_COVARIATE_HPP_INCLUDED
#define EMPHASIS_ED_COVARIATE_HPP_INCLUDED

#include <vector>
#include <algorithm>
#include <limits>
#include <cstddef>

namespace emphasis {
  namespace ed {

    constexpr int no_parent_idx = -1;

    // The tree-only structure: birth order and children lists (CSR layout,
    // children of lineage i are child_list[child_begin[i] .. child_begin[i+1])
    // in birth order).
    struct forest_index {
      std::vector<int>    order;        // lineage indices sorted by birth (stable)
      std::vector<double> order_birth;  // birth[order[k]]
      std::vector<int>    child_begin;  // size n + 1
      std::vector<int>    child_list;

      void build(const std::vector<int>& parent, const std::vector<double>& birth) {
        const std::size_t n = parent.size();
        order.resize(n);
        for (std::size_t i = 0; i < n; ++i) order[i] = static_cast<int>(i);
        std::stable_sort(order.begin(), order.end(),
                         [&](int a, int b) { return birth[static_cast<std::size_t>(a)] < birth[static_cast<std::size_t>(b)]; });
        order_birth.resize(n);
        for (std::size_t k = 0; k < n; ++k) order_birth[k] = birth[static_cast<std::size_t>(order[k])];
        child_begin.assign(n + 1, 0);
        for (std::size_t i = 0; i < n; ++i)
          if (parent[i] >= 0) ++child_begin[static_cast<std::size_t>(parent[i]) + 1];
        for (std::size_t i = 0; i < n; ++i) child_begin[i + 1] += child_begin[i];
        child_list.assign(n, -1);
        std::vector<int> fill(child_begin.begin(), child_begin.end() - 1);
        // walking the birth order fills each lineage's children in birth order
        for (int i : order)
          if (parent[static_cast<std::size_t>(i)] >= 0)
            child_list[static_cast<std::size_t>(fill[static_cast<std::size_t>(parent[static_cast<std::size_t>(i)])]++)] = i;
      }

      // Number of lineages born at or before t: a prefix of `order`.
      std::size_t prefix(double t) const {
        return static_cast<std::size_t>(std::upper_bound(order_birth.begin(), order_birth.end(), t) - order_birth.begin());
      }
    };

    // Per-call outputs and reusable scratch.  c[i], ts[i] are NaN for a
    // lineage that is not alive (or not yet born); n_desc[i] is the number of
    // alive lineages in i's subtree, i included.
    struct result_t {
      std::vector<double> c;
      std::vector<double> ts;
      std::vector<int>    n_desc;
      std::vector<double> inherit;   // scratch
    };

    // The pass, over the lineages born by t, with the index built once.
    inline void fair_proportion(const forest_index& ix,
                                const std::vector<int>&    parent,
                                const std::vector<double>& birth,
                                const std::vector<char>&   alive,
                                double t,
                                result_t& out)
    {
      (void)parent;
      const std::size_t n = birth.size();
      const double nan = std::numeric_limits<double>::quiet_NaN();
      out.c.assign(n, nan);
      out.ts.assign(n, nan);
      out.n_desc.assign(n, 0);
      out.inherit.assign(n, 0.0);
      const std::size_t m = ix.prefix(t);

      // Alive descendants, bottom-up: children before parents.  A child born
      // after t is outside the prefix and contributes nothing.
      for (std::size_t k = m; k-- > 0;) {
        const int i = ix.order[k];
        int s = alive[static_cast<std::size_t>(i)] ? 1 : 0;
        for (int p = ix.child_begin[static_cast<std::size_t>(i)]; p < ix.child_begin[static_cast<std::size_t>(i) + 1]; ++p) {
          const int c = ix.child_list[static_cast<std::size_t>(p)];
          if (birth[static_cast<std::size_t>(c)] <= t) s += out.n_desc[static_cast<std::size_t>(c)];
        }
        out.n_desc[static_cast<std::size_t>(i)] = s;
      }

      // Top-down: what each lineage inherits at its birth, then its own
      // branches up to each child's birth and up to t.  The children are in
      // birth order, so the suffix sums of their subtrees are formed by one
      // pass from the last child backwards, stored in place of a scratch
      // array by walking twice over a short list.
      for (std::size_t k = 0; k < m; ++k) {
        const std::size_t i = static_cast<std::size_t>(ix.order[k]);
        const int b0 = ix.child_begin[i], b1 = ix.child_begin[i + 1];
        // children born by t are a prefix of [b0, b1)
        int e = b0;
        while (e < b1 && birth[static_cast<std::size_t>(ix.child_list[static_cast<std::size_t>(e)])] <= t) ++e;
        const int self = alive[i] ? 1 : 0;
        // total alive below all children born by t, then peel from the front
        int suf = 0;
        for (int p = b0; p < e; ++p) suf += out.n_desc[static_cast<std::size_t>(ix.child_list[static_cast<std::size_t>(p)])];
        double pre = 0.0;
        double from = birth[i];
        const double inh = out.inherit[i];
        for (int p = b0; p < e; ++p) {
          const std::size_t c = static_cast<std::size_t>(ix.child_list[static_cast<std::size_t>(p)]);
          const double to = birth[c];
          const int desc = self + suf;                 // alive lineages below this branch
          if (desc > 0) pre += (to - from) / static_cast<double>(desc);
          from = to;
          out.inherit[c] = inh + pre;                  // the child inherits every branch above its birth
          suf -= out.n_desc[c];
        }
        if (alive[i]) {
          out.c[i]  = inh + pre;
          out.ts[i] = from;                            // the most recent split, or the birth
        }
      }
    }

    // Convenience form: one evaluation, index built here.
    inline void fair_proportion(const std::vector<int>&    parent,
                                const std::vector<double>& birth,
                                const std::vector<char>&   alive,
                                double t,
                                result_t& out)
    {
      forest_index ix;
      ix.build(parent, birth);
      fair_proportion(ix, parent, birth, alive, t, out);
    }

  }
}

#endif

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
// split from (no_parent for a crown lineage), birth[i] the time of that split,
// and alive[i] whether it is alive on the segment the caller is evaluating.
// A lineage dead by t still carries the branches its living descendants
// inherit, so it is kept; a lineage born after t does not exist yet and is
// skipped.
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
// One top-down pass: for each lineage in birth order, the branches between
// its consecutive splits, each divided by the alive lineages below it (the
// lineage itself if alive, plus the subtrees of the children born at or after
// the branch's end), accumulate into a prefix that every child inherits at
// its own birth.  O(n log n) for the child sort, O(n) otherwise.
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

    // c[i], ts[i] for every lineage; c is NaN for a lineage that is not alive
    // (or not yet born).  n_desc[i] is the number of alive lineages in i's
    // subtree, i included, which the simulator reports and tests check.
    struct result_t {
      std::vector<double> c;
      std::vector<double> ts;
      std::vector<int>    n_desc;
    };

    inline void fair_proportion(const std::vector<int>&    parent,
                                const std::vector<double>& birth,
                                const std::vector<char>&   alive,
                                double t,
                                result_t& out)
    {
      const std::size_t n = parent.size();
      const double nan = std::numeric_limits<double>::quiet_NaN();
      out.c.assign(n, nan);
      out.ts.assign(n, nan);
      out.n_desc.assign(n, 0);

      // Lineages in birth order; a parent is born before its children, so
      // this order is also topological (crowns first).
      std::vector<std::size_t> order;
      order.reserve(n);
      for (std::size_t i = 0; i < n; ++i) if (birth[i] <= t) order.push_back(i);
      std::stable_sort(order.begin(), order.end(),
                       [&](std::size_t a, std::size_t b) { return birth[a] < birth[b]; });

      // children[i]: indices of i's children born at or before t, in birth
      // order (order is sorted by birth, so appending in that order sorts them).
      std::vector<std::vector<std::size_t>> children(n);
      for (std::size_t i : order) {
        if (parent[i] >= 0) children[static_cast<std::size_t>(parent[i])].push_back(i);
      }

      // Alive descendants, bottom-up: children before parents.
      for (auto it = order.rbegin(); it != order.rend(); ++it) {
        const std::size_t i = *it;
        int s = alive[i] ? 1 : 0;
        for (std::size_t c : children[i]) s += out.n_desc[c];
        out.n_desc[i] = s;
      }

      // Top-down: what each lineage inherits from its ancestry at its birth,
      // then its own branches up to each child's birth and up to t.
      std::vector<double> inherit(n, 0.0);
      for (std::size_t i : order) {
        const auto& ch = children[i];
        const std::size_t m = ch.size();
        // suffix sums of the children's alive subtrees: suf[k] = sum_{j >= k} n_desc[ch[j]]
        std::vector<int> suf(m + 1, 0);
        for (std::size_t k = m; k-- > 0;) suf[k] = suf[k + 1] + out.n_desc[ch[k]];
        const int self = alive[i] ? 1 : 0;
        double pre = 0.0;          // sum over the closed branches so far
        double from = birth[i];    // start of the current branch
        for (std::size_t k = 0; k < m; ++k) {
          const double to = birth[ch[k]];
          const int desc = self + suf[k];        // alive lineages below this branch
          if (desc > 0) pre += (to - from) / static_cast<double>(desc);
          from = to;
          inherit[ch[k]] = inherit[i] + pre;      // the child inherits every branch above its birth
        }
        if (alive[i]) {
          out.c[i]  = inherit[i] + pre;
          out.ts[i] = from;                       // the most recent split, or the birth
        }
      }
    }

  }
}

#endif

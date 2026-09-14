#ifndef EMPHASIS_UNPACK_H_INCLUDED
#define EMPHASIS_UNPACK_H_INCLUDED

// [[Rcpp::plugins(cpp14)]]

#include <Rcpp.h>
#include "emphasis.hpp"


using namespace Rcpp;


DataFrame unpack(const emphasis::tree_t& tree)
{
NumericVector brts, n, t_ext, pd, tip_start, focal_tip_start;
IntegerVector clade, id, parent_id;
for (const emphasis::node_t& node : tree) {
    brts.push_back(node.brts);
    n.push_back(node.n);
    t_ext.push_back(node.t_ext);
    pd.push_back(node.pd);
    tip_start.push_back(node.tip_start);
    focal_tip_start.push_back(node.focal_tip_start);
    // Which pendant-age convention the tree was built under.  It is carried so
    // that a replay of the forward sweep over the returned data frame is the
    // sweep the sampler itself ran, rather than always the legacy one.
    clade.push_back(node.clade);
    id.push_back(node.id);
    parent_id.push_back(node.parent_id);
}
return DataFrame::create(Named("brts") = brts,
                            Named("n") = n,
                            Named("t_ext") = t_ext,
                            Named("pd") = pd,
                            Named("tip_start") = tip_start,
                            Named("focal_tip_start") = focal_tip_start,
                            Named("clade") = clade,
                            Named("id") = id,
                            Named("parent_id") = parent_id);
}


#endif

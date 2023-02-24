#ifndef EMPHASIS_UNPACK_H_INCLUDED
#define EMPHASIS_UNPACK_H_INCLUDED

// [[Rcpp::plugins(cpp14)]]

#include <Rcpp.h>
#include "emphasis.hpp"


using namespace Rcpp;


DataFrame unpack(const emphasis::tree_t& tree)
{
NumericVector brts, n, t_ext, pd;
for (const emphasis::node_t& node : tree) {
    brts.push_back(node.brts);
    n.push_back(node.n);
    t_ext.push_back(node.t_ext);
    pd.push_back(node.pd);
}
return DataFrame::create(Named("brts") = brts, 
                            Named("n") = n, 
                            Named("t_ext") = t_ext,
                            Named("pd") = pd);
}


#endif

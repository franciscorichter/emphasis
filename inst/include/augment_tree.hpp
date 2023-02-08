#ifndef EMP_AUGMENT_TREE_HPP_INCLUDED
#define EMP_AUGMENT_TREE_HPP_INCLUDED

#include <vector>
#include "emphasis.hpp"

namespace emphasis {

  // thrown if missing branches exceeded
  class augmentation_overrun : public std::runtime_error
  {
  public:
    augmentation_overrun() : std::runtime_error("number of missing branches exceeded") {}
  };


  // thrown if lambda exceeded
  class augmentation_lambda : public std::runtime_error
  {
  public:
    augmentation_lambda() : std::runtime_error("lambda exceeded") {}
  };


  void augment_tree(const param_t& pars,
                    const tree_t& input_tree,
                    const Model& model,
                    int max_missing,
                    double max_lambda,
                    tree_t& out);

}

#endif

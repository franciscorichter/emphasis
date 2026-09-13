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


  // Number of thinning candidates whose acceptance probability exceeded 1
  // (the envelope did not dominate nh at the candidate) since the last reset.
  long long thinning_envelope_violations(bool reset);


  void augment_tree(const param_t& pars,
                    const tree_t& input_tree,
                    const Model& model,
                    int max_missing,
                    double max_lambda,
                    tree_t& out);


  // returns augmented tree per vpars
  // failures results in empty tree
  std::vector<tree_t> augment_trees(const std::vector<param_t>& vpars, 
                                    const tree_t& input_tree, 
                                    const Model& model, 
                                    int max_missing, 
                                    double max_lambda,
                                    int num_threads);

}

#endif

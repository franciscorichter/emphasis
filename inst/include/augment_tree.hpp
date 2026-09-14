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


  // What the proposal had to choose from at one augmented birth: the lineages
  // alive there, the labelled attachments they carry, and the parent it
  // recorded.  Reported by replaying the finished tree through the same sweep
  // the sampler carried forward, so a test can hold the sampler's candidate set
  // against the alive count node.n and against the 2*tips + Ne the density
  // charges.
  struct attachment_t
  {
    double brts;          // the birth
    double n;             // node.n there: the lineages alive
    long long candidates; // lineages the proposal could draw from
    long long attachments;// labelled attachments: 2 per observed, 1 per augmented
    int parent_id;        // the parent it recorded
    bool parent_alive;    // was that parent one of the candidates
  };

  std::vector<attachment_t> attachment_report(const tree_t& tree);


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

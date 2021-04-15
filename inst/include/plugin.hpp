/* Copyright (c) 2007-2014 Massachusetts Institute of Technology
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

//
// C++-API for plug-in diversification model
// Hanno 2020
//

#ifndef EMPHASIS_PLUGIN_HPP_INCLUDED
#define EMPHASIS_PLUGIN_HPP_INCLUDED

#include <memory>
#include <vector>
#include <functional>
#include "plugin.h"   // emp_node_t


namespace emphasis {


  using param_t = std::vector<double>;                  // unspecific parameters
  using node_t = emp_node_t;                            // brts, n, t_ext, pd
  using tree_t = std::vector<node_t>;                   // tree, sorted by note_t::brts
  constexpr double t_ext_tip = emp_t_ext_tip;           // t_ext for present nodes
  constexpr double t_ext_extinct = emp_t_ext_extinct;   // t_ext for extinction nodes

  
  

  // abstract diversification model
  class Model
  {
  public:
    Model() = default;
    virtual ~Model() = default;

    // optional 
    virtual const char* description() const { return "not set"; }   // textual description of the model
    virtual bool is_threadsafe() const { return false; }            // is this implementation thread-save?
    virtual bool numerical_max_lambda() const { return true; }
    virtual int nparams() const = 0;                                // number of parameters

    // diversification model
    virtual double extinction_time(double t_speciations, const param_t& pars, const tree_t& tree) const = 0;
    virtual double nh_rate(double t, const param_t& pars, const tree_t& tree) const = 0;
    virtual double sampling_prob(const param_t& pars, const tree_t& tree) const = 0;
    virtual double loglik(const param_t& pars, const tree_t& tree) const = 0;

    // optional hints for optimization step
    virtual param_t lower_bound() const { return param_t(); }
    virtual param_t upper_bound() const { return param_t(); }
  };

}

#endif

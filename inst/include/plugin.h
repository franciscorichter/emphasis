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

/***********************************************************/
/* C-API for plug-in diversification models                */
/* Hanno 2020                                              */
/***********************************************************/


#ifndef EMPHASIS_PLUGIN_H_INCLUDED
#define EMPHASIS_PLUGIN_H_INCLUDED

#include <stddef.h>


#if (defined(_WIN32) || defined(__WIN32__)) && !defined(__LCC__)
# ifdef EMP_BUILD_STANDALONE_CPP
#   define EMP_EXTERN(T) extern "C" __declspec(dllexport) T
# else
#   define EMP_EXTERN(T) extern "C" T
# endif
#else
# define EMP_EXTERN(T) extern "C" T
#endif


/* tree node */
struct emp_node_t
{
  double brts;
  double n;         /* n[i] = number of species in [time_i-1, time_i) */
  double t_ext;     /* emp_t_ext_tip for present-day species;  emp_t_ext_extinct for extinction nodes */
  double pd;
  int clade; // required for sim_tree
};


#define emp_t_ext_tip 10e10     /* t_ext for present nodes */
#define emp_t_ext_extinct 0.0   /* t_ext for extinction nodes */


typedef const char* (*emp_description_func)();
typedef bool (*emp_is_threadsafe_func)();
typedef bool (*emp_numerical_max_lambda_func)();
typedef int (*emp_nparams_func)();


/* diversification model */
typedef double (*emp_extinction_time_func)(double, const double*, unsigned, const emp_node_t*);
typedef double (*emp_nh_rate_func)(double, const double*, unsigned, const emp_node_t*);
typedef double (*emp_sampling_prob_func)(const double*, unsigned, const emp_node_t*);
typedef double (*emp_loglik_func)(const double*, unsigned, const emp_node_t*);


/* optional hints for optimizer */
typedef void (*emp_lower_bound_func)(double*);
typedef void (*emp_upper_bound_func)(double*);


#endif

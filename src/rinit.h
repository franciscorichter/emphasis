#ifndef EMPHASIS_RINIT_H_INCLUDED
#define EMPHASIS_RINIT_H_INCLUDED

#include <nlopt.h>


#if (defined(_WIN32) || defined(__WIN32__)) && !defined(__LCC__)
#  define REMP_EXPORT extern "C" __declspec(dllexport)
#else
#  define REMP_EXPORT extern "C"
#endif


#ifdef __cplusplus
extern "C" {
#endif


/* some function pointers into nlopt to be filled in by R */


REMP_EXPORT nlopt_opt(*remp_create)(nlopt_algorithm, unsigned);
REMP_EXPORT void(*remp_destroy)(nlopt_opt);
REMP_EXPORT nlopt_result(*remp_optimize)(nlopt_opt, double *, double *);
REMP_EXPORT nlopt_result(*remp_set_min_objective)(nlopt_opt, nlopt_func, void *);
REMP_EXPORT nlopt_result(*remp_set_max_objective)(nlopt_opt, nlopt_func, void *);
REMP_EXPORT nlopt_result(*remp_set_lower_bounds)(nlopt_opt, const double *);
REMP_EXPORT nlopt_result(*remp_set_lower_bounds1)(nlopt_opt, double);
REMP_EXPORT nlopt_result(*remp_set_upper_bounds)(nlopt_opt, const double *);
REMP_EXPORT nlopt_result(*remp_set_upper_bounds1)(nlopt_opt, double);
REMP_EXPORT nlopt_result(*remp_set_xtol_rel)(nlopt_opt, double);
REMP_EXPORT nlopt_result(*remp_set_xtol_abs)(nlopt_opt, double);


#ifdef __cplusplus
}
#endif


#endif

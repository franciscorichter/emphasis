#include <Rcpp.h>
#include "rinit.h"


nlopt_opt(*remp_create)(nlopt_algorithm, unsigned) = NULL;
void(*remp_destroy)(nlopt_opt) = NULL;
nlopt_result(*remp_optimize)(nlopt_opt, double *, double *) = NULL;
nlopt_result(*remp_set_min_objective)(nlopt_opt, nlopt_func, void *) = NULL;
nlopt_result(*remp_set_max_objective)(nlopt_opt, nlopt_func, void *) = NULL;
nlopt_result(*remp_set_lower_bounds)(nlopt_opt, const double *) = NULL;
nlopt_result(*remp_set_lower_bounds1)(nlopt_opt, double) = NULL;
nlopt_result(*remp_set_upper_bounds)(nlopt_opt, const double *) = NULL;
nlopt_result(*remp_set_upper_bounds1)(nlopt_opt, double) = NULL;
nlopt_result(*remp_set_xtol_rel)(nlopt_opt, double) = NULL;
nlopt_result(*remp_set_xtol_abs)(nlopt_opt, double) = NULL;


// [[Rcpp::init]]
void remphasis_init(DllInfo *dll)
{
  remp_create = (nlopt_opt(*)(nlopt_algorithm, unsigned)) R_GetCCallable("nloptr","nlopt_create");
  remp_destroy = (void(*)(nlopt_opt)) R_GetCCallable("nloptr","nlopt_destroy");
  remp_optimize = (nlopt_result(*)(nlopt_opt, double *, double *)) R_GetCCallable("nloptr","nlopt_optimize");
  remp_set_min_objective = (nlopt_result(*)(nlopt_opt, nlopt_func, void *)) R_GetCCallable("nloptr","nlopt_set_min_objective");
  remp_set_max_objective = (nlopt_result(*)(nlopt_opt, nlopt_func, void *)) R_GetCCallable("nloptr","nlopt_set_max_objective");
  remp_set_lower_bounds = (nlopt_result(*)(nlopt_opt, const double *)) R_GetCCallable("nloptr","nlopt_set_lower_bounds");
  remp_set_lower_bounds1 = (nlopt_result(*)(nlopt_opt, double)) R_GetCCallable("nloptr","nlopt_set_lower_bounds1");
  remp_set_upper_bounds = (nlopt_result(*)(nlopt_opt, const double *)) R_GetCCallable("nloptr","nlopt_set_upper_bounds");
  remp_set_upper_bounds1 = (nlopt_result(*)(nlopt_opt, double)) R_GetCCallable("nloptr","nlopt_set_upper_bounds1");
  remp_set_xtol_rel = (nlopt_result(*)(nlopt_opt, double)) R_GetCCallable("nloptr","nlopt_set_xtol_rel");
  remp_set_xtol_abs = (nlopt_result(*)(nlopt_opt, double)) R_GetCCallable("nloptr","nlopt_set_xtol_abs");
}

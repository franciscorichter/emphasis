#include "emphasis.hpp"
#include "sbplx.hpp"
#include "rinit.h"


namespace emphasis {


  sbplx_1::sbplx_1() : lower_(-std::numeric_limits<double>::max()),
    upper_(+std::numeric_limits<double>::max())
  {
    nlopt_ = remp_create(NLOPT_LN_SBPLX, 1);
    if (nullptr == nlopt_) {
      throw emphasis_error("nlopt_create failed");
    }
  }


  sbplx_1::~sbplx_1()
  {
    if (nlopt_) remp_destroy(nlopt_);
  }


  void sbplx_1::set_xtol_rel(double val)
  {
    if (NLOPT_SUCCESS > (result_ = remp_set_xtol_rel(nlopt_, val))) {
      throw emphasis::emphasis_error("nlopt_set_xtol_failed");
    }
  }


  void sbplx_1::set_lower_bounds(double val)
  {
    lower_ = val;
    if (NLOPT_SUCCESS > (result_ = remp_set_lower_bounds1(nlopt_, lower_))) {
      throw emphasis_error("nlopt_set_lower_bounds failed");
    }
  }


  void sbplx_1::set_upper_bounds(double val)
  {
    upper_ = val;
    if (NLOPT_SUCCESS > (result_ = remp_set_upper_bounds1(nlopt_, upper_))) {
      throw emphasis_error("nlopt_set_upper_bounds failed");
    }
  }


  void sbplx_1::set_min_objective(nlopt_func dx, void* fdata)
  {
    if (NLOPT_SUCCESS > (result_ = remp_set_min_objective(nlopt_, dx, fdata))) {
      throw emphasis_error("nlopt_set_min_objective failed");
    }
  }


  void sbplx_1::set_max_objective(nlopt_func dx, void* fdata)
  {
    if (NLOPT_SUCCESS > (result_ = remp_set_max_objective(nlopt_, dx, fdata))) {
      throw emphasis_error("nlopt_set_max_objective failed");
    }
  }


  double sbplx_1::optimize(double& x)
  {
    double fmin = 0.0;
    if (NLOPT_SUCCESS > (result_ = remp_optimize(nlopt_, &x, &fmin))) {
      throw emphasis_error("optimize failed");
    }
    return fmin;
  }


  nlopt_result sbplx_1::result()
  {
    return result_;
  }


  sbplx::sbplx(size_t nparams)
    : lower_(nparams, -std::numeric_limits<double>::max()),
    upper_(nparams, +std::numeric_limits<double>::max())
  {
    nlopt_ = remp_create(NLOPT_LN_SBPLX, static_cast<unsigned>(nparams));
    if (nullptr == nlopt_) {
      throw emphasis_error("nlopt_create failed");
    }
  }


  sbplx::~sbplx()
  {
    if (nlopt_) remp_destroy(nlopt_);
  }


  void sbplx::set_xtol_rel(double val)
  {
    if (NLOPT_SUCCESS > (result_ = remp_set_xtol_rel(nlopt_, val))) {
      throw emphasis::emphasis_error("nlopt_set_xtol_failed");
    }
  }


  void sbplx::set_lower_bounds(const std::vector<double>& val)
  {
    lower_ = val;
    if (NLOPT_SUCCESS > (result_ = remp_set_lower_bounds(nlopt_, lower_.data()))) {
      throw emphasis_error("nlopt_set_lower_bounds failed");
    }
  }


  void sbplx::set_upper_bounds(const std::vector<double>& val)
  {
    upper_ = val;
    if (NLOPT_SUCCESS > (result_ = remp_set_upper_bounds(nlopt_, upper_.data()))) {
      throw emphasis_error("nlopt_set_upper_bounds failed");
    }
  }

  void sbplx::set_min_objective(nlopt_func dx, void* fdata)
  {
    if (NLOPT_SUCCESS > (result_ = remp_set_min_objective(nlopt_, dx, fdata))) {
      throw emphasis_error("nlopt_set_min_objective failed");
    }
  }


  void sbplx::set_max_objective(nlopt_func dx, void* fdata)
  {
    if (NLOPT_SUCCESS > (result_ = remp_set_max_objective(nlopt_, dx, fdata))) {
      throw emphasis_error("nlopt_set_max_objective failed");
    }
  }


  double sbplx::optimize(std::vector<double>& x)
  {
    double fmin = 0.0;
    if (NLOPT_SUCCESS > (result_ = remp_optimize(nlopt_, x.data(), &fmin))) {
      throw emphasis_error("optimize failed");
    }
    return fmin;
  }


  nlopt_result sbplx::result()
  {
    return result_;
  }

}
#ifndef EMPHASIS_SBPLX_HPP_INCLUDED
#define EMPHASIS_SBPLX_HPP_INCLUDED

#include <memory>
#include <vector>
#include <nlopt.h>


namespace emphasis {


  class sbplx_1
  {
  public:
    sbplx_1(const sbplx_1&) = delete;
    sbplx_1& operator=(const sbplx_1&) = delete;

    sbplx_1();
    ~sbplx_1();
    void set_xtol_rel(double);
    void set_lower_bounds(double);
    void set_upper_bounds(double);
    void set_min_objective(nlopt_func, void*);
    void set_max_objective(nlopt_func, void*);
    double optimize(double&);
    nlopt_result result();

  private:
    nlopt_opt nlopt_ = nullptr;
    nlopt_result result_ = nlopt_result::NLOPT_FAILURE;
    double lower_, upper_;
  };


  class sbplx
  {
  public:
    sbplx(const sbplx&) = delete;
    sbplx& operator=(const sbplx&) = delete;

    sbplx(size_t nparams);
    ~sbplx();
    void set_xtol_rel(double);
    void set_lower_bounds(const std::vector<double>&);
    void set_upper_bounds(const std::vector<double>&);
    // The step the first simplex is built with, per coordinate.  Left unset,
    // NLopt takes 0.75 of the distance to the nearer bound, so a start a hair
    // inside a bound gets a step of that hair and the optimiser stops on its
    // x-tolerance before it has tried anything (see M_step).
    void set_initial_step(const std::vector<double>&);
    void set_min_objective(nlopt_func, void*);
    void set_max_objective(nlopt_func, void*);
    double optimize(std::vector<double>&);
    nlopt_result result();

  private:
    nlopt_opt nlopt_ = nullptr;
    nlopt_result result_ = nlopt_result::NLOPT_FAILURE;
    std::vector<double> lower_, upper_, step_;
  };

}

#endif

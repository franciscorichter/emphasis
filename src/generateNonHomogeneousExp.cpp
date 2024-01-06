#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector generateNonHomogeneousExpCpp(int num_variates, NumericMatrix covariates, NumericVector parameters, double start_time, double max_time) {
  NumericVector variates;
  
  // Define the exponential rate function
  auto exponentialRate = [&](double t) -> double {
    double rate = exp(parameters[0]); // Initialize with the bias term
    for (int i = 0; i < covariates.ncol(); ++i) {
      rate += parameters[i + 1] * covariates(0, i); // Apply parameters to covariates
    }
    rate = exp(rate); // Apply exponential
    return rate;
  };
  
  while (variates.size() < num_variates) {
    // Find the maximum rate over the interval
    double max_rate = exponentialRate(start_time);
    for (double t = start_time; t <= max_time; t += (max_time - start_time) / 100.0) {
      double current_rate = exponentialRate(t);
      if (current_rate > max_rate) {
        max_rate = current_rate;
      }
    }
    
    // Generate a candidate event time
    double candidate_time = start_time + rexp(1, max_rate)[0];
    
    if (candidate_time < max_time) {
      // Accept or reject the candidate time
      double true_rate = exponentialRate(candidate_time);
      if (R::runif(0, 1) < true_rate / max_rate) {
        variates.push_back(candidate_time); // Append the accepted time to the vector
      }
    }
  }
  
  return variates;
}

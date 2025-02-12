// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <boost/math/distributions/non_central_t.hpp>
#include <cmath> // For std::isnan

using namespace boost::math;
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector dnct(const NumericVector x, const double df, const NumericVector ncp) {
  
  R_xlen_t nx = x.length();
  R_xlen_t nncp = ncp.length();
  
  NumericVector y(std::max(nx, nncp));

  // Single x and multiple ncp
  if (nx == 1) {
    double single_x = x[0];

    for (R_xlen_t i = 0; i < nncp; ++i) {
      try {
        non_central_t dist(df, ncp[i]);
        y[i] = std::isnan(pdf(dist, single_x)) ? 0.0 : pdf(dist, single_x);
      } catch (...) {
        y[i] = 0.0; // **Suppresses ALL errors without slowing down**
      }
    }
  }
  // Single ncp and multiple x
  else if (nncp == 1) {
    try {
      non_central_t dist(df, ncp[0]);
      for (R_xlen_t i = 0; i < nx; ++i) {
        y[i] = std::isnan(pdf(dist, x[i])) ? 0.0 : pdf(dist, x[i]);
      }
    } catch (...) {
      std::fill(y.begin(), y.end(), 0.0); // **Assigns all zeros if a single error occurs**
    }
  }
  // Multiple x and multiple ncp with the same length
  else if (nx == nncp) {
    for (R_xlen_t i = 0; i < nx; ++i) {
      try {
        non_central_t dist(df, ncp[i]);
        y[i] = std::isnan(pdf(dist, x[i])) ? 0.0 : pdf(dist, x[i]);
      } catch (...) {
        y[i] = 0.0; // **Catches errors PER INDEX, keeping other values safe**
      }
    }
  } else {
    stop("Lengths of x and ncp must either be the same or one of them must be 1.");
  }
  
  return y;
}

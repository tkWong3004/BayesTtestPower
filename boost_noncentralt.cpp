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
  
  // Result vector
  NumericVector y(std::max(nx, nncp));
  
  if (nx == 1) {
    // Single x, multiple ncp
    double single_x = x[0];
    for (R_xlen_t i = 0; i < nncp; ++i) {
      non_central_t dist(df, ncp[i]);  // Create distribution once per ncp
      double pdf_value = pdf(dist, single_x);
      y[i] = std::isnan(pdf_value) ? 0.0 : pdf_value;  // Replace NaN with 0
    }
  } 
  else if (nncp == 1) {
    // Multiple x, single ncp
    double single_ncp = ncp[0];
    non_central_t dist(df, single_ncp);  // Create distribution once per ncp
    for (R_xlen_t i = 0; i < nx; ++i) {
      double pdf_value = pdf(dist, x[i]);
      y[i] = std::isnan(pdf_value) ? 0.0 : pdf_value;  // Replace NaN with 0
    }
  } 
  else if (nx == nncp) {
    // Multiple x and ncp with the same length
    for (R_xlen_t i = 0; i < nx; ++i) {
      non_central_t dist(df, ncp[i]);  // Create distribution once per ncp
      double pdf_value = pdf(dist, x[i]);
      y[i] = std::isnan(pdf_value) ? 0.0 : pdf_value;  // Replace NaN with 0
    }
  } 
  else {
    stop("Lengths of x and ncp must either be the same or one of them must be 1.");
  }
  
  return y;
}

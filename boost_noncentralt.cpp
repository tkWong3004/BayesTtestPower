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
      non_central_t dist(df, ncp[i]);
      double pdf_value = pdf(dist, single_x);
      y[i] = (pdf_value > 0.0) ? pdf_value : 0.0; // Replace NaN/negative values with 0
    }
  }
  // Single ncp and multiple x
  else if (nncp == 1) {
    double single_ncp = ncp[0];
    non_central_t dist(df, single_ncp);
    for (R_xlen_t i = 0; i < nx; ++i) {
      double pdf_value = pdf(dist, x[i]);
      y[i] = (pdf_value > 0.0) ? pdf_value : 0.0; // Replace NaN/negative values with 0
    }
  }
  // Multiple x and multiple ncp with the same length
  else if (nx == nncp) {
    for (R_xlen_t i = 0; i < nx; ++i) {
      non_central_t dist(df, ncp[i]);
      double pdf_value = pdf(dist, x[i]);
      y[i] = (pdf_value > 0.0) ? pdf_value : 0.0; // Replace NaN/negative values with 0
    }
  } else {
    stop("Lengths of x and ncp must either be the same or one of them must be 1.");
  }
  
  return y;
}

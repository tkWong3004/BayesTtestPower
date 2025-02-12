// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <boost/math/distributions/non_central_t.hpp>
#include <cmath> // For std::isnan

using namespace boost::math;
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector dnct(const double x, const double df, const NumericVector ncp) {
  
  R_xlen_t n = ncp.length();
  NumericVector y(n);
  
  for (R_xlen_t i = 0; i < n; ++i) {
    try {
      non_central_t dist(df, ncp[i]);
      double pdf_value = pdf(dist, x);
      y[i] = std::isnan(pdf_value) ? 0.0 : pdf_value; // Replace NaN with 0
    } catch (std::exception &e) {
      // Suppress the message to avoid flooding the console
      // Rcpp::Rcout << "Boost exception for ncp = " << ncp[i] << ": " << e.what() << std::endl;
      y[i] = 0.0; // Assign 0 on exception
    }
  }
  
  return y;
}
// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <boost/math/distributions/non_central_t.hpp> 

using namespace boost::math;
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector dnct_n(const double x, const double df, const NumericVector ncp) {
  
  R_xlen_t n = ncp.length();
  NumericVector y(n);
  
  for (R_xlen_t i = 0; i < n; ++i) {
    non_central_t dist(df, ncp[i]);
    y[i] = pdf(dist, x);
  }
  
  return y;
}
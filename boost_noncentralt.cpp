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
  
  if (nx == 1) {
    // Single x, multiple ncp
    double single_x = x[0];
    for (R_xlen_t i = 0; i < nncp; ++i) {
      try {
        non_central_t dist(df, ncp[i]);
        double pdf_value = pdf(dist, single_x);
        y[i] = std::isnan(pdf_value) ? 0.0 : pdf_value; // Replace NaN with 0
      } catch (std::exception &e) {
        // Suppress the message to avoid flooding the console
        // Rcpp::Rcout << "Boost exception for ncp = " << ncp[i] << ": " << e.what() << std::endl;
        y[i] = 0.0; // Assign 0 on exception
      }
    }
  } else if (nncp == 1) {
    // Multiple x, single ncp
    double single_ncp = ncp[0];
    non_central_t dist(df, single_ncp);
    for (R_xlen_t i = 0; i < nx; ++i) {
      try {
        double pdf_value = pdf(dist, x[i]);
        y[i] = std::isnan(pdf_value) ? 0.0 : pdf_value; // Replace NaN with 0
      } catch (std::exception &e) {
        // Suppress the message to avoid flooding the console
        // Rcpp::Rcout << "Boost exception for x = " << x[i] << ": " << e.what() << std::endl;
        y[i] = 0.0; // Assign 0 on exception
      }
    }
  } else if (nx == nncp) {
    // Multiple x and ncp with same length
    for (R_xlen_t i = 0; i < nx; ++i) {
      try {
        non_central_t dist(df, ncp[i]);
        double pdf_value = pdf(dist, x[i]);
        y[i] = std::isnan(pdf_value) ? 0.0 : pdf_value; // Replace NaN with 0
      } catch (std::exception &e) {
        // Suppress the message to avoid flooding the console
        // Rcpp::Rcout << "Boost exception for x = " << x[i] << ", ncp = " << ncp[i] << ": " << e.what() << std::endl;
        y[i] = 0.0; // Assign 0 on exception
      }
    }
  } else {
    stop("Lengths of x and ncp must either be the same or one of them must be 1.");
  }
  
  return y;
}

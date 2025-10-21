#ifndef TENSOR_FUNCTIONS_H
#define TENSOR_FUNCTIONS_H

#include <Rcpp.h>

// Declare functions here
Rcpp::NumericVector mn_mode_prod_cpp(Rcpp::NumericVector A, Rcpp::NumericVector B, int n, int m);
Rcpp::NumericVector n_mode_prod_cpp(Rcpp::NumericVector tensor, Rcpp::NumericMatrix mat, int n);

#endif

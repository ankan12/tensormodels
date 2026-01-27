#ifndef TENSOR_FUNCTIONS_H
#define TENSOR_FUNCTIONS_H

#include <Rcpp.h>

// Declare functions here
Rcpp::NumericVector nm_prod(Rcpp::NumericVector A, Rcpp::NumericVector B, int n, int m);
Rcpp::NumericVector n_prod(Rcpp::NumericVector tensor, Rcpp::NumericMatrix mat, int n);

#endif

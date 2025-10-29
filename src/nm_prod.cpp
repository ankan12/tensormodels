#include <Rcpp.h>
#include "tensor_functions.h"

using namespace Rcpp;
//' nm-mode product
//'
//' Compute the nm-mode product of a tensor with a tensor
//'
//' @param tensor_A An array representing the tensor A.
//' @param tensor_B An array representing tensor B.
//' @param n An integer specifying which mode of tensor A to multiply across.
//' @param m An integer specifying which mode of tensor B to multiply across.
//'
//' @return An array: the result of the nm-mode product between A and B.
//'
//' @examples
//' A <- matrix(c(1, 2, 3, 4), nrow = 2)
//' x <- matrix(c(5, 6), nrow = 2)
//' nm_prod(A, x, 1, 1)
//'
// [[Rcpp::export]]
NumericVector nm_prod(NumericVector A, NumericVector B, int n, int m) {
  // --- Dimensions and checks ---
  IntegerVector dimsA = A.attr("dim");
  IntegerVector dimsB = B.attr("dim");
  int ndA = dimsA.size();
  int ndB = dimsB.size();

  int In = dimsA[n-1]; // size of mode n in A
  int Jm = dimsB[m-1]; // size of mode m in B
  if (In != Jm) {
    stop("Dimension mismatch: A dim[%d]=%d but B dim[%d]=%d", n, In, m, Jm);
  }

  // --- Build result dims (drop contracted dimensions) ---
  IntegerVector res_dims;
  for (int i = 0; i < ndA; i++) {
    if (i != n-1) res_dims.push_back(dimsA[i]);
  }
  for (int j = 0; j < ndB; j++) {
    if (j != m-1) res_dims.push_back(dimsB[j]);
  }

  // --- Strides (column-major, like R) ---
  std::vector<int> stridesA(ndA);
  stridesA[0] = 1;
  for (int i = 1; i < ndA; i++) stridesA[i] = stridesA[i-1] * dimsA[i-1];

  std::vector<int> stridesB(ndB);
  stridesB[0] = 1;
  for (int j = 1; j < ndB; j++) stridesB[j] = stridesB[j-1] * dimsB[j-1];

  // --- Prepare result container ---
  int total = 1;
  for (int d : res_dims) total *= d;
  NumericVector C(total);

  double* A_ptr = REAL(A);
  double* B_ptr = REAL(B);

  // --- Iterate over result indices (flattened) ---
  int res_order = res_dims.size();
  std::vector<int> multi_idx(res_order);

  for (int flat = 0; flat < total; flat++) {
    // Decode flat index into multi-index
    int tmp = flat;
    for (int i = 0; i < res_order; i++) {
      multi_idx[i] = tmp % res_dims[i];
      tmp /= res_dims[i];
    }

    // Map multi-index into A and B indices
    std::vector<int> idxA(ndA), idxB(ndB);

    int ia = 0, ib = 0;
    for (int i = 0; i < ndA; i++) {
      if (i == n-1) continue;
      idxA[i] = multi_idx[ia++];
    }
    for (int j = 0; j < ndB; j++) {
      if (j == m-1) continue;
      idxB[j] = multi_idx[ia + ib++];
    }

    // Perform contraction along k = 0..K-1
    double sum_val = 0.0;
    for (int k = 0; k < In; k++) {
      idxA[n-1] = k;
      idxB[m-1] = k;

      // linear offset for A
      int offsetA = 0;
      for (int i = 0; i < ndA; i++) offsetA += idxA[i] * stridesA[i];
      // linear offset for B
      int offsetB = 0;
      for (int j = 0; j < ndB; j++) offsetB += idxB[j] * stridesB[j];

      sum_val += A_ptr[offsetA] * B_ptr[offsetB];
    }

    C[flat] = sum_val;
  }

  C.attr("dim") = res_dims;
  return C;
}

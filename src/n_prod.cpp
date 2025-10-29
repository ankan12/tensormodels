#include <Rcpp.h>
#include "tensor_functions.h"
//' n-mode product
//'
//' Compute the n-mode product of a tensor with a matrix.
//'
//' @param tensor An array representing the tensor.
//' @param mat A matrix of size \eqn{m \times n} to multiply with.
//' @param n  An integer specifying the mode to multiply across.
//'
//' @return An array: the result of the n-mode product.
//'
//' @examples
//' a <- array(1:3, dim = c(3, 1, 1))
//' x <- matrix(4:9, nrow = 2, ncol = 3)
//' n_prod(a, x, 1)
//'
// [[Rcpp::export]]
Rcpp::NumericVector n_prod(Rcpp::NumericVector tensor,
                                    Rcpp::NumericMatrix mat,
                                    int n) {
  using namespace Rcpp;

  // --- Step 1: Extract and validate tensor dimensions ---
  IntegerVector dimsA = tensor.attr("dim");
  int NdA = dimsA.size();
  if (n < 1 || n > NdA)
    stop("Invalid mode index: n must be between 1 and %d", NdA);

  // --- Step 2: Check dimension compatibility ---
  // The contracted dimension (mode n of A) must equal ncol(B)
  if (mat.ncol() != dimsA[n - 1]) {
    stop("Dimension mismatch: A dim[%d]=%d but B ncol=%d",
         n, dimsA[n - 1], mat.ncol());
  }

  // --- Step 3: Perform contraction along mode n ---
  // mn_mode_prod_cpp returns a tensor with dimensions:
  //   (A_except_n, B_except_m)
  // For a matrix B (m = 2), B_except_m = (mat.nrow())
  NumericVector res = nm_prod(tensor, mat, n, /*m=*/2);

  // --- Step 4: Verify the resulting tensor rank ---
  IntegerVector dimsRes = res.attr("dim");
  if (dimsRes.size() != NdA)
    stop("Unexpected result rank: got %d, expected %d.",
         dimsRes.size(), NdA);

  // --- Step 5: Build the correct permutation order ---
  // Current dims: (A₁..Aₙ₋₁, Aₙ₊₁..AₙdA, J₁)
  // Desired dims: (A₁..Aₙ₋₁, J₁, Aₙ₊₁..AₙdA)
  // → permutation = c(1:(n-1), NdA, n:(NdA-1))

  IntegerVector perm(NdA);
  int pos = 0;

  // (1) A₁..Aₙ₋₁
  for (int i = 1; i <= n - 1; ++i)
    perm[pos++] = i;

  // (2) J₁ (currently last axis, index NdA)
  perm[pos++] = NdA;

  // (3) Remaining Aₙ₊₁..AₙdA
  for (int i = n; i <= NdA - 1; ++i)
    perm[pos++] = i;

  // --- Step 6: Apply permutation and return result ---
  Function aperm("aperm");
  NumericVector out = aperm(res, perm);

  return out;
}

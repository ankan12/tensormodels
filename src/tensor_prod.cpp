#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector tensor_prod_cpp(NumericVector A, NumericVector B, bool simplify = true) {
  // --- Get dims and orders ---
  IntegerVector dimsA = A.hasAttribute("dim") ? A.attr("dim") : IntegerVector::create(A.size());
  IntegerVector dimsB = B.hasAttribute("dim") ? B.attr("dim") : IntegerVector::create(B.size());
  int ndA = dimsA.size();
  int ndB = dimsB.size();

  // --- Full concatenated dims (for indexing) ---
  IntegerVector full_dims(ndA + ndB);
  for (int i = 0; i < ndA; i++) full_dims[i] = dimsA[i];
  for (int j = 0; j < ndB; j++) full_dims[ndA + j] = dimsB[j];

  // --- Result dims (possibly simplified) ---
  IntegerVector res_dims = clone(full_dims);
  if (simplify) {
    IntegerVector simplified;
    for (int d : res_dims) if (d != 1) simplified.push_back(d);
    if (simplified.size() == 0) simplified = IntegerVector::create(1);
    res_dims = simplified;
  }

  // --- Compute total size (always from full_dims!) ---
  int total = 1;
  for (int d : full_dims) total *= d;

  // --- Precompute strides for A and B ---
  std::vector<int> stridesA(ndA);
  stridesA[0] = 1;
  for (int i = 1; i < ndA; i++) stridesA[i] = stridesA[i-1] * dimsA[i-1];

  std::vector<int> stridesB(ndB);
  stridesB[0] = 1;
  for (int j = 1; j < ndB; j++) stridesB[j] = stridesB[j-1] * dimsB[j-1];

  // --- Prepare result ---
  NumericVector C(total);
  double* A_ptr = REAL(A);
  double* B_ptr = REAL(B);

  std::vector<int> multi_idx(full_dims.size());

  // --- Loop over all entries of result tensor ---
  for (int flat = 0; flat < total; flat++) {
    // Decode flat index -> multi-index in full_dims
    int tmp = flat;
    for (int k = 0; k < full_dims.size(); k++) {
      multi_idx[k] = tmp % full_dims[k];
      tmp /= full_dims[k];
    }

    // Split into A and B indices
    std::vector<int> idxA(ndA), idxB(ndB);
    for (int i = 0; i < ndA; i++) idxA[i] = multi_idx[i];
    for (int j = 0; j < ndB; j++) idxB[j] = multi_idx[ndA + j];

    // Compute offsets
    int offsetA = 0;
    for (int i = 0; i < ndA; i++) offsetA += idxA[i] * stridesA[i];
    int offsetB = 0;
    for (int j = 0; j < ndB; j++) offsetB += idxB[j] * stridesB[j];

    C[flat] = A_ptr[offsetA] * B_ptr[offsetB];
  }

  C.attr("dim") = res_dims;
  return C;
}

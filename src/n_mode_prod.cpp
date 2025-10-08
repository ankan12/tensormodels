#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector n_mode_prod_cpp(NumericVector tensor, NumericMatrix mat, int n) {
  IntegerVector dims = tensor.attr("dim"); // get modes of tensor
  int nd = dims.size(); // order of tensor
  int mode_dim = dims[n-1]; //size of nth mode

  if (mode_dim != mat.ncol()) //check if nth mode = cols of matrix
    stop("Dimension mismatch: tensor mode %d is %d but matrix has %d columns",
         n, mode_dim, mat.ncol());

  IntegerVector res_dims = clone(dims); //build result dims
  res_dims[n-1] = mat.nrow(); //replace nth dim with rows of mat

  int total = 1;
  for (int i = 0; i < res_dims.size(); i++) {
    total *= res_dims[i]; //compute total number of elements in result
  }

  std::vector<int> strides(nd); //strides for original tensor in column-major order
  strides[0] = 1;
  for (int i = 1; i < nd; i++){
    strides[i] = strides[i-1] * dims[i-1];
  }

  NumericVector res(total); //flat vector
  double *tensor_ptr = REAL(tensor); //pointers for speed
  double *mat_ptr = REAL(mat);
  int mat_ncol = mat.ncol(); //cols of mat
  int mat_nrow = mat.nrow(); //rows of mat

  int slice_total = total / res_dims[n-1]; // number of combinations for the fixed indices

  // loop over all combinations of indices for dimensions other than n
  for (int slice = 0; slice < slice_total; slice++){
    int tmp = slice;
    int base_tensor_offset = 0;

    std::vector<int> multi_idx(nd, 0); //store multi-index except at dim n
    for (int k = 0; k < nd; k++){
      if (k == n-1) { //skip mode n
        continue;
      }
      int idx = tmp % dims[k];
      tmp /= dims[k];
      multi_idx[k] = idx;
      base_tensor_offset += idx * strides[k];
    }

    // This slice is the set of elements tensor[base_tensor_offset + j * strides[n-1]]
    // for j = 0, ..., mode_dim - 1.
    std::vector<double> tensor_slice(mode_dim); //cache tensor slice for this index
    for (int j = 0; j < mode_dim; j++){
      tensor_slice[j] = tensor_ptr[ base_tensor_offset + j * strides[n-1] ];
    }

    // Now compute the dot product for each row r in the matrix.
    for (int r = 0; r < mat_nrow; r++){
      double sum_val = 0.0;
      // Multiply the matrix row by the cached tensor slice.
      for (int j = 0; j < mode_dim; j++){
        double mval = mat_ptr[r + j * mat_nrow];
        sum_val += mval * tensor_slice[j];
      }

      // Reconstruct the full multi-index for the result.
      // Insert the current r into the n-th dimension.
      int res_index = 0;
      int multiplier = 1;
      for (int k = 0; k < nd; k++){
        int idx = (k == n-1) ? r : multi_idx[k];
        res_index += idx * multiplier;
        multiplier *= res_dims[k];
      }
      res[res_index] = sum_val;
    }
  }

  res.attr("dim") = res_dims;
  return res;
}

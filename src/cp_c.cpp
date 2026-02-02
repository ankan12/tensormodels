// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <vector>
#include <numeric>
#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;
using Eigen::Ref;

// ---------- helpers ----------
static inline std::size_t prod_dims(const std::vector<int>& d) {
  std::size_t p = 1;
  for (int x : d) p *= (std::size_t)x;
  return p;
}

// Convert linear index (column-major) to multi-index
// idx = i1 + I1*(i2 + I2*(i3 + ...))
// returns indices[0..N-1] (0-based)
static inline void lin_to_sub_colmajor(std::size_t idx,
                                       const std::vector<int>& dims,
                                       std::vector<int>& subs) {
  const int N = (int)dims.size();
  subs.resize(N);
  for (int n = 0; n < N; ++n) {
    const std::size_t In = (std::size_t)dims[n];
    subs[n] = (int)(idx % In);
    idx /= In;
  }
}

// MTTKRP for dense tensor X (stored in R column-major order):
// M = X_(mode) * KRP(factors_except_mode)
// Returns (I_mode x R)
static MatrixXd mttkrp_dense(const double* X,
                             const std::vector<int>& dims,
                             const std::vector<MatrixXd>& A,
                             int mode) {

  const int N = (int)dims.size();
  const int R = (int)A[0].cols();
  const int Im = dims[mode];
  const std::size_t P = prod_dims(dims);

  MatrixXd M = MatrixXd::Zero(Im, R);
  std::vector<int> subs;
  subs.reserve(N);

  // Stream through tensor entries
  for (std::size_t lin = 0; lin < P; ++lin) {
    const double x = X[lin];
    if (x == 0.0) continue; // harmless, can remove if you don't want branch

    lin_to_sub_colmajor(lin, dims, subs);
    const int i_mode = subs[mode];

    // For each rank component r, multiply factor rows across other modes
    // and accumulate into M(i_mode, r)
    for (int r = 0; r < R; ++r) {
      double pr = 1.0;
      for (int n = 0; n < N; ++n) {
        if (n == mode) continue;
        pr *= A[n](subs[n], r);
      }
      M(i_mode, r) += x * pr;
    }
  }

  return M;
}

// Compute ||X||_F^2 once
static double norm_frob_sq(const double* X, std::size_t P) {
  double s = 0.0;
  for (std::size_t i = 0; i < P; ++i) {
    const double v = X[i];
    s += v * v;
  }
  return s;
}

// Compute ||Xhat||_F^2 given lambda and Gram matrices Gk = A_k^T A_k
// ||Xhat||^2 = sum_{r,s} lambda_r lambda_s * Π_k Gk(r,s)
static double norm_xhat_sq(const VectorXd& lambda,
                           const std::vector<MatrixXd>& G) {
  const int N = (int)G.size();
  MatrixXd H = G[0];
  for (int k = 1; k < N; ++k) {
    H.array() *= G[k].array();
  }
  // weighted sum: sum_{r,s} (lambda_r lambda_s) * H_{r,s}
  MatrixXd W = lambda * lambda.transpose();
  return (W.array() * H.array()).sum();
}

// Inner product <X, Xhat> computed from MTTKRP at any mode m:
// <X, Xhat> = sum_r lambda_r * (A_m(:,r)^T * MTTKRP_m(:,r))
static double innerprod_x_xhat(const double* X,
                               const std::vector<int>& dims,
                               const std::vector<MatrixXd>& A,
                               const VectorXd& lambda,
                               int mode_for_mttkrp) {

  MatrixXd M = mttkrp_dense(X, dims, A, mode_for_mttkrp);
  const MatrixXd& Am = A[mode_for_mttkrp];

  double ip = 0.0;
  const int R = (int)lambda.size();
  for (int r = 0; r < R; ++r) {
    ip += lambda[r] * Am.col(r).dot(M.col(r));
  }
  return ip;
}

// ---------- main CP-ALS ----------

// [[Rcpp::export]]
Rcpp::List cp_c(Rcpp::NumericVector X,
                      int R = 2,
                      int max_iter = 2000,
                      double tol = 1e-3,
                      Rcpp::Nullable<Rcpp::List> init = R_NilValue,
                      bool verbose = false) {

  // dims
  Rcpp::IntegerVector rdims = X.attr("dim");
  if (rdims.size() < 2) Rcpp::stop("X must be an array with dim attribute.");
  const int N = rdims.size();

  std::vector<int> dims(N);
  for (int i = 0; i < N; ++i) dims[i] = rdims[i];

  for (int i = 0; i < N; ++i) {
    if (R > dims[i]) Rcpp::stop("R must be <= every tensor dimension.");
  }

  const double* xptr = REAL(X);
  const std::size_t P = prod_dims(dims);
  const double normX2 = norm_frob_sq(xptr, P);
  const double normX  = std::sqrt(normX2);

  // initialize factors A_k (I_k x R)
  std::vector<MatrixXd> A(N);
  if (init.isNotNull()) {
    Rcpp::List L(init);
    if (L.size() != N) Rcpp::stop("init must be a list of length = tensor order.");
    for (int k = 0; k < N; ++k) {
      Rcpp::NumericMatrix Ak = L[k];
      if (Ak.nrow() != dims[k] || Ak.ncol() != R) {
        Rcpp::stop("init[[%d]] must be (%d x %d).", k+1, dims[k], R);
      }
      // copy into Eigen
      A[k] = Map<MatrixXd>(REAL(Ak), Ak.nrow(), Ak.ncol());
    }
  } else {
    // simple random init (similar spirit to many CP codes)
    for (int k = 0; k < N; ++k) {
      A[k] = MatrixXd::Random(dims[k], R);
      // normalize columns
      for (int r = 0; r < R; ++r) {
        double nrm = A[k].col(r).norm();
        if (nrm == 0.0) nrm = 1.0;
        A[k].col(r) /= nrm;
      }
    }
  }

  VectorXd lambda = VectorXd::Ones(R);

  // Gram caches Gk = A_k^T A_k
  std::vector<MatrixXd> G(N);
  for (int k = 0; k < N; ++k) {
    G[k].noalias() = A[k].transpose() * A[k];
  }

  double conv = std::numeric_limits<double>::infinity();

  for (int it = 1; it <= max_iter; ++it) {

    for (int mode = 0; mode < N; ++mode) {

      // V = Hadamard product of Gram matrices excluding 'mode'
      MatrixXd V = MatrixXd::Ones(R, R);
      for (int k = 0; k < N; ++k) {
        if (k == mode) continue;
        V.array() *= G[k].array();
      }

      // MTTKRP for this mode
      MatrixXd M = mttkrp_dense(xptr, dims, A, mode);

      // Solve for A_mode: A = M * V^{-1}
      // Do it as: A^T = V^{-1} * M^T (since V is small SPD-ish)
      Eigen::LDLT<MatrixXd> ldlt(V);
      MatrixXd At = ldlt.solve(M.transpose()); // (R x I_mode)
      A[mode] = At.transpose();                // (I_mode x R)

      // Column norms -> update lambda and normalize columns
      for (int r = 0; r < R; ++r) {
        double nrm = A[mode].col(r).norm();
        if (nrm == 0.0) nrm = 1.0;
        lambda[r] = nrm;
        A[mode].col(r) /= nrm;
      }

      // Update Gram cache for this mode
      G[mode].noalias() = A[mode].transpose() * A[mode];
    }

    // Convergence without full reconstruction:
    // ||X - Xhat|| / ||X||
    const double ip   = innerprod_x_xhat(xptr, dims, A, lambda, 0);
    const double nxh2 = norm_xhat_sq(lambda, G);
    const double res2 = std::max(0.0, normX2 + nxh2 - 2.0 * ip);
    conv = std::sqrt(res2) / normX;

    if (verbose && (it % 10 == 0 || it == 1)) {
      Rcpp::Rcout << "iter " << it << " conv=" << conv << "\n";
    }

    if (conv < tol) {
      if (verbose) Rcpp::Rcout << "Converged at iter " << it << "\n";
      break;
    }

    if (it == max_iter && verbose) {
      Rcpp::Rcout << "Reached max_iter " << max_iter << "\n";
    }
  }

  // return as list(lambda=..., mats=...)
  Rcpp::List mats(N);
  for (int k = 0; k < N; ++k) {
    Rcpp::NumericMatrix Ak(dims[k], R);
    Map<MatrixXd>(REAL(Ak), dims[k], R) = A[k];
    mats[k] = Ak;
  }

  Rcpp::NumericVector lambda_out(R);
  for (int r = 0; r < R; ++r) lambda_out[r] = lambda[r];

  return Rcpp::List::create(
    Rcpp::Named("lambda") = lambda_out,
    Rcpp::Named("mats")   = mats,
    Rcpp::Named("conv")   = conv
  );
}

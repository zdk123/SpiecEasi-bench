## Which step of Armadillo's sqrtmat_sympd fails on this BLAS/LAPACK (seen with MKL).
Rcpp::cppFunction(depends="RcppArmadillo", '
Rcpp::List diag_sqrtmat(arma::mat K, double mu) {
  arma::mat K2 = K*K + 4*mu*arma::eye(K.n_rows, K.n_rows), V, KI;
  arma::vec d;
  bool dc  = arma::eig_sym(d, V, K2, "dc");
  double dmin_dc = dc ? d.min() : NA_REAL;
  bool std = arma::eig_sym(d, V, K2, "std");
  double dmin_std = std ? d.min() : NA_REAL;
  bool sq = arma::sqrtmat_sympd(KI, K2);
  return Rcpp::List::create(Rcpp::_["finite"]=K2.is_finite(),
    Rcpp::_["asym"]=arma::abs(K2 - K2.t()).max(), Rcpp::_["eig_dc"]=dc,
    Rcpp::_["min_dc"]=dmin_dc, Rcpp::_["eig_std"]=std, Rcpp::_["min_std"]=dmin_std,
    Rcpp::_["sqrtmat"]=sq);
}')
for (p in c(10, 100)) {
  set.seed(p); X <- matrix(rnorm(p * p), p); M <- (X + t(X)) / 2
  cat("p =", p, "\n"); str(diag_sqrtmat(M, 1))
}

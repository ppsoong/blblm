#include <RcppArmadillo.h>

using namespace Rcpp;

//[[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
List fast_lm(const arma::vec &y, const arma::mat &X, const arma::vec &w) {

  int n = X.n_rows, k = X.n_cols;
  arma::mat wt = diagmat(w);
  arma::mat coef = solve((X.t() * wt) * X, (X.t() * wt) * y);
  arma::colvec residuals = y - X * coef;
  arma::vec fitted_val = X * coef;
  double sigma_sq = arma::as_scalar(residuals.t() * residuals / (n - k));
  arma::colvec stderr_bar = arma::sqrt(sigma_sq * arma::diagvec((X.t() * wt * X).i()) );

  return List::create(Named("coefficients") = coef.t(),
                      Named("stderr") = stderr_bar,
                      Named("fitted_vals") = fitted_val,
                      Named("residuals") = residuals,
                      Named("rank") = arma::rank(X),
                      Named("response") = y,
                      Named("weights") = w);
}

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector sample_intC(DataFrame df, int m){

  int n = df.nrow();
  IntegerVector subs = seq(1, m);
  IntegerVector sub_samples = sample(subs, n, true);

  return sub_samples;
}
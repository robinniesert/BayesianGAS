#include <RcppArmadillo.h>
using namespace Rcpp;

inline bool isNaInt(int x){ return all(is_na(IntegerVector::create(x))); }
inline bool isNaDouble(double x){ return all(is_na(NumericVector::create(x))); }
NumericVector armaToNumVec(arma::vec x);
arma::mat numericVecToArmaMatByRef(NumericVector x, int &nRows, int &nCols);

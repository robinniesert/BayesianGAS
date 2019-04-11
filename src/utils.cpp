#include "utils.h"
#include <RcppArmadillo.h>
using namespace Rcpp;

NumericVector armaToNumVec(arma::vec x){
  return as<NumericVector>(wrap(x));
}

arma::mat numericVecToArmaMatByRef(NumericVector x, int &nRows, int &nCols){
  if (x.hasAttribute("dim")){
    NumericVector dims = x.attr("dim");
    if (dims.size() == 1){
      nRows = dims.at(0);
      nCols = 1;
    }else if (dims.size() == 2){
      nRows = dims.at(0);
      nCols = dims.at(1);
    }else{
      stop("x should not have more than 2 dimensions");
    }
  }else{
    nRows = x.size();
    nCols = 1;
  }
  arma::mat mx(x.begin(), nRows, nCols, false);
  return mx;
}

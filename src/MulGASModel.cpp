#include "MulGASModel.h"
#include "utils.h"
#include <RcppArmadillo.h>
using namespace Rcpp;

MulGASModel::MulGASModel()
{
  Omega.zeros(1);
  A.zeros(1, 1);
  B.zeros(1, 1);
}

MulGASModel::~MulGASModel()
{
}

arma::vec MulGASModel::ScaledScore(arma::vec yt, arma::vec ft)
{
  return ScoreScale(yt, ft) * Score(yt, ft);
}

arma::mat MulGASModel::ScoreScale(arma::vec yt, arma::vec ft)
{
  return arma::eye(NumF, NumF);
}

void MulGASModel::CalculateScaledScore(arma::vec yt, arma::vec ft,
                                       arma::vec &st, arma::vec &score_t)
{
  score_t = Score(yt, ft);
  st = ScoreScale(yt, ft) * score_t;
}

void MulGASModel::CalculateDerrScaledScore(arma::vec yt, arma::vec ft,
    arma::vec dft, arma::vec &dst, arma::vec &dsft)
{
  // Not split further into seperate derrivatives for scale and score,
  // because the joint derrivative implementation can often be a lot simpler.
  stop("Derrivative s is not implemented, Only implement for Density Params.");
}

double MulGASModel::LogL(NumericVector y_, RObject f1)
{
  SetParamsCheck();

  double logL;

  if (not ParamsValid()){
    logL = -INFINITY;
  }
  else{
    int i;
    int numObs;
    int numCols;
    arma::mat y = numericVecToArmaMatByRef(y_, numObs, numCols);
    arma::inplace_trans(y);
    arma::vec ft;
    arma::vec st;

    ft = as<arma::vec>(f1);
    SetTimeTInputs(y.col(0), ft);
    logL = (numObs * LogConstant()) + UpdateLogL(y.col(0), ft);
    for (i = 1; i < numObs; i++){
      st = ScaledScore(y.col(i - 1), ft);
      ft = FilterUpdate(st, ft);

      SetTimeTInputs(y.col(i), ft);
      logL += UpdateLogL(y.col(i), ft);
    }

    arma::inplace_trans(y);

    if (logL != logL){
      logL = -INFINITY;
    }
  }

  return logL;
}

NumericVector MulGASModel::GradLogL(NumericVector y, RObject f1)
{
  // SetParamsCheck();
  //
  // int i;
  // int numObs = y.size();
  // arma::vec ft;
  // arma::vec st;
  // arma::vec score_t;
  // arma::vec dsf;
  // NumericVector grad(NumParams);
  // NumericVector df(NumParams);
  // NumericVector ds(NumParams);
  // NumericVector dLogConst(NumParams);
  //
  // ft = as<arma::vec>(f1);
  // CalculateScaledScore(y[0], ft, st, score_t);
  //
  // // compute derivative components that are constant over time
  // CalculateDerrLogConstant(dLogConst);
  // grad += numObs * dLogConst;
  //
  // // compute derrivative component assosciated to the log likelihood w fixed f
  // UpdateGradLogLFixedF(y[0], ft, grad);
  //
  // for(i = 1; i < numObs; i++){
  //   // compute derrivatives w.r.t s at time t-1
  //   CalculateDerrScaledScore(y[i - 1], ft, df, ds, dsf);
  //
  //   // compute derrivatives w.r.t f at time t
  //   CalculateDerrF(y[i - 1], ft, st, ds, dsf, df);
  //
  //   // update f and s to time t
  //   ft = FilterUpdate(score_t, ft);
  //   CalculateScaledScore(y[i], ft, st, score_t);
  //
  //   // update gradient
  //   UpdateGradWScore(score_t, df, grad);
  //   UpdateGradLogLFixedF(y[i], ft, grad);
  // }

  return NumericVector(NumParams);
}

NumericVector MulGASModel::Filter(NumericVector y_, RObject f1)
{
  SetParamsCheck();

  int i;
  int numObs;
  int numCols;
  arma::mat y = numericVecToArmaMatByRef(y_, numObs, numCols);
  arma::inplace_trans(y);
  arma::vec st;
  arma::mat f(NumF, numObs);

  f.col(0) = as<arma::vec>(f1);

  for(i = 1; i < numObs; i++){
    SetTimeTInputs(y.col(i-1), f.col(i-1));
    st = ScaledScore(y.col(i-1), f.col(i-1));
    f.col(i) = FilterUpdate(st, f.col(i-1));
  }

  arma::inplace_trans(f);
  arma::inplace_trans(y);

  return as<NumericVector>(wrap(f));
}

arma::vec MulGASModel::FilterUpdate(arma::vec st, arma::vec fPrev)
{
  if (DiagCoeffMat){
    return Omega % (arma::ones(NumF) - B) + A % st + B % fPrev;
  }else{
    return Omega * (arma::eye(NumF, NumF) - B) + A * st + B * fPrev;
  }
}

// void MulGASModel::UpdateGradWScore(arma::vec score_t, NumericVector dft,
//                                    NumericVector &grad)
// {
//   grad += score_t * dft;
// }
//
// void MulGASModel::CalculateDerrF(arma::vec yt, arma::vec ft, arma::vec st, NumericVector dst,
//                                  arma::vec dsft, NumericVector &df)
// {
//   arma::vec dsOmega = dsf * df[0];
//   arma::vec dsA = dsf * df[1];
//   arma::vec dsB = dsf * df[2];
//   df[0] = (1 - B) + A * dsOmega + B * df[0];
//   df[1] = st + A * dsA + B * df[1];
//   df[2] = - Omega + ft + A * dsB + B * df[2];
//   for (int i = 3; i < NumParams; i++){
//     df[i] = DerrFDensityParamUpdate(ds[i], df[i]);
//   }
// }
//
// double MulGASModel::DerrFDensityParamUpdate(arma::vec dsParam, arma::vec dfParam)
// {
//   return A * dsParam + B * dfParam;
// }

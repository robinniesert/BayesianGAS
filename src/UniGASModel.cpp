#include "UniGASModel.h"
#include <Rcpp.h>
using namespace Rcpp;

UniGASModel::UniGASModel()
{
  Omega = NA_REAL;
  A = NA_REAL;
  B = NA_REAL;
}

UniGASModel::~UniGASModel()
{
}

double UniGASModel::ScaledScore(double yt, double ft)
{
  return ScoreScale(yt, ft) * Score(yt, ft);
}

double UniGASModel::ScoreScale(double yt, double ft)
{
  return 1.;
}

void UniGASModel::CalculateScaledScore(double yt, double ft, double &st,
                                       double &score_t)
{
  score_t = Score(yt, ft);
  st = ScoreScale(yt, ft) * score_t;
}

void UniGASModel::CalculateDerrScaledScore(double yt, double ft,
    NumericVector dft, NumericVector &dst, double &dsft)
{
  // Not split further into seperate derrivatives for scale and score_t,
  // because the joint derrivative implementation can often be a lot simpler.
  stop("Derrivative s is not implemented, Only implement for Density Params.");
}

double UniGASModel::LogL(NumericVector y, RObject f1)
{
  SetParamsCheck();

  double logL;

  if (not ParamsValid()){
    logL = -INFINITY;
  }
  else{
    int i;
    int numObs = y.size();
    double ft;
    double st;

    ft = as<double>(f1);
    logL = (numObs * LogConstant()) + UpdateLogL(y[0], ft);

    for (i = 1; i < numObs; i++){
      st = ScaledScore(y[i-1], ft);
      ft = FilterUpdate(st, ft);
      logL += UpdateLogL(y[i], ft);
    }

    if (logL != logL){
      logL = -INFINITY;
    }
  }

  return logL;
}

NumericVector UniGASModel::GradLogL(NumericVector y, RObject f1)
{
  SetParamsCheck();

  int i;
  int numObs = y.size();
  double ft;
  double st;
  double score_t;
  double dsft;
  NumericVector grad(NumParams);
  NumericVector dft(NumParams);
  NumericVector dst(NumParams);
  NumericVector dLogConst(NumParams);

  ft = as<double>(f1);
  CalculateScaledScore(y[0], ft, st, score_t);

  // compute derivative components that are constant over time
  CalculateDerrLogConstant(dLogConst);
  grad += numObs * dLogConst;

  // compute derrivative component assosciated to the log likelihood w fixed f
  UpdateGradLogLFixedF(y[0], ft, grad);

  for(i = 1; i < numObs; i++){
    // compute derrivatives w.r.t s at time t-1
    CalculateDerrScaledScore(y[i - 1], ft, dft, dst, dsft);

    // compute derrivatives w.r.t f at time t
    CalculateDerrF(y[i - 1], ft, st, dst, dsft, dft);

    // update f and s to time t
    ft = FilterUpdate(score_t, ft);
    CalculateScaledScore(y[i], ft, st, score_t);

    // update gradient
    UpdateGradWScore(score_t, dft, grad);
    UpdateGradLogLFixedF(y[i], ft, grad);
  }

  return grad;
}

NumericVector UniGASModel::Filter(NumericVector y, RObject f1)
{
  SetParamsCheck();

  int i;
  int numObs = y.size();
  double st;
  NumericVector f(y.size());

  f[0]= as<double>(f1);

  for(i = 1; i < numObs; i++){
    st = ScaledScore(y[i-1], f[i-1]);
    f[i] = FilterUpdate(st, f[i-1]);
  }
  return f;
}

void UniGASModel::UpdateGradWScore(double score_t, NumericVector dft,
    NumericVector &grad)
{
  grad += score_t * dft;
}

void UniGASModel::CalculateDerrF(double yt, double ft, double st,
    NumericVector dst, double dsft, NumericVector &dft)
{
  double dsOmega = dsft * dft[0];
  double dsA = dsft * dft[1];
  double dsB = dsft * dft[2];
  dft[0] = (1 - B) + A * dsOmega + B * dft[0];
  dft[1] = st + A * dsA + B * dft[1];
  dft[2] = - Omega + ft + A * dsB + B * dft[2];
  for (int i = 3; i < NumParams; i++){
    dft[i] = DerrFDensityParamUpdate(dst[i], dft[i]);
  }
}

double UniGASModel::DerrFDensityParamUpdate(double dsParam, double dfParam)
{
  return A * dsParam + B * dfParam;
}

double UniGASModel::FilterUpdate(double st, double fPrev)
{
  return Omega * (1 - B) + A * st + B * fPrev;
}

#include "GASModel.h"
#include "PriorStack.h"
#include <Rcpp.h>
using namespace Rcpp;

GASModel::GASModel()
{
  Params = NumericVector::create(NA_REAL);
  ParamsML = NumericVector::create(NA_REAL);
  StdsML = NumericVector::create(NA_REAL);
  NumParams = NA_INTEGER;
  LogLValML = NA_REAL;
  Omega = NA_REAL;
  A = NA_REAL;
  B = NA_REAL;

  Name = NA_STRING;
  paramsNamed = false;
}

GASModel::~GASModel()
{
}

void GASModel::SetParams(NumericVector initParams)
{
  if (initParams.size() != NumParams){
    stop("Parameter vector is not of correct length.");
  }
  Params = initParams;
}

bool GASModel::ParamsValid()
{
  return true;
}

double GASModel::ScaledScore(double y, double f)
{
  return ScoreScale(y, f) * Score(y, f);
}

double GASModel::ScoreScale(double y, double f)
{
  return 1.;
}

void GASModel::CalculateScaledScore(double y, double f, double &s,
    double &score)
{
  score = Score(y, f);
  s = ScoreScale(y, f) * score;
}
void GASModel::CalculateDerrScaledScore(double y, double f, NumericVector df,
    NumericVector &ds, double &dsf)
{
  // Not split further into seperate derrivatives for scale and score,
  // because the joint derrivative implementation can often be a lot simpler.
  stop("Derrivative s is not implemented, Only implement for Density Params.");
}

double GASModel::LogL(NumericVector y, double f1)
{
  SetParamsCheck();

  double logL;

  if (not ParamsValid()){
    logL = -INFINITY;
  }
  else{
    int i;
    int numObs = y.size();
    double f;
    double s;

    f = f1;
    logL = (numObs * LogConstant()) + UpdateLogL(y[0], f);

    for (i = 1; i < numObs; i++){
      s = ScaledScore(y[i-1], f);
      f = FilterUpdate(s, f);
      logL += UpdateLogL(y[i], f);
    }

    if (logL != logL){
      logL = -INFINITY;
    }
  }

  return logL;
}

double GASModel::LogLWPar(NumericVector params, NumericVector y, double f1)
{
  SetParams(params);
  return LogL(y, f1);
}

double GASModel::LogPosterior(NumericVector y, double f1)
{
  if (not PriorStack_.SpecsSet){
    stop("PriorStack not specified.");
  }
  return LogL(y, f1) + PriorStack_.LogPriors(Params);
}

double GASModel::LogPosteriorWPar(NumericVector params, NumericVector y,
    double f1)
{
  SetParams(params);
  return LogPosterior(y, f1);
}

NumericVector GASModel::GradLogL(NumericVector y, double f1)
{
  SetParamsCheck();

  int i;
  int numObs = y.size();
  double f;
  double s;
  double score;
  double dsf;
  NumericVector grad(NumParams);
  NumericVector df(NumParams);
  NumericVector ds(NumParams);
  NumericVector dLogConst(NumParams);

  f = f1;
  CalculateScaledScore(y[0], f, s, score);

  // compute derivative components that are constant over time
  CalculateDerrLogConstant(dLogConst);
  grad += numObs * dLogConst;

  // compute derrivative component assosciated to the log likelihood w fixed f
  UpdateGradLogLFixedF(y[0], f, grad);

  for(i = 1; i < numObs; i++){
    // compute derrivatives w.r.t s at time t-1
    CalculateDerrScaledScore(y[i - 1], f, df, ds, dsf);

    // compute derrivatives w.r.t f at time t
    CalculateDerrF(y[i - 1], f, s, ds, dsf, df);

    // update f and s to time t
    f = FilterUpdate(score, f);
    CalculateScaledScore(y[i], f, s, score);

    // update gradient
    UpdateGradWScore(score, df, grad);
    UpdateGradLogLFixedF(y[i], f, grad);
  }

  return grad;
}

NumericVector GASModel::GradLogLWPar(NumericVector params, NumericVector y,
    double f1)
{
  SetParams(params);
  return GradLogL(y, f1);
}

NumericVector GASModel::GradLogPosterior(NumericVector y, double f1)
{
  if (not PriorStack_.SpecsSet){
    stop("PriorStack not specified.");
  }
  return GradLogL(y, f1) + PriorStack_.GradLogPriors(Params);
}

NumericVector GASModel::GradLogPosteriorWPar(NumericVector params,
    NumericVector y, double f1)
{
  SetParams(params);
  return GradLogPosterior(y, f1);
}

double GASModel::PosteriorWPar(NumericVector params, NumericVector y, double f1,
    double logOffset)
{
  return ::exp(LogPosteriorWPar(params, y, f1) + logOffset);
}

NumericVector GASModel::Filter(NumericVector y, double f1)
{
  SetParamsCheck();

  int i;
  int numObs = y.size();
  double s;
  NumericVector f(y.size());

  f[0] = f1;

  for(i = 1; i < numObs; i++){
    s = ScaledScore(y[i-1], f[i-1]);
    f[i] = FilterUpdate(s, f[i-1]);
  }
  return f;
}

void GASModel::SetParamsCheck()
{
  if (not paramsNamed){
    if (is_true(any(is_na(Params)))){
      stop("Params attribute not properly set yet.");
    }
    SetParams(Params);
  }
}

double GASModel::FilterUpdate(double s, double fPrev)
{
  return Omega * (1 - B) + A * s + B * fPrev;
}

void GASModel::UpdateGradWScore(double score, NumericVector df,
    NumericVector &grad)
{
  grad += score * df;
}

void GASModel::CalculateDerrF(double y, double f, double s, NumericVector ds,
    double dsf, NumericVector &df)
{
  double dsOmega = dsf * df[0];
  double dsA = dsf * df[1];
  double dsB = dsf * df[2];
  df[0] = (1 - B) + A * dsOmega + B * df[0];
  df[1] = s + A * dsA + B * df[1];
  df[2] = - Omega + f + A * dsB + B * df[2];
  for (int i = 3; i < NumParams; i++){
    df[i] = DerrFDensityParamUpdate(ds[i], df[i]);
  }
}

double GASModel::DerrFDensityParamUpdate(double dsParam, double dfParam)
{
  return A * dsParam + B * dfParam;
}

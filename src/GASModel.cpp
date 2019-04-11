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

double GASModel::LogLWPar(NumericVector params, NumericVector y, RObject f1)
{
  SetParams(params);
  return LogL(y, f1);
}

double GASModel::LogPosterior(NumericVector y, RObject f1)
{
  if (not PriorStack_.SpecsSet){
    stop("PriorStack not specified.");
  }
  return LogL(y, f1) + PriorStack_.LogPriors(Params);
}

double GASModel::LogPosteriorWPar(NumericVector params, NumericVector y,
    RObject f1)
{
  SetParams(params);
  return LogPosterior(y, f1);
}

NumericVector GASModel::GradLogLWPar(NumericVector params, NumericVector y,
    RObject f1)
{
  SetParams(params);
  return GradLogL(y, f1);
}

NumericVector GASModel::GradLogPosterior(NumericVector y, RObject f1)
{
  if (not PriorStack_.SpecsSet){
    stop("PriorStack not specified.");
  }
  return GradLogL(y, f1) + PriorStack_.GradLogPriors(Params);
}

NumericVector GASModel::GradLogPosteriorWPar(NumericVector params,
    NumericVector y, RObject f1)
{
  SetParams(params);
  return GradLogPosterior(y, f1);
}

double GASModel::PosteriorWPar(NumericVector params, NumericVector y, RObject f1,
    double logOffset)
{
  return ::exp(LogPosteriorWPar(params, y, f1) + logOffset);
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

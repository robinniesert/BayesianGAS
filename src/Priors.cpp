#include "Priors.h"
#include <Rcpp.h>
using namespace Rcpp;


Prior::Prior()
{
  PriorParams = NumericVector::create(NA_REAL);
  Proper = NA_LOGICAL;
}

Prior::~Prior()
{
}

double Prior::LogVal(NumericVector params, NumericVector priorParams)
{
  SetParams(priorParams);
  return LogVal(priorParams);
}

ImproperUniform::ImproperUniform()
{
  lowerBound_ = R_NegInf;
  upperBound_ = R_PosInf;
  PriorParams = NumericVector::create(lowerBound_, upperBound_);
  Proper = false;
}

ImproperUniform::ImproperUniform(NumericVector priorParams)
{
  SetParams(priorParams);
  Proper = false;
}

void ImproperUniform::SetParams(NumericVector priorParams)
{
  lowerBound_ = priorParams[0];
  upperBound_ = priorParams[1];
  PriorParams = priorParams;
}

double ImproperUniform::LogVal(NumericVector params)
{
  return
    sum(ifelse((params > lowerBound_) & (params < upperBound_), 0.0, R_NegInf));
}

NumericVector ImproperUniform::GradLogVal(NumericVector params)
{
  NumericVector grad(params.size(), 0.);
  return grad;
}




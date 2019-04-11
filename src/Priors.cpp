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

Normal::Normal()
{
  mu_ = 0.;
  sigma_ = R_PosInf;
  PriorParams = NumericVector::create(mu_, sigma_);
  Proper = true;
}

Normal::Normal(NumericVector priorParams)
{
  SetParams(priorParams);
  Proper = true;
}

void Normal::SetParams(NumericVector priorParams)
{
  mu_ = priorParams[0];
  sigma_ = priorParams[1];
  PriorParams = priorParams;
}

double Normal::LogVal(NumericVector params)
{
  return sum(dnorm(params, mu_, sigma_, true));
}

NumericVector Normal::GradLogVal(NumericVector params)
{
  stop("Not implemented yet");
}


TruncatedNormal::TruncatedNormal()
{
  mu_ = 0.;
  sigma_ = R_PosInf;
  lowerBound_ = R_NegInf;
  upperBound_ = R_PosInf;
  PriorParams = NumericVector::create(mu_, sigma_);
  Proper = true;
}

TruncatedNormal::TruncatedNormal(NumericVector priorParams)
{
  SetParams(priorParams);
  Proper = true;
}

void TruncatedNormal::SetParams(NumericVector priorParams)
{
  mu_ = priorParams[0];
  sigma_ = priorParams[1];
  lowerBound_ = priorParams[2];
  upperBound_ = priorParams[3];
  PriorParams = priorParams;
}

double TruncatedNormal::LogVal(NumericVector params)
{
  NumericVector logPr = (
    dnorm(params, mu_, sigma_, true) - ::log(sigma_) -
      ::log(R::pnorm(upperBound_, mu_, sigma_, true, false) -
        R::pnorm(lowerBound_, mu_, sigma_, true, false))
  );
  return sum(
    ifelse((params > lowerBound_) & (params < upperBound_), logPr, R_NegInf));
}

NumericVector TruncatedNormal::GradLogVal(NumericVector params)
{
  stop("Not implemented yet");
}


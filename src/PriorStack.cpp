#include "Priors.h"
#include "PriorStack.h"
#include "PriorFactory.h"
#include <Rcpp.h>
using namespace Rcpp;

// [[plugins(cpp11)]]

PriorStack::PriorStack() : priorToParamIndex_(2, 2)
{
  priorStrs_ = StringVector::create(NA_STRING);
  PriorParams = List::create(NA_REAL);
  numPriors_ = NA_INTEGER;
  oneToOne_ = NA_LOGICAL;
  SpecsSet = false;
}

PriorStack::PriorStack(StringVector priorStrs, List priorParams) :
  priorToParamIndex_(2, 2)
{
  priorStrs_ = priorStrs;
  PriorParams = priorParams;
  numPriors_ = priorStrs_.size();
  if (priorParams.size() != numPriors_){
    stop("Num of priors is not equal to num of prior parameter sets.");
  }
  //SetPriors();
  oneToOne_ = true;
  SpecsSet = true;
}

PriorStack::PriorStack(StringVector priorStrs, List priorParams,
    IntegerMatrix priorToParamIndex) : PriorStack(priorStrs, priorParams)
{
  if ((priorToParamIndex.ncol() != 2)
        || (priorToParamIndex.nrow() != numPriors_)){
    stop("Specify proper format for param index, see doc.");
  }
  priorToParamIndex_ = priorToParamIndex;
  oneToOne_ = false;
}

double PriorStack::LogPriors(NumericVector params)
{
  if (not ParamsValid(params)){
    stop("Params for evaluation is not congruent with the PriorStack spec.");
  }
  double logPrior = 0.;

  for (int i = 0; i < numPriors_; i++){
    std::unique_ptr<Prior> prior = PriorFactory::BuildPrior(
      priorStrs_[i], PriorParams[i]);
    if (oneToOne_){
      logPrior += prior -> LogVal(NumericVector::create(params[i]));
    }else{
      logPrior += prior -> LogVal(
          params[Range(priorToParamIndex_.at(i, 0),
                       priorToParamIndex_.at(i, 1))]);

    }
  }
  return logPrior;
}

double PriorStack::LogPriorsWPar(NumericVector params, List priorParams)
{
  PriorParams = priorParams;
  //SetPriors();
  return LogPriors(params);
}

NumericVector PriorStack::GradLogPriors(NumericVector params)
{
  if (not ParamsValid(params)){
    stop("Params for evaluation is not congruent with the PriorStack spec.");
  }
  NumericVector grad(params.size(), 0.);

  for (int i = 0; i < numPriors_; i++){
    std::unique_ptr<Prior> prior = PriorFactory::BuildPrior(
      priorStrs_[i], PriorParams[i]);
    if (oneToOne_){
      grad[Range(i, i)] = prior -> GradLogVal(NumericVector::create(params[i]));
    }else{
      grad[Range(priorToParamIndex_.at(i, 0),
                 priorToParamIndex_.at(i, 1))] = prior -> GradLogVal(
        params[Range(priorToParamIndex_.at(i, 0),
                     priorToParamIndex_.at(i, 1))]);

    }
  }
  return grad;
}

// void PriorStack::SetPriors()
// {
//   priors.reserve(numPriors_);
//   priors.resize(numPriors_);
//   for (int i = 0; i < numPriors_; i++){
//     priors[0] = PriorFactory::BuildPrior(priorStrs_[i], PriorParams[i]);
//   }
// }

bool PriorStack::ParamsValid(NumericVector params)
{
  if (oneToOne_){
    return params.size() == numPriors_;
  }else{
    return (params.size() - 1) == priorToParamIndex_(numPriors_ - 1, 1);
  }
}

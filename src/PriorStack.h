#ifndef PriorStack_h
#define PriorStack_h 1

#include "Priors.h"
#include <RcppArmadillo.h>
using namespace Rcpp;

class PriorStack{
public:
  PriorStack();
  PriorStack(StringVector priorStrs, List priorParams);
  PriorStack(StringVector priorStrs, List priorParams,
             IntegerMatrix priorToParamIndex);

  double LogPriors(NumericVector params);
  double LogPriorsWPar(NumericVector params, List priorParams);
  NumericVector GradLogPriors(NumericVector params);

  List PriorParams;
  bool SpecsSet;

private:
  //void SetPriors();
  bool ParamsValid(NumericVector params);

  //std::vector<std::unique_ptr<Prior>> priors;
  StringVector priorStrs_;
  IntegerMatrix priorToParamIndex_;
  int numPriors_;
  bool oneToOne_;

};

#endif

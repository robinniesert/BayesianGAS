#ifndef PriorFactory_h
#define PriorFactory_h 1

#include "Priors.h"
#include <Rcpp.h>
using namespace Rcpp;

class PriorFactory
{
public:
  static std::unique_ptr<Prior> BuildPrior(String priorStr);
  static std::unique_ptr<Prior> BuildPrior(
      String priorStr, NumericVector initParams);
};

#endif

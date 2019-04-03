#include "PriorFactory.h"
#include "Priors.h"

#include <Rcpp.h>
using namespace Rcpp;

// [[plugins(cpp11)]]

std::unique_ptr<Prior> PriorFactory::BuildPrior(String priorStr)
{
  if (priorStr=="ImproperUniform") {
    return std::unique_ptr<Prior>( new ImproperUniform());
  }else{
    stop("Specify an implemented prior, see doc for available priors.");
    return NULL;
  }
}

std::unique_ptr<Prior> PriorFactory::BuildPrior(String priorStr,
    NumericVector initParams)
{
  if (priorStr=="ImproperUniform") {
    return std::unique_ptr<Prior>( new ImproperUniform(initParams));;
  }else{
    stop("Specify an implemented prior, see doc for available priors.");
    return NULL;
  }
}

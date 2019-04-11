#include "PriorFactory.h"
#include "Priors.h"

#include <Rcpp.h>
using namespace Rcpp;

// [[plugins(cpp11)]]

std::unique_ptr<Prior> PriorFactory::BuildPrior(String priorStr)
{
  if (priorStr == "ImproperUniform") {
    return std::unique_ptr<Prior>( new ImproperUniform());
  }else if (priorStr == "Normal") {
    return std::unique_ptr<Prior>( new Normal());
  }else if (priorStr == "TruncatedNormal") {
    return std::unique_ptr<Prior>( new TruncatedNormal());
  }else{
    stop("Specify an implemented prior, see doc for available priors.");
    return NULL;
  }
}

std::unique_ptr<Prior> PriorFactory::BuildPrior(String priorStr,
    NumericVector initParams)
{
  if (priorStr == "ImproperUniform") {
    return std::unique_ptr<Prior>( new ImproperUniform(initParams));
  }else if (priorStr == "Normal") {
    return std::unique_ptr<Prior>( new Normal(initParams));
  }else if (priorStr == "TruncatedNormal") {
    return std::unique_ptr<Prior>( new TruncatedNormal(initParams));
  }else{
    stop("Specify an implemented prior, see doc for available priors.");
    return NULL;
  }
}

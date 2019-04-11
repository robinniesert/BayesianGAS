#include "GASModelFactory.h"
#include "PriorStack.h"
#include <Rcpp.h>
using namespace Rcpp;

NumericMatrix ConvertVectorToMatrix(NumericVector x){
  x.attr("dim") = Dimension(1, x.size());
  return as<NumericMatrix>(x);
}

NumericVector VectorizedPosterior(RObject params, String modelStr,
    PriorStack priorStack, NumericMatrix y, RObject f1, double logOffset,
    bool log){
  int i;
  NumericMatrix params_;

  if (is<NumericVector>(params)){
    if (not Rf_isMatrix(params)){
      // Params is a vector.
      params_ = ConvertVectorToMatrix(as<NumericVector>(params));
    }else{
      params_ = as<NumericMatrix>(params);
    }
  }else{
    stop("Params of incorrect type.");
  }

  NumericVector initParams = params_(0, _);
  int numDraws = params_.nrow();
  NumericVector logPosteriorVals(numDraws);

  GASModel *model = GASModelFactory::BuildGASModelWParWPrior(
    modelStr, initParams, priorStack);

  for(i = 0; i < numDraws; i++){
    logPosteriorVals[i] = model -> LogPosteriorWPar(params_(i, _), y, f1);
  }

  logPosteriorVals = logPosteriorVals + logOffset;

  if (log){
    return logPosteriorVals;
  }else{
    return exp(logPosteriorVals);
  }
}

RCPP_EXPOSED_CLASS(PriorStack);

RCPP_MODULE(VectorizedPosterior) {
  function("VectorizedPosterior", &VectorizedPosterior,
           List::create(_["params"], _["modelStr"], _["priorStack"], _["y"],
                        _["f1"], _["logOffset"] = 0, _["log"] = true),
                        "Vectorized GAS model posterior formatted for AdMit");
}

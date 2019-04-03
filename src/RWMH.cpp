#include "GASModelFactory.h"
#include "PriorStack.h"
#include <Rcpp.h>
using namespace Rcpp;

NumericMatrix RWMH(String modelStr, PriorStack priorStack, NumericMatrix y,
    double f1, NumericVector initParams, NumericMatrix sigma, int iter,
    double stepsize, int df, bool verbose, int printIter)
{
  int i;
  int accept = 0;
  double u;
  double logPostVal;
  double logPostProp;
  NumericVector propParams;

  GASModel *model = GASModelFactory::BuildGASModelWParWPrior(
    modelStr, initParams, priorStack);
  NumericMatrix draws(iter, model -> NumParams);
  draws(0, _) = model -> Params;

  Environment mvtnormPkg = Environment::namespace_env("mvtnorm");
  Function fRmvt = mvtnormPkg["rmvt"];
  NumericMatrix eps = fRmvt(iter, Named("sigma") = sigma, Named("df") = df);
  eps = stepsize * eps;

  logPostVal = model -> LogPosterior(y, f1);

  for(i = 1; i < iter; i++){

    propParams = draws(i - 1, _) + eps(i, _);
    logPostProp = model -> LogPosteriorWPar(propParams, y, f1);

    u = R::runif(0.0, 1.0);

    if (u < exp(logPostProp - logPostVal)){
      accept++;
      draws(i, _) = propParams;
      logPostVal = logPostProp;
    }else{
      draws(i, _) = draws(i - 1, _);
    }

    if ((i % 100) == 0){
      checkUserInterrupt();
    }

    if (verbose){
      if ((i > 0) && ((i % printIter) == 0)){
        Rprintf("iter %i\n", i);
      }
    }

  }

  if (verbose){
    Rprintf("RWMH - Accept ratio is: %.3f \n", (1. * accept / iter));
  }

  return draws;
}

RCPP_EXPOSED_CLASS(PriorStack);

RCPP_MODULE(RWMH) {
  function("RWMH", &RWMH,
           List::create(_["modelStr"], _["priorStack"], _["y"],  _["f1"],
                        _["initParams"], _["sigma"], _["iter"] = 10000,
                        _["stepsize"] = 0.001, _["df"] = 3, _["verbose"] = true,
                        _["printIter"] = 1000),
           "Random walk Metropolis Hasting sampler for GAS models");
}

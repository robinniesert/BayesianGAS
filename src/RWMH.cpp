#include "GASModelFactory.h"
#include "PriorStack.h"
#include <Rcpp.h>
using namespace Rcpp;

NumericMatrix RWMH(String modelStr, PriorStack priorStack, NumericMatrix y,
    RObject f1, NumericVector initParams, NumericMatrix sigma, int iter,
    double stepsize, int df, bool verbose, int printIter, int thinning)
{
  int i;
  int draw = 0;
  int accept = 0;
  int numDraws = iter / thinning;
  double u;
  double logPostVal;
  double logPostProp;
  NumericVector propParams;
  NumericVector currentParams;

  GASModel *model = GASModelFactory::BuildGASModelWParWPrior(
    modelStr, initParams, priorStack);
  NumericMatrix draws(numDraws, model -> NumParams);
  currentParams = clone(model -> Params);
  draws(draw, _) = clone(model -> Params);

  Environment mvtnormPkg = Environment::namespace_env("mvtnorm");
  Function fRmvt = mvtnormPkg["rmvt"];
  NumericMatrix eps = fRmvt(iter, Named("sigma") = sigma, Named("df") = df);
  eps = stepsize * eps;

  logPostVal = model -> LogPosterior(y, f1);

  for(i = 1; i < iter; i++){
    propParams = currentParams + eps(i, _);
    try{
      logPostProp = model -> LogPosteriorWPar(propParams, y, f1);
    }catch(std::runtime_error &ex) {
      String e("Warning: ");
      //e += ex.what();
      //e += "\n ";
      warning(e.push_back(ex.what()).push_back("\n "));
      logPostProp = R_NegInf;
    }

    u = R::runif(0.0, 1.0);
    if (u < exp(logPostProp - logPostVal)){
      accept++;
      currentParams = clone(propParams);
      logPostVal = logPostProp;
    }

    if ((i % thinning) == 0){
      draw ++;
      draws(draw, _) = currentParams;
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
                        _["printIter"] = 1000, _["thinning"] = 1),
           "Random walk Metropolis Hasting sampler for GAS models");
}

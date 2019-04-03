#include "GASModelFactory.h"
#include "PriorStack.h"
#include <Rcpp.h>
using namespace Rcpp;

NumericMatrix GGS(String modelStr, PriorStack priorStack, NumericMatrix y,
                  double f1, NumericVector initParams, NumericMatrix grid,
                  int iter, double logOffset, bool verbose, int printIter){
  int i;
  int j;
  int g;
  int drawIdx;
  int numGridPoints = grid.nrow();
  double u;
  double remainderRatio;
  NumericVector propParams;
  NumericVector conditionalPostGrid(numGridPoints);
  NumericVector cumConditionalPostProb(numGridPoints);
  cumConditionalPostProb[0] = 0.;

  GASModel *model = GASModelFactory::BuildGASModelWParWPrior(
      modelStr, initParams, priorStack);
  NumericMatrix draws(iter, model -> NumParams);
  draws(0, _) = initParams;
  propParams = model -> Params;

  for(i = 1; i < iter; i++){
    for(j = 0; j < model -> NumParams; j++){

      for(g = 0; g < numGridPoints; g++){
        propParams[j] = grid.at(g, j);
        conditionalPostGrid[g] = model -> PosteriorWPar(
          propParams, y, f1, logOffset);

        if(conditionalPostGrid[g] != conditionalPostGrid[g]){
          Rprintf("NaNs produced.\n");
        }

        if (g > 0){
          // compute cumulative kernel using trapezoid rule
          cumConditionalPostProb[g] = (cumConditionalPostProb[g - 1] +
            (grid.at(g, j) - grid.at(g - 1, j)) *
            (conditionalPostGrid[g] + conditionalPostGrid[g - 1]) * 0.5);
        }
      }

      // invert aproximate CDF using lineair interpolation
      u = R::runif(0, cumConditionalPostProb[numGridPoints - 1]);
      drawIdx = 1;
      while (cumConditionalPostProb[drawIdx] < u){
        drawIdx++;
      }

      remainderRatio = ((u - cumConditionalPostProb[drawIdx - 1]) /
        (cumConditionalPostProb[drawIdx] - cumConditionalPostProb[drawIdx -1]));
      propParams[j] = (remainderRatio *
        (grid.at(drawIdx, j) - grid.at(drawIdx - 1, j)) +
        grid.at(drawIdx - 1, j));
    }

    draws(i, _) = propParams;

    if ((i % 100) == 0){
      checkUserInterrupt();
    }

    if (verbose){
      if ((i > 0) && ((i % printIter) == 0)){
        Rprintf("iter %i\n", i);
      }
    }

  }

  return draws;
}

RCPP_EXPOSED_CLASS(PriorStack);

RCPP_MODULE(GGS) {
  function("GGS", &GGS,
           List::create(_["modelStr"], _["priorStack"], _["y"],  _["f1"],
                        _["initParams"], _["grid"], _["iter"] = 10000,
                        _["logOffset"] = 0, _["verbose"] = true,
                        _["printIter"] = 1000),
          "Griddy Gibbs sampler for GAS models");
}

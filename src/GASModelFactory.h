#ifndef GASModelFactory_h
#define GASModelFactory_h 1

#include "PriorStack.h"
#include "GASModel.h"
#include <Rcpp.h>
using namespace Rcpp;

class GASModelFactory
{
  public:
    static GASModel* BuildGASModel(String modelStr);
    static GASModel* BuildGASModelWPar(
        String modelStr, NumericVector initParams);
    static GASModel* BuildGASModelWParWPrior(
        String modelStr, NumericVector initParams, PriorStack priorStack);
};

#endif

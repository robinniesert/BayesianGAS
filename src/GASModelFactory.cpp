#include "PriorStack.h"
#include "GASModelFactory.h"
#include "GASModel.h"
#include "BetaGenTEGARCH.h"
#include "BetaTEGARCH.h"
#include "DPMP.h"

#include <Rcpp.h>
using namespace Rcpp;

GASModel* GASModelFactory::BuildGASModel(String modelStr)
{
  if (modelStr=="BetaGenTEGARCH") {
    return new BetaGenTEGARCH();
  }else if (modelStr=="BetaTEGARCH") {
    return new BetaTEGARCH();
  }else if (modelStr=="DPMP") {
    return new DPMP();
  }else if (modelStr=="DPMP1-I") {
    return new DPMP(1, 0.);
  }else if (modelStr=="DPMP1-H") {
    return new DPMP(1, -.5);
  }else if (modelStr=="DPMP1-Inv") {
    return new DPMP(1, -1.);
  }else if (modelStr=="DPMP2-I") {
    return new DPMP(2, 0.);
  }else if (modelStr=="DPMP2-H") {
    return new DPMP(2, -.5);
  }else if (modelStr=="DPMP2-Inv") {
    return new DPMP(2, -1.);
  }else if (modelStr=="DPMP3-I") {
    return new DPMP(3, 0.);
  }else if (modelStr=="DPMP3-H") {
    return new DPMP(3, -.5);
  }else if (modelStr=="DPMP3-Inv") {
    return new DPMP(3, -1.);
  }else{
    stop("Specify an implemented model, see doc for available models.");
    return NULL;
  }
}

GASModel* GASModelFactory::BuildGASModelWPar(
    String modelStr, NumericVector initParams)
{
  if (modelStr=="BetaGenTEGARCH") {
    return new BetaGenTEGARCH(initParams);
  }else if (modelStr=="BetaTEGARCH") {
    return new BetaTEGARCH(initParams);
  }else if (modelStr=="DPMP") {
    return new DPMP(initParams);
  }else if (modelStr=="DPMP1-I") {
    return new DPMP(initParams, 1, 0.);
  }else if (modelStr=="DPMP1-H") {
    return new DPMP(initParams, 1, -.5);
  }else if (modelStr=="DPMP1-Inv") {
    return new DPMP(initParams, 1, -1.);
  }else if (modelStr=="DPMP2-I") {
    return new DPMP(initParams, 2, 0.);
  }else if (modelStr=="DPMP2-H") {
    return new DPMP(initParams, 2, -.5);
  }else if (modelStr=="DPMP2-Inv") {
    return new DPMP(initParams, 2, -1.);
  }else if (modelStr=="DPMP3-I") {
    return new DPMP(initParams, 3, 0.);
  }else if (modelStr=="DPMP3-H") {
    return new DPMP(initParams, 3, -.5);
  }else if (modelStr=="DPMP3-Inv") {
    return new DPMP(initParams, 3, -1.);
  }else{
    stop("Specify an implemented model, see doc for available models.");
    return NULL;
  }
}



GASModel* GASModelFactory::BuildGASModelWParWPrior(
    String modelStr, NumericVector initParams, PriorStack priorStack)
{
  if (modelStr=="BetaGenTEGARCH") {
    return new BetaGenTEGARCH(initParams, priorStack);
  }else if (modelStr=="BetaTEGARCH") {
    return new BetaTEGARCH(initParams, priorStack);
  }else if (modelStr=="DPMP") {
    return new DPMP(initParams, priorStack);
  }else if (modelStr=="DPMP1-I") {
    return new DPMP(initParams, priorStack, 1, 0.);
  }else if (modelStr=="DPMP1-H") {
    return new DPMP(initParams, priorStack, 1, -.5);
  }else if (modelStr=="DPMP1-Inv") {
    return new DPMP(initParams, priorStack, 1, -1.);
  }else if (modelStr=="DPMP2-I") {
    return new DPMP(initParams, priorStack, 2, 0.);
  }else if (modelStr=="DPMP2-H") {
    return new DPMP(initParams, priorStack, 2, -.5);
  }else if (modelStr=="DPMP2-Inv") {
    return new DPMP(initParams, priorStack, 2, -1.);
  }else if (modelStr=="DPMP3-I") {
    return new DPMP(initParams, priorStack, 3, 0.);
  }else if (modelStr=="DPMP3-H") {
    return new DPMP(initParams, priorStack, 3, -.5);
  }else if (modelStr=="DPMP3-Inv") {
    return new DPMP(initParams, priorStack, 3, -1.);
  }else{
    stop("Specify an implemented model, see doc for available models.");
    return NULL;
  }
}

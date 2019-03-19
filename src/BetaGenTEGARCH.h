#ifndef BetaGenTEGARCH_h
#define BetaGenTEGARCH_h 1

#include "GASModel.h"
#include <Rcpp.h>
using namespace Rcpp;

class BetaGenTEGARCH: public GASModel{
public:
  BetaGenTEGARCH(double initMu, double initOmega, double initA, double initB, 
                 double initEtaBar, double initUpsilon);
  
  virtual void NameParams(NumericVector initParams);
  virtual double LogL(NumericVector y, double f1);
  virtual NumericVector Filter(NumericVector y, double f1);
  
  double mu;
  double etaBar;
  double upsilon;
  
private:
  double LogLUpdate(double y, double f);
  double CalculateScore(double y, double f);
  bool CheckParamValidity();
};

#endif

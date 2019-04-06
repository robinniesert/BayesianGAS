#ifndef BetaTEGARCH_h
#define BetaTEGARCH_h 1

#include "GASModel.h"
#include "PriorStack.h"
#include <Rcpp.h>
using namespace Rcpp;

class BetaTEGARCH: public GASModel{
public:
  BetaTEGARCH();
  BetaTEGARCH(NumericVector initParams);
  BetaTEGARCH(NumericVector initParams, PriorStack priorStack);
  BetaTEGARCH(double initMu, double initOmega, double initA, double initB,
                 double initEtaBar);

  virtual void SetParams(NumericVector initParams);
  virtual bool ParamsValid();

  virtual double UpdateLogL(double y, double f);
  virtual double ScaledScore(double y, double f);
  virtual double ScoreScale(double y, double f);
  virtual double Score(double y, double f);
  virtual double LogConstant();

  virtual void CalculateScaledScore(
      double y, double f, double &s, double &score);
  virtual void CalculateDerrLogConstant(NumericVector &dLogConst);
  virtual void CalculateDerrScaledScore(
      double y, double f, NumericVector df, NumericVector &ds, double &dsf);
  virtual void UpdateGradLogLFixedF(double y, double f, NumericVector &grad);

  NumericVector VolFilter(NumericVector y, double f1);

  double Mu;
  double NuBar;

};

#endif

#ifndef BetaGenTEGARCH_h
#define BetaGenTEGARCH_h 1

#include "UniGASModel.h"
#include <Rcpp.h>
using namespace Rcpp;

class BetaGenTEGARCH: public UniGASModel{
public:
  BetaGenTEGARCH();
  BetaGenTEGARCH(NumericVector initParams);
  BetaGenTEGARCH(NumericVector initParams, PriorStack priorStack);
  BetaGenTEGARCH(double initOmega, double initA, double initB, double initMu,
                 double initEtaBar, double initUpsilon);

  virtual void SetParams(NumericVector initParams);
  virtual bool ParamsValid();

  virtual double UpdateLogL(double yt, double ft);
  virtual double ScaledScore(double yt, double ft);
  virtual double ScoreScale(double yt, double ft);
  virtual double Score(double yt, double ft);
  virtual double LogConstant();

  virtual void CalculateScaledScore(
      double yt, double ft, double &st, double &score_t);
  virtual void CalculateDerrLogConstant(NumericVector &dLogConst);
  virtual void CalculateDerrScaledScore(
      double yt, double ft, NumericVector dft, NumericVector &dst,
      double &dsft);
  virtual void UpdateGradLogLFixedF(double yt, double ft, NumericVector &grad);

  NumericVector VolFilter(NumericVector y, RObject f1);

  double Mu;
  double EtaBar;
  double Upsilon;

};

#endif

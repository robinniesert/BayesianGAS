#ifndef DPMP_h
#define DPMP_h 1

#include "MulGASModel.h"
#include <RcppArmadillo.h>
using namespace Rcpp;

const int kNumTransitions = 4;
const int kNumParamsOneF = 9;
const int kNumParamsTwoF = 10;
const int kNumParamsThreeF = 12;
const double kEpslion = 1e-8;

class DPMP: public MulGASModel{
public:
  DPMP(int numF = 1, double fisherInfPower = 0.);
  DPMP(NumericVector initParams, int numF = 1, double fisherInfPower = 0.);
  DPMP(NumericVector initParams, PriorStack priorStack, int numF = 1,
       double fisherInfPower = 0.);

  virtual void SetParams(NumericVector initParams);
  virtual bool ParamsValid();

  virtual double LogConstant();
  virtual double UpdateLogL(arma::vec yt, arma::vec ft);
  virtual arma::vec ScaledScore(arma::vec yt, arma::vec ft);
  virtual arma::mat ScoreScale(arma::vec yt, arma::vec ft);
  virtual arma::vec Score(arma::vec yt, arma::vec ft);

  virtual void CalculateScaledScore(
      arma::vec yt, arma::vec ft, arma::vec &st, arma::vec &score_t);
  virtual void CalculateDerrLogConstant(arma::vec &dLogConst);
  virtual void CalculateDerrScaledScore(
      arma::vec yt, arma::vec ft, arma::vec dft, arma::vec &dst,
      arma::vec &dsft);
  virtual void UpdateGradLogLFixedF(
      arma::vec yt, arma::vec ft, arma::vec &grad);
  virtual void SetTimeTInputs(arma::vec yt, arma::vec ft);

  arma::mat IntensityFilter(NumericVector y, RObject f1, bool log = true);

  void setOmega(arma::vec newOmega);
  void setA(arma::mat newA);
  void setB(arma::mat newB);
  arma::vec getOmega();
  arma::mat getA();
  arma::mat getB();

  arma::mat C;
  arma::vec W;

private:;
  arma::mat FisherInf();

  double fisherInfPower_;
  double dtau_t_;
  arma::vec logLambda_t_;
  arma::vec yt_;
  arma::vec Kt_;
};

#endif

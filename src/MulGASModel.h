#ifndef MulGASModel_h
#define MulGASModel_h 1

#include "GASModel.h"
#include <RcppArmadillo.h>
using namespace Rcpp;

class MulGASModel: public GASModel{
public:
  MulGASModel();
  virtual ~MulGASModel();

  virtual double UpdateLogL(arma::vec yt, arma::vec ft) = 0;
  virtual arma::vec ScaledScore(arma::vec yt, arma::vec ft) = 0;
  virtual arma::mat ScoreScale(arma::vec yt, arma::vec ft) = 0;
  virtual arma::vec Score(arma::vec yt, arma::vec ft) = 0;

  virtual void CalculateScaledScore(
      arma::vec yt, arma::vec ft, arma::vec &st, arma::vec &score_t) = 0;
  virtual void CalculateDerrLogConstant(arma::vec &dLogConst) = 0;
  virtual void CalculateDerrScaledScore(
      arma::vec yt, arma::vec ft, arma::vec dft, arma::vec &dst,
      arma::vec &dsf) = 0;
  virtual void UpdateGradLogLFixedF(
      arma::vec yt, arma::vec ft, arma::vec &grad) = 0;
  virtual void SetTimeTInputs(arma::vec yt, arma::vec ft) = 0;

  virtual double LogL(NumericVector y, RObject f1);
  virtual NumericVector GradLogL(NumericVector y, RObject f1);
  virtual NumericVector Filter(NumericVector y, RObject f1);

  arma::vec Omega;
  arma::mat A;
  arma::mat B;

  int NumF;
  bool DiagCoeffMat;

protected:
  // void UpdateGradWScore(arma::vec score_t, arma::vec dft, arma::vec &grad);
  // void CalculateDerrF(arma::vec yt, arma::vec ft, arma::vec st, arma::vec dst,
  //                     arma::vec dsft, arma::vec &dft);
  // double DerrFDensityParamUpdate(arma::vec dsParam, arma::vec dfParam);
  arma::vec FilterUpdate(arma::vec st, arma::vec fPrev);
};

#endif

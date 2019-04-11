#ifndef UniGASModel_h
#define UniGASModel_h 1

#include "GASModel.h"
#include <Rcpp.h>
using namespace Rcpp;

class UniGASModel: public GASModel{
public:
  UniGASModel();
  virtual ~UniGASModel();

  virtual double UpdateLogL(double yt, double ft) = 0;
  virtual double ScaledScore(double yt, double ft) = 0;
  virtual double ScoreScale(double yt, double ft) = 0;
  virtual double Score(double yt, double ft) = 0;

  virtual void CalculateScaledScore(
      double yt, double ft, double &st, double &score_t) = 0;
  virtual void CalculateDerrLogConstant(NumericVector &dLogConst) = 0;
  virtual void CalculateDerrScaledScore(
      double yt, double ft, NumericVector dft, NumericVector &dst,
      double &dsft) = 0;
  virtual void UpdateGradLogLFixedF(
      double yt, double ft, NumericVector &grad) = 0;

  virtual double LogL(NumericVector y, RObject f1);
  virtual NumericVector GradLogL(NumericVector y, RObject f1);
  virtual NumericVector Filter(NumericVector y, RObject f1);

  double Omega;
  double A;
  double B;

protected:
  void UpdateGradWScore(double score_t, NumericVector dft, NumericVector &grad);
  void CalculateDerrF(double yt, double ft, double st, NumericVector dst,
                      double dsft, NumericVector &dft);
  double DerrFDensityParamUpdate(double dsParam, double dfParam);
  double FilterUpdate(double st, double fPrev);
};

#endif

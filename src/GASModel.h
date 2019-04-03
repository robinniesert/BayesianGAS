#ifndef GASModel_h
#define GASModel_h 1

#include "PriorStack.h"
#include <Rcpp.h>
using namespace Rcpp;

class GASModel{
public:
  GASModel();
  virtual ~GASModel();

  virtual void SetParams(NumericVector initParams) = 0;
  virtual bool ParamsValid() = 0;

  virtual double UpdateLogL(double y, double f) = 0;
  virtual double ScaledScore(double y, double f) = 0;
  virtual double ScoreScale(double y, double f) = 0;
  virtual double Score(double y, double f) = 0;
  virtual double LogConstant() = 0;

  virtual void CalculateScaledScore(
      double y, double f, double &s, double &score) = 0;
  virtual void CalculateDerrLogConstant(NumericVector &dLogConst) = 0;
  virtual void CalculateDerrScaledScore(
      double y, double f, NumericVector df, NumericVector &ds, double &dsf) = 0;
  virtual void UpdateGradLogLFixedF(
      double y, double f, NumericVector &grad) = 0;

  double LogL(NumericVector y, double f1);
  double LogLWPar(NumericVector params, NumericVector y, double f1);
  double LogPosterior(NumericVector y, double f1);
  double LogPosteriorWPar(NumericVector params, NumericVector y, double f1);
  NumericVector GradLogL(NumericVector y, double f1);
  NumericVector GradLogLWPar(NumericVector params, NumericVector y, double f1);
  NumericVector GradLogPosterior(NumericVector y, double f1);
  NumericVector GradLogPosteriorWPar(NumericVector params, NumericVector y,
      double f1);
  double PosteriorWPar(NumericVector params, NumericVector y, double f1,
      double logOffset = 0);
  NumericVector Filter(NumericVector y, double f1);

  PriorStack PriorStack_;

  NumericVector Params;
  NumericVector ParamsML;
  NumericVector StdsML;

  String Name;
  int NumParams;
  double LogLValML;
  double Omega;
  double A;
  double B;

protected:
  void SetParamsCheck();

  void UpdateGradWScore(double score, NumericVector df, NumericVector &grad);
  void CalculateDerrF(double y, double f, double s, NumericVector ds,
      double dsf, NumericVector &df);
  double DerrFDensityParamUpdate(double dsParam, double dfParam);
  double FilterUpdate(double s, double fPrev);

  bool paramsNamed;
};

#endif

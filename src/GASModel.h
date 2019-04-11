#ifndef GASModel_h
#define GASModel_h 1

#define ARMA_DONT_PRINT_ERRORS
#include "PriorStack.h"
#include <RcppArmadillo.h>
using namespace Rcpp;

class GASModel{
public:
  GASModel();
  virtual ~GASModel();

  virtual void SetParams(NumericVector initParams) = 0;
  virtual bool ParamsValid() = 0;

  virtual double LogL(NumericVector y, RObject f1) = 0;
  virtual NumericVector GradLogL(NumericVector y, RObject f1) = 0;
  virtual NumericVector Filter(NumericVector y, RObject f1) = 0;
  virtual double LogConstant() = 0;

  double LogLWPar(NumericVector params, NumericVector y, RObject f1);
  double LogPosterior(NumericVector y, RObject f1);
  double LogPosteriorWPar(NumericVector params, NumericVector y, RObject f1);
  NumericVector GradLogLWPar(NumericVector params, NumericVector y, RObject f1);
  NumericVector GradLogPosterior(NumericVector y, RObject f1);
  NumericVector GradLogPosteriorWPar(NumericVector params, NumericVector y,
      RObject f1);
  double PosteriorWPar(NumericVector params, NumericVector y, RObject f1,
      double logOffset = 0);

  PriorStack PriorStack_;

  NumericVector Params;
  NumericVector ParamsML;
  NumericVector StdsML;

  String Name;
  int NumParams;
  double LogLValML;

protected:
  void SetParamsCheck();

  bool paramsNamed;
};

#endif

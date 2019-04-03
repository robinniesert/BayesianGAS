#ifndef Priors_h
#define Priors_h 1

#include <Rcpp.h>
using namespace Rcpp;

class Prior{
public:
  Prior();
  virtual ~Prior();

  virtual void SetParams(NumericVector priorParams) = 0;
  virtual double LogVal(NumericVector params) = 0;
  virtual NumericVector GradLogVal(NumericVector params) = 0;
  double LogVal(NumericVector params, NumericVector priorParams);

  NumericVector PriorParams;
  bool Proper;
};


class ImproperUniform: public Prior{
public:
  ImproperUniform();
  ImproperUniform(NumericVector priorParams);

  virtual void SetParams(NumericVector priorParams);
  virtual double LogVal(NumericVector params);
  virtual NumericVector GradLogVal(NumericVector params);

private:
  double lowerBound_;
  double upperBound_;
};


#endif

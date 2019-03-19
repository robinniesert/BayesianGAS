#ifndef GASModel_h
#define GASModel_h 1

#include <Rcpp.h>
using namespace Rcpp;

class GASModel{
public:
  GASModel();
  GASModel(NumericVector params);
  virtual ~GASModel();
  
  virtual void NameParams(NumericVector initParams) = 0;
  virtual double LogL(NumericVector y, double f1) = 0;
  virtual NumericVector Filter(NumericVector y, double f1) = 0;
  
  NumericVector params;
  NumericVector paramsML;
  NumericVector stdsML;
  int numParams;
  double logLValML;
  double omega;
  double A;
  double B;
  bool paramsNamed;

protected:
  double FilterUpdate(double score, double fPrev);
};
  
#endif

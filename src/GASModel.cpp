#include "GASModel.h"
#include <Rcpp.h>
using namespace Rcpp;

GASModel::GASModel() 
{
  params = NA_REAL;
  paramsML = NA_REAL;
  stdsML = NA_REAL;
  numParams = NA_INTEGER;
  logLValML = NA_REAL;
  omega = NA_REAL;
  A = NA_REAL;
  B = NA_REAL;
  paramsNamed = false;
}

GASModel::GASModel(NumericVector initParams)
{
  params = initParams;
  paramsML = NA_REAL;
  stdsML = NA_REAL;
  numParams = params.size();
  logLValML = NA_REAL;
  paramsNamed = false;
  omega = NA_REAL;
  A = NA_REAL;
  B = NA_REAL;
}

GASModel::~GASModel()
{
}

double GASModel::FilterUpdate(double score, double fPrev)
{
  return omega * (1 - B) + A * score + B * fPrev;
}
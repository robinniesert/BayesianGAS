#include "BetaGenTEGARCH.h"
#include <Rcpp.h>
using namespace Rcpp;

BetaGenTEGARCH::BetaGenTEGARCH(double initMu, double initOmega, double initA, 
                               double initB, double initEtaBar, 
                               double initUpsilon)
{
  mu = initMu;
  omega = initOmega;
  A = initA;
  B = initB;
  etaBar = initEtaBar;
  upsilon = initUpsilon;
  paramsNamed = true;
}

void BetaGenTEGARCH::nameParams(NumericVector initParams)
{
  mu = initParams[0];
  omega = initParams[1];
  A = initParams[2]; 
  B = initParams[3];
  etaBar = initParams[4];
  upsilon = initParams[5];
  paramsNamed = true;
}

double BetaGenTEGARCH::LogL(NumericVector y, double f1)
{
  double logL;
  
  if (not paramsNamed){
    NameParams(params);
  }
  
  bool paramsValid = CheckParamValidity();
  if (not paramsValid){
    logL = -INFINITY; 
  }
  else{
    int i;
    int numObs = y.size();
    double f;
    double K;
    double nom;
    double score;
    
    f = f1;
    K = (0.5 * upsilon * pow(etaBar , 1 / upsilon) / 
      R::beta(1 / (upsilon * etaBar), 1 / upsilon));
    logL = (numObs * log(K)) + LogLUpdate(y[0], f);
    
    for (i = 1; i < numObs; i++){
      score = CalculateScore(y[i-1], f);
      f = FilterUpdate(score, f);
      logL += LogLUpdate(y[i], f);
    }
    
    if (logL != logL || logL == INFINITY){
      logL = -INFINITY;
    }
  }
  
  return logL;
}

NumericVector BetaGenTEGARCH::Filter(NumericVector y, double f1)
{
  int i;
  int numObs = y.size();
  double score;
  NumericVector f (y.size());
  
  if (not paramsNamed){
    NameParams(params);
  }
  
  f[0] = f1;
  
  for(i = 1; i < numObs; i++){
    score = CalculateScore(y[i-1], f[i-1]);
    f[i] = FilterUpdate(score, f[i-1]);
  }
  return f;
}

double BetaGenTEGARCH::LogLUpdate(double y, double f)
{
  return -(f + (((1 / etaBar) + 1) / upsilon) * 
           log( 1 + etaBar * pow(fabs((y - mu) * exp(-f)), upsilon)));
}

double BetaGenTEGARCH::CalculateScore(double y, double f)
{
  double nom = pow(fabs(y - mu) * exp(-f), upsilon) * etaBar;
  return ((1 / etaBar) + 1) * ( nom / (nom + 1) ) - 1;
}

bool BetaGenTEGARCH::CheckParamValidity()
{
  return (fabs(B) >= 1 || etaBar <= 0 || etaBar >= 0.5 || upsilon <= 0);
}

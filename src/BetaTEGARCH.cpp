#include "BetaTEGARCH.h"
#include <Rcpp.h>
using namespace Rcpp;

// [[plugins(cpp11)]]

BetaTEGARCH::BetaTEGARCH()
{
  Mu = NA_REAL;
  NuBar = NA_REAL;
  NumParams = 5;
  Name = "BetaTEGARCH";
}

BetaTEGARCH::BetaTEGARCH(NumericVector initParams) : BetaTEGARCH()
{
  SetParams(initParams);
}


BetaTEGARCH::BetaTEGARCH(NumericVector initParams, PriorStack priorStack)
  : BetaTEGARCH(initParams)
{
  PriorStack_ = priorStack;
}

BetaTEGARCH::BetaTEGARCH(double initMu, double initOmega, double initA,
  double initB, double initNuBar) : BetaTEGARCH(
      NumericVector::create(initMu, initOmega, initA, initB, initNuBar))
{
}

void BetaTEGARCH::SetParams(NumericVector initParams)
{
  GASModel::SetParams(initParams);
  Omega = initParams[0];
  A = initParams[1];
  B = initParams[2];
  Mu = initParams[3];
  NuBar = initParams[4];
  paramsNamed = true;
}

bool BetaTEGARCH::ParamsValid()
{
  return (::fabs(B) <= 1. || NuBar >= 0. || NuBar <= 0.5);
}

double BetaTEGARCH::UpdateLogL(double y, double f)
{
  return -(f + (((1 / NuBar) + 1) / 2) *
           log( 1 + NuBar * ::pow(fabs((y - Mu) * ::exp(-f)), 2)));
}

double BetaTEGARCH::ScaledScore(double y, double f)
{
  return GASModel::ScaledScore(y, f);
}

double BetaTEGARCH::ScoreScale(double y, double f)
{
  return GASModel::ScoreScale(y, f);
}

double BetaTEGARCH::Score(double y, double f)
{
  double nom = ::pow(::fabs(y - Mu) * ::exp(-f), 2) * NuBar;
  return ((1 / NuBar) + 1) * (nom / (nom + 1)) - 1;
}

double BetaTEGARCH::LogConstant()
{
  double K = (::tgamma((1 + NuBar) / (NuBar * 2))  * ::sqrt(NuBar /  M_PI) /
              ::tgamma(1 / (NuBar * 2)));
  return ::log(K);
}

void BetaTEGARCH::CalculateScaledScore(double y, double f, double &s,
                                       double &score)
{
  GASModel::CalculateScaledScore(y, f, s, score);
}

void BetaTEGARCH::CalculateDerrLogConstant(NumericVector &dLogConst)
{
  dLogConst[4] = (
    (R::digamma(1 / (NuBar * 2)) - R::digamma((1 + NuBar) / (NuBar * 2)) +
      NuBar) / (NuBar * NuBar * 2)
  );
}

void BetaTEGARCH::UpdateGradLogLFixedF(double y, double f, NumericVector &grad)
{
  double deMeanedY = y - Mu;
  double nom = ::pow(::fabs(deMeanedY) * ::exp(-f), 2.0) * NuBar;
  double b = nom / (nom + 1);

  grad[3] += ((1 / NuBar) + 1) * b / deMeanedY;
  grad[4] += ((- ::log(1 - b) - (NuBar + 1) * b) / (2 * NuBar * NuBar));
}

void BetaTEGARCH::CalculateDerrScaledScore(double y, double f, NumericVector df,
                                           NumericVector &ds, double &dsf)
{
  double dsb = (NuBar + 1) / NuBar;
  double deMeanedY = y - Mu;
  double nom = ::pow(::fabs(deMeanedY) * ::exp(-f), 2) * NuBar;
  double b = nom / (nom + 1);
  double bOneMinB = b * (1 - b);

  double dbf = - 2 * bOneMinB;
  double dbMu = dbf * ((1 / deMeanedY) + df[3]);
  double dbNuB = (bOneMinB / NuBar) + dbf * df[4];

  // compute derrivatives w.r.t s at time t
  dsf = dsb * dbf;
  ds[3] = dsb * dbMu;
  ds[4] = dsb * dbNuB - b / (NuBar * NuBar);
}

NumericVector BetaTEGARCH::VolFilter(NumericVector y, double f1)
{
  NumericVector scales = Filter(y, f1);
  NumericVector vols = exp(scales) / ::sqrt(1 - 2 * NuBar);
  return vols;
}


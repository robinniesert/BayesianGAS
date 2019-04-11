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

BetaTEGARCH::BetaTEGARCH(double initOmega, double initA, double initB,
    double initMu, double initNuBar) : BetaTEGARCH(
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

double BetaTEGARCH::UpdateLogL(double yt, double ft)
{
  return -(ft + (((1 / NuBar) + 1) / 2) *
           log( 1 + NuBar * ::pow(fabs((yt - Mu) * ::exp(-ft)), 2)));
}

double BetaTEGARCH::ScaledScore(double yt, double ft)
{
  return UniGASModel::ScaledScore(yt, ft);
}

double BetaTEGARCH::ScoreScale(double yt, double ft)
{
  return UniGASModel::ScoreScale(yt, ft);
}

double BetaTEGARCH::Score(double yt, double ft)
{
  double nom = ::pow(::fabs(yt - Mu) * ::exp(-ft), 2) * NuBar;
  return ((1 / NuBar) + 1) * (nom / (nom + 1)) - 1;
}

double BetaTEGARCH::LogConstant()
{
  double K = (::tgamma((1 + NuBar) / (NuBar * 2))  * ::sqrt(NuBar /  M_PI) /
              ::tgamma(1 / (NuBar * 2)));
  return ::log(K);
}

void BetaTEGARCH::CalculateScaledScore(double yt, double ft, double &st,
                                       double &score_t)
{
  UniGASModel::CalculateScaledScore(yt, ft, st, score_t);
}

void BetaTEGARCH::CalculateDerrLogConstant(NumericVector &dLogConst)
{
  dLogConst[4] = (
    (R::digamma(1 / (NuBar * 2)) - R::digamma((1 + NuBar) / (NuBar * 2)) +
      NuBar) / (NuBar * NuBar * 2)
  );
}

void BetaTEGARCH::UpdateGradLogLFixedF(double yt, double ft, NumericVector &grad)
{
  double deMeanedY = yt - Mu;
  double nom = ::pow(::fabs(deMeanedY) * ::exp(-ft), 2.0) * NuBar;
  double b = nom / (nom + 1);

  grad[3] += ((1 / NuBar) + 1) * b / deMeanedY;
  grad[4] += ((- ::log(1 - b) - (NuBar + 1) * b) / (2 * NuBar * NuBar));
}

void BetaTEGARCH::CalculateDerrScaledScore(double yt, double ft,
    NumericVector dft, NumericVector &dst, double &dsft)
{
  double dsb = (NuBar + 1) / NuBar;
  double deMeanedY = yt - Mu;
  double nom = ::pow(::fabs(deMeanedY) * ::exp(-ft), 2) * NuBar;
  double b = nom / (nom + 1);
  double bOneMinB = b * (1 - b);

  double dbf = - 2 * bOneMinB;
  double dbMu = dbf * ((1 / deMeanedY) + dft[3]);
  double dbNuB = (bOneMinB / NuBar) + dbf * dft[4];

  // compute derrivatives w.r.t s at time t
  dsft = dsb * dbf;
  dst[3] = dsb * dbMu;
  dst[4] = dsb * dbNuB - b / (NuBar * NuBar);
}

NumericVector BetaTEGARCH::VolFilter(NumericVector y, RObject f1)
{
  NumericVector scales = Filter(y, f1);
  NumericVector vols = exp(scales) / ::sqrt(1 - 2 * NuBar);
  return vols;
}


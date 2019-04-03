#include "BetaGenTEGARCH.h"
#include <Rcpp.h>
using namespace Rcpp;

// [[plugins(cpp11)]]

BetaGenTEGARCH::BetaGenTEGARCH()
{
  Mu = NA_REAL;
  EtaBar = NA_REAL;
  Upsilon = NA_REAL;
  NumParams = 6;
  Name = "BetaGenTEGARCH";
}

BetaGenTEGARCH::BetaGenTEGARCH(NumericVector initParams) : BetaGenTEGARCH()
{
  SetParams(initParams);
}


BetaGenTEGARCH::BetaGenTEGARCH(NumericVector initParams, PriorStack priorStack)
  : BetaGenTEGARCH(initParams)
{
  PriorStack_ = priorStack;
}

BetaGenTEGARCH::BetaGenTEGARCH(
  double initMu, double initOmega, double initA, double initB,
  double initEtaBar, double initUpsilon) : BetaGenTEGARCH(
      NumericVector::create(
        initMu, initOmega, initA, initB, initEtaBar, initUpsilon))
{
}

void BetaGenTEGARCH::SetParams(NumericVector initParams)
{
  GASModel::SetParams(initParams);
  Omega = initParams[0];
  A = initParams[1];
  B = initParams[2];
  Mu = initParams[3];
  EtaBar = initParams[4];
  Upsilon = initParams[5];
  paramsNamed = true;
}

bool BetaGenTEGARCH::ParamsValid()
{
  return (::fabs(B) <= 1 || EtaBar >= 0 || EtaBar <= 0.5 || Upsilon >= 0);
}

double BetaGenTEGARCH::UpdateLogL(double y, double f)
{
  return -(f + (((1 / EtaBar) + 1) / Upsilon) *
           log( 1 + EtaBar * ::pow(fabs((y - Mu) * ::exp(-f)), Upsilon)));
}

double BetaGenTEGARCH::ScaledScore(double y, double f)
{
  return GASModel::ScaledScore(y, f);
}

double BetaGenTEGARCH::ScoreScale(double y, double f)
{
  return GASModel::ScoreScale(y, f);
}

double BetaGenTEGARCH::Score(double y, double f)
{
  double nom = ::pow(::fabs(y - Mu) * ::exp(-f), Upsilon) * EtaBar;
  return ((1 / EtaBar) + 1) * ( nom / (nom + 1) ) - 1;
}

double BetaGenTEGARCH::LogConstant()
{
  double K = (0.5 * Upsilon * ::pow(EtaBar, 1 / Upsilon) /
    R::beta(1 / (Upsilon * EtaBar), 1 / Upsilon));
  return ::log(K);
}

void BetaGenTEGARCH::CalculateScaledScore(double y, double f, double &s,
    double &score)
{
  GASModel::CalculateScaledScore(y, f, s, score);
}

void BetaGenTEGARCH::CalculateDerrLogConstant(NumericVector &dLogConst)
{
  dLogConst[4] = (
    (R::digamma(1 / (EtaBar * Upsilon)) -
      R::digamma((1 + EtaBar) / (EtaBar * Upsilon)) + EtaBar) /
      (Upsilon * EtaBar * EtaBar)
  );
  dLogConst[5] = (
    (EtaBar * R::digamma(1 / Upsilon) + R::digamma(1 / (EtaBar * Upsilon)) -
      (EtaBar + 1) * R::digamma((1 + EtaBar) / (EtaBar * Upsilon)) +
      (EtaBar * Upsilon) - ::log(EtaBar) * EtaBar) /
      (Upsilon * Upsilon * EtaBar)
  );
}

void BetaGenTEGARCH::UpdateGradLogLFixedF(double y, double f,
    NumericVector &grad)
{
  double etaBRatio = (EtaBar + 1) / EtaBar;
  double deMeanedY = y - Mu;
  double nom = ::pow(::fabs(deMeanedY) * ::exp(-f), Upsilon) * EtaBar;
  double b = nom / (nom + 1);

  grad[3] += etaBRatio * b / deMeanedY;
  grad[4] += ((- ::log(1 - b) - (EtaBar + 1) * b) /
    (Upsilon * EtaBar * EtaBar));
  grad[5] += (- (EtaBar + 1) *
    (log(1 - b) + Upsilon * b * (::log(fabs(deMeanedY)) - f)) /
    (Upsilon * Upsilon * EtaBar));
}

void BetaGenTEGARCH::CalculateDerrScaledScore(double y, double f,
    NumericVector df, NumericVector &ds, double &dsf)
{
  double dsb = (EtaBar + 1) / EtaBar;
  double deMeanedY = y - Mu;
  double nom = ::pow(::fabs(deMeanedY) * ::exp(-f), Upsilon) * EtaBar;
  double b = nom / (nom + 1);
  double bOneMinB = b * (1 - b);

  double dbf     = - Upsilon * bOneMinB;
  double dbMu    = dbf * ((1 / deMeanedY) + df[3]);
  double dbEtaB  = (bOneMinB / EtaBar) + dbf * df[4];
  double dbUps   = bOneMinB * (::log(::fabs(deMeanedY)) - f) + dbf * df[5];

  // compute derrivatives w.r.t s at time t
  dsf = dsb * dbf;
  ds[3] = dsb * dbMu;
  ds[4] = dsb * dbEtaB - b / (EtaBar * EtaBar);
  ds[5] = dsb * dbUps;
}

NumericVector BetaGenTEGARCH::VolFilter(NumericVector y, double f1)
{
  NumericVector scales = Filter(y, f1);
  NumericVector vols = (
    exp(scales) * ::pow(EtaBar, -(1 / Upsilon)) *
    sqrt(::tgamma(3 / Upsilon) * ::tgamma((1 / (EtaBar * Upsilon)) - (2 / Upsilon)) /
         (::tgamma(1 / (EtaBar * Upsilon)) * ::tgamma(1 / Upsilon)))
  );
  return vols;
}


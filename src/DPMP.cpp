#include "DPMP.h"
#include "utils.h"
using namespace Rcpp;

// [[plugins(cpp11)]]


DPMP::DPMP(int numF, double fisherInfPower)
{
  NumF = numF;
  fisherInfPower_ = fisherInfPower;
  switch(NumF){
  case 1:
    NumParams = kNumParamsOneF;
    Name = String("DPMP1");
    break;
  case 2:
    NumParams = kNumParamsTwoF;
    Name = String("DPMP2");
    break;
  case 3:
    NumParams = kNumParamsThreeF;
    Name = String("DPMP3");
    break;
  default:
    stop("Invalid number of factors specifed.");
    break;
  }
  if (::fabs(fisherInfPower_) < kEpslion){
    Name += "-I";
  }else if(::fabs(fisherInfPower_ + .5) < kEpslion){
    Name += "-H";
  }else if(::fabs(fisherInfPower_ + 1.) < kEpslion){
    Name += "-Inv";
  }else{
    stop("Power for Fisher Information matrix is not set to a valid value.");
  }
  C.zeros(kNumTransitions, NumF);
  W.zeros(kNumTransitions);
  DiagCoeffMat = true;
}

DPMP::DPMP(NumericVector initParams, int numF, double fisherInfPower)
    : DPMP(numF, fisherInfPower)
{
  NumParams = initParams.size();
  SetParams(initParams);
}


DPMP::DPMP(NumericVector initParams, PriorStack priorStack, int numF,
    double fisherInfPower) : DPMP(initParams, numF, fisherInfPower)
{
  PriorStack_ = priorStack;
}

void DPMP::SetParams(NumericVector initParams)
{
  if (isNaInt(NumParams)){
    NumParams = initParams.size();
  }
  GASModel::SetParams(initParams);

  int tmpIdx = 2 * NumF;
  switch(NumF){
  case 1:
    if (NumParams != kNumParamsOneF){
      stop("Supplied number of parameters incorrect for number of factors.");
    }
    C = arma::ones(kNumTransitions, NumF);
    C.head_rows(3) = arma::vec(initParams[Range(tmpIdx, tmpIdx + 2)]);
    W = arma::mat(initParams[Range(tmpIdx + 3, tmpIdx + 6)]);
    break;
  case 2:
    if (NumParams != kNumParamsTwoF){
      stop("Supplied number of parameters incorrect for number of factors.");
    }
    C = arma::zeros(kNumTransitions, NumF);
    C(0, 0) = initParams[tmpIdx];
    C(1, 0) = initParams[tmpIdx + 1];
    C(2, 1) = 1.;
    C(3, 0) = 1.;
    W = arma::mat(initParams[Range(tmpIdx + 2, tmpIdx + 5)]);
    break;
  case 3:
    if (NumParams != kNumParamsThreeF){
      stop("Supplied number of parameters incorrect for number of factors.");
    }
    C = arma::zeros(kNumTransitions, NumF);
    C(0, 0) = 1.;
    C(1, 0) = initParams[tmpIdx];
    C(1, 2) = initParams[tmpIdx + 1];
    C(2, 1) = 1.;
    C(3, 2) = 1.;
    W = arma::mat(initParams[Range(tmpIdx + 2, tmpIdx + 5)]);
    break;
  default:
    stop("Invalid number of factors specifed.");
    break;
  }
  Omega = arma::zeros(NumF);
  A = arma::mat(initParams[Range(0, NumF - 1)]);
  B = arma::mat(initParams[Range(NumF, tmpIdx - 1)]);
  paramsNamed = true;
}

bool DPMP::ParamsValid()
{
  return arma::all(arma::vectorise(arma::abs(B)) < 1.);
}

double DPMP::UpdateLogL(arma::vec yt, arma::vec ft)
{
  return (arma::dot(yt_, logLambda_t_) -
          dtau_t_ * arma::dot(Kt_, arma::exp(logLambda_t_)));
}

arma::vec DPMP::ScaledScore(arma::vec yt, arma::vec ft)
{
  if (::fabs(fisherInfPower_) < kEpslion){
    return Score(yt, ft);
  }else if(NumF == 1){
    return arma::pow(FisherInf(), fisherInfPower_) * Score(yt, ft);
  }else{
    if(::fabs(fisherInfPower_ + .5) < kEpslion){
      return arma::solve(arma::trimatu(arma::chol(FisherInf())), Score(yt, ft));
    }else if(::fabs(fisherInfPower_ + 1.) < kEpslion){
      return arma::solve(FisherInf(), Score(yt, ft), arma::solve_opts::fast);
    }else{
      stop("Power for Fisher Information matrix is not set to a valid value.");
    }
  }
}

arma::mat DPMP::ScoreScale(arma::vec yt, arma::vec ft)
{
  stop("Score scale is not implemented - calculation is integrated into scaled "
         " score utilizing speed ups with solve");
}

arma::vec DPMP::Score(arma::vec yt, arma::vec ft)
{
  return C.t() * (yt_ - dtau_t_ * Kt_ % arma::exp(logLambda_t_));
}

double DPMP::LogConstant()
{
  return 0.;
}

void DPMP::CalculateScaledScore(arma::vec yt, arma::vec ft, arma::vec &st,
    arma::vec &score_t)
{
  score_t = Score(yt, ft);
  if (::fabs(fisherInfPower_) < kEpslion){
    st = score_t;
  }else if(NumF == 1){
    arma::pow(FisherInf(), fisherInfPower_) * score_t;
  }else{
    if(::fabs(fisherInfPower_ + .5) < kEpslion){
      st = arma::solve(arma::trimatu(arma::chol(FisherInf())), score_t);
    }else if(::fabs(fisherInfPower_ + 1.) < kEpslion){
      st = arma::solve(FisherInf(), score_t, arma::solve_opts::fast);
    }else{
      stop("Power for Fisher Information matrix is not set to a valid value.");
    }
  }
}

void DPMP::CalculateDerrLogConstant(arma::vec &dLogConst)
{
  stop("Gradient not implemented for DPMP models.");
}

void DPMP::CalculateDerrScaledScore(arma::vec yt, arma::vec ft, arma::vec dft,
    arma::vec &dst,arma::vec &dsft)
{
  stop("Gradient not implemented for DPMP models.");
}

void DPMP::UpdateGradLogLFixedF(arma::vec yt, arma::vec ft, arma::vec &grad)
{
  stop("Gradient not implemented for DPMP models.");
}

void DPMP::SetTimeTInputs(arma::vec yt, arma::vec ft)
{
  logLambda_t_ = W + C * ft;
  yt_ = yt.head(kNumTransitions);
  dtau_t_ = yt(kNumTransitions);
  Kt_ = yt.tail(kNumTransitions);
}

arma::mat DPMP::IntensityFilter(NumericVector y, RObject f1, bool log)
{
  NumericVector dims = y.attr("dim");
  arma::mat f = Filter(y, f1);
  f.set_size(dims[0], NumF);
  arma::mat logLambda = C * f.t();
  logLambda.each_col() += W;
  if (log){
    return (logLambda).t();
  }else{
    return (arma::exp(logLambda)).t();
  }
}

void DPMP::setOmega(arma::vec newOmega)
{
  Omega = newOmega;
}

void DPMP::setA(arma::mat newA)
{
  A = newA;
}

void DPMP::setB(arma::mat newB)
{
  B = newB;
}


arma::vec DPMP::getOmega()
{
  return Omega;
}

arma::mat DPMP::getA()
{
  return A;
}

arma::mat DPMP::getB()
{
  return B;
}

arma::mat DPMP::FisherInf()
{
  return (C.t() * arma::diagmat(Kt_ % arma::exp(logLambda_t_) /
          dot(Kt_, arma::exp(logLambda_t_))) * C);
}


#include <RcppArmadillo.h>
#include "GASModelFactory.h"
#include "PriorStack.h"
using namespace Rcpp;
using namespace arma;

// [[plugins(cpp11)]]

NumericVector armaToNumVec(vec x){
  return as<NumericVector>(wrap(x));
}

vec HMCdraw(vec currentParams, vec currentMomentum,
    std::function<double(vec)> potentialEnergy,
    std::function<vec(vec)> gradPotentialEnergy,
    double stepsize, int leapFrogSteps, mat &mass, mat &invMass, vec &lb,
    vec &ub, int &accept, ivec &hits, bool &stuck){

  int l;
  int count = 0;
  bool constrained;
  vec params = currentParams;
  vec momentum = currentMomentum;

  // half step momentum at start
  momentum -= .5 * stepsize * gradPotentialEnergy(params);

  // perform L leapfrog integration steps
  for (l = 0; l < leapFrogSteps; l++){

    // one step params
    params += stepsize * invMass * momentum;

    // check constraints until all are met
    constrained = true;
    while (constrained){
      constrained = false;
      if (any(params <= lb)){
        uvec lbCrossIds = find(params <= lb);
        momentum.elem(lbCrossIds) *= - 1.;
        params.elem(lbCrossIds) = (lb.elem(lbCrossIds) +
          (lb.elem(lbCrossIds) - params.elem(lbCrossIds)));

        hits.elem(lbCrossIds) += 1;
        constrained = true;
      }
      if (any(params >= ub)){
        uvec ubCrossIds = find(params >= ub);
        momentum.elem(ubCrossIds) *= - 1.;
        params.elem(ubCrossIds) = (ub.elem(ubCrossIds) -
          (params.elem(ubCrossIds) - ub.elem(ubCrossIds)));

        hits.elem(ubCrossIds) += 1;
        constrained = true;
      }

      // check if stuck
      count++;
      if(count > 20){
        stuck = true;
        break;
      }
    }
    if (stuck){
      break;
    }

    // reset count
    count = 0;

    if(l < (leapFrogSteps - 1)){
      // one step momentum
      momentum -= stepsize * gradPotentialEnergy(params);
    }
  }

  // one more half step for momentum at end
  momentum -= .5 * stepsize * gradPotentialEnergy(params);

  // negate momentum to make proposal symmetric
  momentum = -momentum;

  // Evaluate current and proposed Hamiltonians
  double kineticEnergyCurr = .5 * as_scalar(
    trans(currentMomentum) * invMass * currentMomentum);
  double hamCurr = potentialEnergy(currentParams) + kineticEnergyCurr;
  double kineticEnergyProp = .5 * as_scalar(
    trans(momentum) * invMass * momentum);
  double hamProp = potentialEnergy(params) + kineticEnergyProp;

  // reject if NAs produced or got stuck
  if (hamProp != hamProp || stuck){
    warning("error: Hamiltonian NaN or HMC stuck between bounds.");
    params = currentParams;
  }else{
    // Accept or reject proposal
    double u = R::runif(0., 1.);
    if (u < ::exp(hamCurr - hamProp)){
      accept++;
    }else{
      params = currentParams;
    }
  }

  return params;
}

mat HMC(String modelStr, PriorStack priorStack, NumericMatrix y,
    double f1, NumericVector initParams, int iter, mat mass, double stepsize,
    double integrationTime, vec lb, vec ub, double stepReductionFactor,
    bool verbose, int printIter){
  int i;
  int accept = 0;
  int stuckCount = 0;
  bool stuck = false;
  int leapFrogSteps = ::floor(integrationTime / stepsize);

  GASModel *model = GASModelFactory::BuildGASModelWParWPrior(
    modelStr, initParams, priorStack);
  ivec hits = zeros<ivec>(model -> NumParams);
  mat draws(model -> NumParams, iter, fill::zeros);
  draws.col(0) = as<vec>(model -> Params);

  mat momenta = mvnrnd(zeros(model -> NumParams), mass, iter - 1);
  mat invMass = inv(mass);

  auto potentialEnergy = [&](vec params){
    return -model -> LogPosteriorWPar(armaToNumVec(params), y, f1);
  };
  auto gradPotentialEnergy = [&](vec params){
    NumericVector tmp = - model -> GradLogPosteriorWPar(armaToNumVec(params), y, f1);
    return as<vec>(tmp);
  };

  for(i = 1; i < iter; i++){
    // perform one draw of HMC
    draws.col(i) = HMCdraw(draws.col(i - 1), momenta.col(i - 1),
              potentialEnergy, gradPotentialEnergy, stepsize, leapFrogSteps, mass,
              invMass, lb, ub, accept, hits, stuck);

    // check if stuck
    if(stuck){
      stuck = false;
      stuckCount++;
      stepsize *= stepReductionFactor;
      if (::fabs(stepReductionFactor - 1.) > 1e-10){
        Rprintf("Reducing stepsize to %.5f\n", stepsize);
      }
      if (stuckCount >= 5)
      {
        stop("HMC is not working."
               " Try changing settings, but posterior might be ill-condidtioned");
      }
    }

    if ((i % 100) == 0){
      checkUserInterrupt();
    }

    if (verbose){
      if ((i > 0) && ((i % printIter) == 0)){
        Rprintf("iter %i\n", i);
      }
    }

  }

  if (verbose){
    Rprintf("HMC - Accept ratio is: %.3f \n", (1. * accept / iter));
  }

  return trans(draws);
}

RCPP_EXPOSED_CLASS(PriorStack);

RCPP_MODULE(HMC) {
  function("HMC", &HMC,
           List::create(_["modelStr"], _["priorStack"], _["y"],  _["f1"],
                        _["initParams"], _["iter"], _["mass"], _["stepsize"],
                        _["integrationTime"], _["lb"], _["ub"],
                        _["stepReductionFactor"] = 1., _["verbose"] = true,
                        _["printIter"] = 1000),
           "Hamiltonian Monte Carlo sampler for GAS models");
}


#include "GASModel.h"
#include "UniGASModel.h"
#include "GASModelFactory.h"
#include "BetaGenTEGARCH.h"
#include "BetaTEGARCH.h"
#include "DPMP.h"
#include <RcppArmadillo.h>
using namespace Rcpp;

RCPP_EXPOSED_CLASS(PriorStack);

RCPP_MODULE(GASModel){

  class_<GASModel>("GASModel")
    .factory<String>(&GASModelFactory::BuildGASModel)
    .factory<String, NumericVector>(&GASModelFactory::BuildGASModelWPar)
    .factory<String, NumericVector, PriorStack>(
        &GASModelFactory::BuildGASModelWParWPrior)
    .method("SetParams", &GASModel::SetParams)
    .method("LogL", &GASModel::LogL)
    .method("LogLWPar", &GASModel::LogLWPar)
    .method("LogPosteriorWPar", &GASModel::LogPosteriorWPar)
    .method("GradLogLWPar", &GASModel::GradLogLWPar)
    .method("Filter", &GASModel::Filter)
    .field("PriorStack", &GASModel::PriorStack_)
    .field("Params", &GASModel::Params)
    .field("ParamsML", &GASModel::ParamsML)
    .field("StdsML", &GASModel::StdsML)
    .field("LogLValML", &GASModel::LogLValML)
    .field_readonly("NumParams", &GASModel::NumParams)
    .field_readonly("Name", &GASModel::Name)
  ;

  class_<BetaGenTEGARCH>("BetaGenTEGARCH")
    .derives<GASModel>("GASModel")
    .constructor()
    .constructor<NumericVector>()
    .constructor<double, double, double, double, double, double>()
    .method("VolFilter", &BetaGenTEGARCH::VolFilter)
    .field("Mu", &BetaGenTEGARCH::Mu)
    .field("EtaBar", &BetaGenTEGARCH::EtaBar)
    .field("Upsilon", &BetaGenTEGARCH::Upsilon)
  ;

  class_<BetaTEGARCH>("BetaTEGARCH")
    .derives<GASModel>("GASModel")
    .constructor()
    .constructor<NumericVector>()
    .constructor<NumericVector, PriorStack>()
    .constructor<double, double, double, double, double>()
    .method("VolFilter", &BetaTEGARCH::VolFilter)
    .field("Mu", &BetaTEGARCH::Mu)
    .field("NuBar", &BetaTEGARCH::NuBar)
  ;

  class_<DPMP>("DPMP")
    .derives<GASModel>("GASModel")
    .constructor()
    .constructor<int, double>()
    .constructor<NumericVector, int, double>()
    .constructor<NumericVector, PriorStack, int, double>()
    .method("IntensityFilter", &DPMP::IntensityFilter)
    .method("Simulate", &DPMP::Simulate)
    .method("setOmega", &DPMP::setOmega)
    .method("setA", &DPMP::setA)
    .method("setB", &DPMP::setB)
    .method("getOmega", &DPMP::getOmega)
    .method("getA", &DPMP::getA)
    .method("getB", &DPMP::getB)
    .field("C", &DPMP::C)
    .field("W", &DPMP::W)
  ;
};

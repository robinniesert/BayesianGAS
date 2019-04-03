
#include "GASModel.h"
#include "GASModelFactory.h"
#include "BetaGenTEGARCH.h"
#include <Rcpp.h>
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
    .method("GradLogLWPar", &GASModel::GradLogLWPar)
    .method("Filter", &GASModel::Filter)
    .field("PriorStack", &GASModel::PriorStack_)
    .field("Params", &GASModel::Params)
    .field("Omega", &GASModel::Omega)
    .field("A", &GASModel::A)
    .field("B", &GASModel::B)
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
    .constructor<double,double,double,double,double,double>()
    .method("VolFilter", &BetaGenTEGARCH::VolFilter)
    .field("Mu", &BetaGenTEGARCH::Mu)
    .field("EtaBar", &BetaGenTEGARCH::EtaBar)
    .field("Upsilon", &BetaGenTEGARCH::Upsilon)
  ;

};

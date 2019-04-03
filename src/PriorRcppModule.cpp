
#include "PriorStack.h"
#include <Rcpp.h>
using namespace Rcpp;


RCPP_MODULE(Prior){

  class_<PriorStack>("PriorStack")
    .constructor<StringVector, List>()
    .constructor<StringVector, List, IntegerMatrix>()
    .method("LogPriors", &PriorStack::LogPriors)
    .method("LogPriorsWPar", &PriorStack::LogPriorsWPar)
    .method("GradLogPriors", &PriorStack::GradLogPriors)
    .field("PriorParams", &PriorStack::PriorParams)
  ;

};

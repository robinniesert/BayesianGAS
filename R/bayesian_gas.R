#' @name GASModelIndex
#'
#' @aliases
#' GASModel       Rcpp_GASModel-class       Rcpp_GASModel
#' BetaGenTEGARCH Rcpp_BetaGenTEGARCH-class Rcpp_BetaGenTEGARCH
#' BetaTEGARCH    Rcpp_BetaTEGARCH-class    Rcpp_BetaTEGARCH
#' DPMP           Rcpp_DPMP-class           Rcpp_DPMP
#'
#'
#'
#' @title Generalized Autoregressive Score (GAS) models
#'
#' @description
#' The GAS framework provides a lot of structure, which all models within the
#' framework share. The implementations of the models in this package reflect
#' this structure through class ineheritance. The parent \code{GASmodel} class
#' implements a significant share of the functionality and also functions as an
#' interface. The specific GAS model classes implement the details of the model
#' specific probability density and other desired specialized functionality.
#'
#' @section Usage:
#' \preformatted{
#' m <- new(GASModel, modelStr)
#' m <- new(GASModel, modelStr, initParams)
#' m <- new(GASModel, modelStr, initParams, priorStack)
#'
#' m <- new(BetaTEGARCH, initParams)
#' m <- new(BetaTEGARCH, initParams, priorStack)
#'
#' m$SetParams(params)
#'
#' m$LogL(y, f1)
#' m$LogLWPar(params, y, f1)
#' m$LogPosteriorWPar(params, y, f1)
#' m$GradLogLWPar(params, y, f1)
#'
#' m$Filter(y, f1)
#' }
#'
#' @section Details:
#'
#' \code{new(GASModel, modelStr)},
#' \code{new(GASModel, modelStr, initParams)},
#' \code{new(GASModel, modelStr, initParams, priorStack)}.
#' Create a new GASModel instance of type \code{GASModel}. The
#' instance has access to the functionality associated with the model type
#' specified in modelStr. The string modelStr has to be one of the models
#' included in the supported model list (see \code{\link{BayesianGAS}}). Optionally
#' supply a parameter vector \code{initParams} and a prior specification
#' \code{priorStack} (see \code{\link{PriorStack}})
#'
#' \code{new(Class, initParams)},
#' \code{new(Class, initParams, priorStack)}
#' Create a new GASModel instance of type \code{Class} where \code{Class} is
#' one of the models included in the supported model list. Optionally
#' supply a parameter vector \code{initParams} and a prior specification
#' \code{priorStack}
#'
#' \code{$SetParams(params)}
#' Sets the model specific Parameters with the numeric vector \code{params}.
#'
#' \code{$LogL(y, f1)}
#' Computes the models log likelihood for data \code{y} and the initial state
#' for the time-varying parameters set to \code{f1}.
#'
#' \code{$LogLWPar(y, f1, params)}
#' Same as \code{$LogL}, but is combined with a call to \code{SetParams} to
#' update the models parameters with \code{params}.
#'
#' \code{$LogPosteriorWPar(y, f1, params)}
#' Computes the models log posterior, which is defined as the sum of log
#' likelihood plus log prior for GAS models. Like for \code{$LogLWPar} the
#' model parameters are updated with \code{params}.
#'
#' \code{$GradLogLWPar(params, y, f1)}
#' Computes the gradient of the models log likelhood, given the data \code{y},
#' the initial time-varying paramater state \code{f1}, and model parameters
#' \code{params}.
#'
#' \code{$Filter(y, f1)}
#' Computes the set of time-varying parameters for each observation, given the
#' data \code{y}, the initial time-varying paramater state \code{f1}. Returns
#' array like of dimension T by d, where T is the number of observations and d
#' is the number time-varying parameters of the model.
#'
#' @export
NULL

#' @name PriorStack
#'
#' @aliases
#' Rcpp_PriorStack-class Rcpp_PriorStack
#'
#' @title An object that fully specifies the prior for a model.
#'
#' @description
#' Allows for flexible specification of the prior over all model parameters.
#' Currently supports Normal, TruncatedNormal & ImproperUniform.
#'
#' @section Usage:
#' \preformatted{
#' p <- new(PriorStack, priorStrs, priorParams)
#' p <- new(PriorStack, priorStrs, priorParams, priorToParamIndex)
#'
#' p$LogPriors(params)
#' p$LogPriorsWPar(params, priorParams)
#' p$GradLogPriors(params)
#'
#' p$Filter(y, f1)
#' }
#'
#' @section Details:
#'
#' @export
NULL

#' @name GGS
#'
#' @title Griddy Gibbs Sampler.
#'
#' @description
#' Function that returns a sample of Monte Carlo draws given a GAS model
#' of type \code{modelStr} and with a prior specified by the \code{priorStack}.
#' Draws are generated using the Griddy Gibbs Sampler as introduced in
#' Ritter, C., & Tanner, M. A. (1992). "Facilitating the Gibbs sampler:
#' the Gibbs stopper and the griddy-Gibbs sampler".
#'
#' @export
NULL

#' @name RWMH
#'
#' @title Random Walk Metropolis Hastings.
#'
#' @description
#' Function of similar form as \code{GGS}, but draws are generated using a
#' Random Walk Metropolis Hasting's algorithm with multivariate Student's-t
#' proposal density.
#'
#' @export
NULL

#' @name HMC
#'
#' @title Hamoltonian Monte Carlo.
#'
#' @description
#' Function of similar form as \code{GGS}, but draws are generated
#' using the Hamiltonian Monte Carlo algorithm. Usage requires that the gradient
#' for the specified GAS model is implemented (this is currently not the case
#' for the \code{DPMP} class of models).
#'
#' @export
NULL

#' @name VectorizedPosterior
#'
#' @title Vectorized Posterior Distribution Function.
#'
#' @description
#' Vectorized version of the posterior for a GAS model
#' of type \code{modelStr} and with a prior specified by the \code{priorStack}.
#' Used for inference in combination with the \code{\link[AdMit]{AdMit}}
#' sampling algorithms.
#'
#' @export
NULL

## ensure modules gets loaded
loadModule("GASModel", TRUE)
loadModule("Prior", TRUE)
loadModule("RWMH", TRUE)
loadModule("GGS", TRUE)
loadModule("HMC", TRUE)
loadModule("VectorizedPosterior", TRUE)

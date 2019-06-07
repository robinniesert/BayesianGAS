#' @title BayesianGAS
#'
#' @description BayesianGAS: Package contains the source code and (most of) the
#' data necessary to reproduce the analysis presented in the thesis "Bayesian
#' Inference for Generalized Autoregressive Score Models" (2017).
#'
#' @section GASModel functions:
#' The package uses a general framework for GAS models, most importantly
#' providing implementations of the following functions:
#' \itemize{
#'  \item{\code{LogL}} {computes the model's log likelihood}
#'  \item{\code{LogPosterior}} {computes the log of the posterior
#'  for a model (log likelihood + log prior)}
#'  \item{\code{GradLogL}} {gradient of the log likelihood}
#'  \item{\code{Filter}} {computes a model's time-varying parameters.}
#' }
#' See \code{\link{GASModel}} for usage and further description.
#'
#' @section Supported GAS models:
#' The following models are currently supported
#' \itemize{
#'  \item{\code{\link{BetaGenTEGARCH}}}
#'  \item{\code{\link{BetaTEGARCH}}}
#'  \item{\code{\link{DPMP}}} {support for the 1 through to 3 factor models,
#'  and for three different scaling factors: Identity, Inverse Fisher
#'  Information (FI) and Inverse square-root FI}
#' }
#'
#' @section Inference functions:
#' \itemize{
#'  \item{\code{\link{FitML}}} {Maximum Likelihood estimation}
#'  \item{\code{\link{RWMH}}} {Random Walk Metropolis Hastings Sampler}
#'  \item{\code{\link{GGS}}} {Griddy Gibbs Sampler}
#'  \item{\code{\link{HMC}}} {Hamiltonian Monte Carlo}
#' }
#'
#' @section Vignettes:
#' The analysis is presented in three seperate vignettes, corresponding to the
#' three sections of chapter 4 from the thesis.
#'
#' @docType package
#' @name BayesianGAS
NULL

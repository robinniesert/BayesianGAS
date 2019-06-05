#' BayesianGAS: Package contains the source code and (most of) the data
#' necessairy to reproduce the analysis presented in my master thesis Bayesian
#' Inference for Generalized Autoregressive Score (GAS) Models  (2017).
#'
#' @section GASModel functions:
#' The package uses a general framework for GAS models, most importantly
#' providing implementations of the following functions:
#' \itemize{
#'  \item{\code{LogL}}{Computes the model's Log likelihood}
#'  \item{\code{LogPosterior}}{Computes the log of the posterior
#'  for a model (log likelihood + log prior)}
#'  \item{\code{GradLogL}}{Gradient of the log likelihood}
#'  \item{\code{Filter}}{Computes a model's time-varying parameters.}
#' }
#'
#' @section Supported GAS models:
#' The following models are currently supported
#' \itemize{
#'  \item{\code{BetaGenTEGARCH}}
#'  \item{\code{BetaTEGARCH}}
#'  \item{\code{DPMP}}{Support for the 1 through to 3 factor models,
#'  and for three different scaling factors: Identity and Inverse Fisher
#'  Information (FI) and Inverse square-root FI}
#' }
#'
#' @section Inference functions:
#' \itemize{
#'  \item{\code{FitML}}{Maximum Likelihood estimation}
#'  \item{\code{RWMH}}{Random Walk Metropolis Hastings Sampler}
#'  \item{\code{GGS}}{Griddy Gibbs Sampler}
#'  \item{\code{HMC}}{Hamiltonian Monte Carlo}
#' }
#'
#' @section Vignettes:
#' The analysis is presented in three seperate vignettes, corresponding to the
#' three sections of chapter 4 from the thesis.
#'
#' @docType package
#' @name BayesianGAS
NULL
#' Fit GAS Models by Maximum Likilihood (ML).
#'
#' Wrapper function for \code{\link[stats]{optim}} to facilitate ML estimation
#' of GAS model parameters. Overrides some defaults of the
#' \code{\link[stats]{optim}} function to make it suited for Log Likelihood
#' maximization.
#'
#' @param model String or \code{\link{GASModel}} object.
#' @param initParams Vector or list of initial parameters.
#' @param y Array like data object passed to the Log Likelihood function.
#' @param f1 Starting values for the GAS models time varying parameters.
#' @param method,control,hessian,verbose See \code{\link[stats]{optim}} for
#' details.
#' @param ... Other arguments passed to \code{\link[stats]{optim}}.
#'
#' @export
FitML <- function(model, initParams, y, f1, method = 'BFGS',
                  control = list(maxit = 1e5), hessian = TRUE, verbose = TRUE,
                  ...){
  model <- CreateModel(model)
  control <- c(control, fnscale = -1)

  startTime <- Sys.time()
  optimModel   <- stats::optim(
    initParams,
    model$LogLWPar,
    y = y,
    f1 = f1,
    method = method,
    control = control,
    hessian = hessian,
    ...
  )

  model$LogLValML <- optimModel$value
  model$ParamsML <- optimModel$par
  if (hessian) {
    model$StdsML <- sqrt(diag(solve(-optimModel$hessian)))
  }

  if (verbose) {
    cat("ML Log-Likelihood: ", model$LogLValML, sprintf("\n"))
    cat("ML parameter estimates: ", model$ParamsML, sprintf("\n"))
    if (hessian) {
      cat("ML standard errors: ", model$StdsML, sprintf("\n"))
    }
  }

  if (optimModel$convergence > 0) {
    print(optimModel$message)
  }
  return(model)
}

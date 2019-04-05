FitML <- function(model, initParams, y, f1, method = 'BFGS',
                  control = list(maxit = 1e5), hessian = TRUE, verbose = TRUE){
  model <- CreateModel(model)
  control <- c(control, fnscale = -1)

  startTime <- Sys.time()
  optimModel   <- optim(
    initParams,
    model$LogLWPar,
    y = y,
    f1 = f1,
    method = method,
    control = control,
    hessian = hessian
  )

  model$LogLValML <- optimModel$value
  model$ParamsML <- optimModel$par
  model$StdsML <- sqrt(diag(solve(-optimModel$hessian)))

  if (verbose) {
    cat("ML Log-Likelihood: ", model$LogLValML, sprintf("\n"))
    cat("ML parameter estimates: ", model$ParamsML, sprintf("\n"))
    cat("ML standard errors: ", model$ParamsML, sprintf("\n"))
  }

  if (optimModel$convergence > 0) {
    print(optimModel$message)
  }
  return(model)
}

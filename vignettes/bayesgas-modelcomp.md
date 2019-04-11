Model Comparisons: Dynamic Pooled Marked Point Process Models
=============================================================

This document contains the R code to reproduce the plots and statistical
analysis presented in section 4.2 of &gt; Niesert, R. "Bayesian
Inference for Generalized Autoregressive Score Models." (2017).

    library(BayesianGAS)

    set.seed(100)
    kTransitionTypes <- c("IGtoSIG", "IGtoD", "SIGtoIG", "SIGtoD")
    kScalings <- c(0., -0.5, -1.0)
    kNumFactorSpecs <- c(1, 2, 3)
    kNumParamsVec <- c(9, 10, 12)
    kNumTransitions <- 4

Load data
---------

    data("CreditRatings", package = "BayesianGAS")
    dates <- unique(transitionData$datadate)
    colnames(transitions) <- kTransitionTypes
    colnames(possibleTranstions) <- colnames(transitions)
    y <- cbind(transitions, diffTau, possibleTranstions)

### Plot Transitions

    par(mfcol = c(2, 2))
    for (tt in kTransitionTypes) {
      leg <- paste(tt, "counts")
      plot(dates, transitions[,tt], type = "l", xaxt = "n", xlab = "Years", 
           ylab = "Intensity")
      legend("topleft", legend = leg, bty = "n", lty = 1)
      axis.Date(1, at = seq(min(dates), max(dates), by = "1 years"), format = "%Y")
    }

<img src="bayesgas-modelcomp_files/figure-markdown_strict/unnamed-chunk-1-1.png" alt="Transitions"  />
<p class="caption">
Transitions
</p>

### Set data attributes

    numObs <- dim(y)[1]

Maximum Likelihood (ML) estimation
----------------------------------

    initParamVecs <- list(
      c(A = 0.05, B = 0.95, C = c(0.5, 0.5, -0.5) , w = c(-5, -10, -5, -5)),
      c(A = rep(0.05, 2), B = rep(0.95, 2), C = rep(0.5, 2) , w = c(-5, -10, -5, -5)),
      c(A = rep(0.025, 3), B = rep(0.95, 3), C = rep(0.5, 2) , w = c(-5, -10, -5, -5))
    )
    parScaleVecs <- list(
      c(A = 0.01, B = 0.01, C = rep(0.05, 3) , w = rep(0.2, 4)),
      c(A = rep(0.01, 2), B = rep(0.02, 2), C = rep(0.05, 2) , w = rep(0.2, 4)),
      c(A = rep(0.01, 3), B = rep(0.01, 3), C = rep(0.05, 2) , w = rep(0.2, 4))
    )
    lowerBoundVecs <- list(
      c(A = -Inf, B = -1, C = rep(-Inf, 3) , w = rep(-Inf, 4)),
      c(A = rep(-Inf, 2), B = rep(-1, 2), C = rep(-Inf, 2), w = rep(-Inf, 4)),
      c(A = rep(-Inf, 3), B = rep(-1, 3), C = rep(-Inf, 2), w = rep(-Inf, 4))
    )
    upperBoundVecs <- list(
      c(A = Inf, B = 0.99, C = rep(Inf, 3) , w = rep(Inf, 4)),
      c(A = rep(Inf, 2), B = rep(0.99, 2), C = rep(Inf, 2), w = rep(Inf, 4)),
      c(A = rep(Inf, 3), B = rep(0.99, 3), C = rep(Inf, 2), w = rep(Inf, 4))
    )
    f1s <- list(0, rep(0, 2), rep(0, 3))
    modelsML <- list()
    for (i in 1:3) {
      numF <- kNumFactorSpecs[i]
      initParams <- initParamVecs[[i]]
      parScales <- parScaleVecs[[i]]
      lb <- lowerBoundVecs[[i]]
      ub <- upperBoundVecs[[i]]
      f1 <- f1s[[i]]
      for (s in kScalings) {
        dpmp <- new(DPMP, numF, s)
        cat("Fitting model: ", dpmp$Name, sprintf("\n"))
        dpmp <- FitML(
          model = dpmp,
          initParams = initParams,
          y = y,
          f1 = f1,
          method = "L-BFGS-B",
          control = list(
            maxit = 1e5, 
            parscale = parScales
          ),
          hessian = TRUE,
          verbose = TRUE,
          lower = lb,
          upper = ub
        )
        modelsML <- c(modelsML, dpmp)
        names(modelsML) <- c(names(modelsML)[1:length(modelsML) - 1], dpmp$Name)
        cat(sprintf("\n"))
      }
    }
    #> Fitting model:  DPMP1-I 
    #> ML Log-Likelihood:  -16045.04 
    #> ML parameter estimates:  0.1110525 0.977751 0.3530145 0.7057907 0.02003932 -5.514057 -9.993641 -5.943646 -6.059941 
    #> ML standard errors:  0.1110525 0.977751 0.3530145 0.7057907 0.02003932 -5.514057 -9.993641 -5.943646 -6.059941 
    #> 
    #> Fitting model:  DPMP1-H 
    #> ML Log-Likelihood:  -16040.18 
    #> ML parameter estimates:  0.0764421 0.9707993 0.3359581 0.7258449 0.006360033 -5.526309 -10.04185 -5.946748 -6.118867 
    #> ML standard errors:  0.0764421 0.9707993 0.3359581 0.7258449 0.006360033 -5.526309 -10.04185 -5.946748 -6.118867 
    #> 
    #> Fitting model:  DPMP1-Inv 
    #> ML Log-Likelihood:  -16038.11 
    #> ML parameter estimates:  0.04742759 0.9651559 0.3341132 0.7240282 0.0006072111 -5.541887 -10.08178 -5.948991 -6.169874 
    #> ML standard errors:  0.04742759 0.9651559 0.3341132 0.7240282 0.0006072111 -5.541887 -10.08178 -5.948991 -6.169874 
    #> 
    #> Fitting model:  DPMP2-I
    #> Warning in sqrt(diag(solve(-optimModel$hessian))): NaNs produced
    #> ML Log-Likelihood:  -16033.76 
    #> ML parameter estimates:  0.1111847 0.0178569 0.9772338 0.99 0.3543198 0.7069058 -5.513105 -9.989364 -6.033986 -6.056889 
    #> ML standard errors:  0.1111847 0.0178569 0.9772338 0.99 0.3543198 0.7069058 -5.513105 -9.989364 -6.033986 -6.056889 
    #> 
    #> Fitting model:  DPMP2-H 
    #> ML Log-Likelihood:  -16028.39 
    #> ML parameter estimates:  0.07634509 0.01243045 0.9721026 0.9799281 0.3360157 0.7268368 -5.53718 -10.06788 -5.940947 -6.150467 
    #> ML standard errors:  0.07634509 0.01243045 0.9721026 0.9799281 0.3360157 0.7268368 -5.53718 -10.06788 -5.940947 -6.150467 
    #> 
    #> Fitting model:  DPMP2-Inv 
    #> ML Log-Likelihood:  -16026.38 
    #> ML parameter estimates:  0.04677349 0.007397072 0.966722 0.9770044 0.3340264 0.7277383 -5.562198 -10.1275 -5.916964 -6.230971 
    #> ML standard errors:  0.04677349 0.007397072 0.966722 0.9770044 0.3340264 0.7277383 -5.562198 -10.1275 -5.916964 -6.230971 
    #> 
    #> Fitting model:  DPMP3-I
    #> Warning in sqrt(diag(solve(-optimModel$hessian))): NaNs produced
    #> ML Log-Likelihood:  -16015.32 
    #> ML parameter estimates:  0.05100254 0.01789107 0.1258818 0.9720402 0.99 0.9765571 0.4453739 0.6319207 -5.573313 -10.01201 -6.034603 -6.006915 
    #> ML standard errors:  0.05100254 0.01789107 0.1258818 0.9720402 0.99 0.9765571 0.4453739 0.6319207 -5.573313 -10.01201 -6.034603 -6.006915 
    #> 
    #> Fitting model:  DPMP3-H 
    #> ML Log-Likelihood:  -16009.47 
    #> ML parameter estimates:  0.03191623 0.01213572 0.08295099 0.9652229 0.9806023 0.9729215 0.3064406 0.6577841 -5.614742 -10.09733 -5.943177 -6.13637 
    #> ML standard errors:  0.03191623 0.01213572 0.08295099 0.9652229 0.9806023 0.9729215 0.3064406 0.6577841 -5.614742 -10.09733 -5.943177 -6.13637 
    #> 
    #> Fitting model:  DPMP3-Inv 
    #> ML Log-Likelihood:  -16010.36 
    #> ML parameter estimates:  0.01841421 0.006768604 0.03508431 0.9683615 0.9799687 0.9716508 0.3059018 0.5991695 -5.637055 -10.12841 -5.933319 -6.217162 
    #> ML standard errors:  0.01841421 0.006768604 0.03508431 0.9683615 0.9799687 0.9716508 0.3059018 0.5991695 -5.637055 -10.12841 -5.933319 -6.217162

MCMC using RWMH
---------------

I deviate slightly here from the analysis presented in the thesis, by
thinning the posterior sample by a factor of 10 (i.e. I keep only 1 out
of 10 draws). This is done to reduce memory usage.

    iter <- 4e5
    thinning <- 10
    numDraws <- floor(iter / thinning)
    warmUpRounds <- c(3, 5, 6)
    priorStacks <- list(
      new(
        PriorStack,
        c("Normal", "TruncatedNormal", rep("Normal", 7)),
        list(
          c(0.05, 1),
          c(0.95, 1, -1, 1),
          c(0.5, 5),
          c(0.5, 5),
          c(-0.5, 5),
          c(-5, 5),
          c(-10, 5),
          c(-5, 5),
          c(-5, 5)
        )
      ),
      new(
        PriorStack,
        c(rep("Normal", 2), rep("TruncatedNormal", 2), rep("Normal", 6)),
        list(
          c(0.05, 1),
          c(0.05, 1),
          c(0.95, 1, -1, 1),
          c(0.95, 1, -1, 1),
          c(0.5, 5),
          c(0.5, 5),
          c(-5, 5),
          c(-10, 5),
          c(-5, 5),
          c(-5, 5)
        )
      ),
      new(
        PriorStack,
        c(rep("Normal", 3), rep("TruncatedNormal", 3), rep("Normal", 6)),
        list(
          c(0.05, 1),
          c(0.05, 1),
          c(0.05, 1),
          c(0.95, 1, -1, 1),
          c(0.95, 1, -1, 1),
          c(0.95, 1, -1, 1),
          c(0.5, 5),
          c(0.5, 5),
          c(-5, 5),
          c(-10, 5),
          c(-5, 5),
          c(-5, 5)
        )
      )
    )

### Run RWMH

    drawsRWMHLst <- list()
    for (i in 1:3) {
      numF <- kNumFactorSpecs[i]
      initParams <- initParamVecs[[i]]
      f1 <- f1s[[i]]
      priorStack <- priorStacks[[i]]
      numParams <- kNumParamsVec[i]
      for (s in kScalings) {
        if ((s  == -1) && (numF > 1)) {
          stepsize1 <- 0.0025
        }else{
          stepsize1 <- 0.006
        }
        dpmp <- new(DPMP, numF, s)
        cat("Running RWMH for model: ", dpmp$Name, sprintf("\n"))
        startTime <- Sys.time()
        cat(sprintf("Warm up 1 \n"))
        warmUpRWMH <- RWMH(
          dpmp$Name,
          priorStack,
          y = y,
          f1 = f1,
          initParams = initParams,
          sigma = diag(numParams),
          iter = 1e4,
          stepsize = stepsize1,
          printIter = 1e5,
          thinning = thinning
        )
        for (round in 2:max(2, warmUpRounds[i])) {
          cat(sprintf("Warm up %i \n", round))
          warmUpRWMH <- RWMH(
            dpmp$Name,
            priorStack,
            y = y,
            f1 = f1,
            initParams = initParams,
            sigma = cov(warmUpRWMH),
            iter = 1e4,
            stepsize = .4,
            printIter = 1e5,
            thinning = thinning
          )
        }
        drawsRWMH <- RWMH(
          dpmp$Name,
          priorStack,
          y = y,
          f1 = f1,
          initParams = initParams,
          sigma = cov(warmUpRWMH),
          iter = iter,
          stepsize = .4,
          printIter = 1e5,
          thinning = thinning
        )
        endTime <- Sys.time()
        timeRWMH <- difftime(endTime, startTime, units = 'secs')
        cat("RWMH Time: ", timeRWMH, sprintf(" seconds\n"))
        colnames(drawsRWMH) <- names(initParams)
        drawsRWMHLst <- append(drawsRWMHLst, list(drawsRWMH))
        names(drawsRWMHLst) <- 
          c(names(drawsRWMHLst)[1:length(drawsRWMHLst) - 1], dpmp$Name)
        
        ESSs <- coda::effectiveSize(drawsRWMH)
        ESSs <- t(data.frame(ESSs))
        colnames(ESSs) <- names(initParams)
        print(round(ESSs, 1))
        cat(sprintf("\n"))
      }
    }

    #> Running RWMH for model:  DPMP1-I 
    #> Warm up 1 
    #> RWMH - Accept ratio is: 0.465 
    #> Warm up 2 
    #> RWMH - Accept ratio is: 0.355 
    #> Warm up 3 
    #> RWMH - Accept ratio is: 0.414 
    #> iter 100000
    #> iter 200000
    #> iter 300000
    #> RWMH - Accept ratio is: 0.397 
    #> RWMH Time:  49.17175  seconds
    #>           A      B     C1   C2     C3     w1     w2     w3     w4
    #> ESSs 7403.7 5944.1 7480.7 6768 7239.4 5974.9 6649.3 9003.2 5989.5
    #> 
    #> Running RWMH for model:  DPMP1-H 
    #> Warm up 1 
    #> RWMH - Accept ratio is: 0.513 
    #> Warm up 2 
    #> RWMH - Accept ratio is: 0.302 
    #> Warm up 3 
    #> RWMH - Accept ratio is: 0.406 
    #> iter 100000
    #> iter 200000
    #> iter 300000
    #> RWMH - Accept ratio is: 0.372 
    #> RWMH Time:  76.20131  seconds
    #>           A    B   C1     C2     C3     w1     w2     w3     w4
    #> ESSs 6574.1 5999 6557 5730.4 7514.2 6544.4 5482.5 8358.8 6570.2
    #> 
    #> Running RWMH for model:  DPMP1-Inv 
    #> Warm up 1 
    #> RWMH - Accept ratio is: 0.467 
    #> Warm up 2 
    #> RWMH - Accept ratio is: 0.271 
    #> Warm up 3 
    #> RWMH - Accept ratio is: 0.431 
    #> iter 100000
    #> iter 200000
    #> iter 300000
    #> RWMH - Accept ratio is: 0.364 
    #> RWMH Time:  75.74495  seconds
    #>           A      B     C1   C2     C3     w1     w2     w3     w4
    #> ESSs 6046.6 5880.1 6538.7 5571 7168.4 4941.3 4433.3 8422.1 4837.9
    #> 
    #> Running RWMH for model:  DPMP2-I 
    #> Warm up 1 
    #> RWMH - Accept ratio is: 0.289 
    #> Warm up 2 
    #> RWMH - Accept ratio is: 0.288 
    #> Warm up 3 
    #> RWMH - Accept ratio is: 0.366 
    #> Warm up 4 
    #> RWMH - Accept ratio is: 0.321 
    #> Warm up 5 
    #> RWMH - Accept ratio is: 0.328 
    #> iter 100000
    #> iter 200000
    #> iter 300000
    #> RWMH - Accept ratio is: 0.318 
    #> RWMH Time:  44.13673  seconds
    #>          A1     A2   B1    B2     C1     C2     w1     w2     w3     w4
    #> ESSs 5512.5 2102.7 4603 391.1 5251.4 5239.7 4922.4 4700.8 2760.6 4795.2
    #> 
    #> Running RWMH for model:  DPMP2-H 
    #> Warm up 1 
    #> RWMH - Accept ratio is: 0.285 
    #> Warm up 2 
    #> RWMH - Accept ratio is: 0.289 
    #> Warm up 3 
    #> RWMH - Accept ratio is: 0.327 
    #> Warm up 4 
    #> RWMH - Accept ratio is: 0.315 
    #> Warm up 5 
    #> RWMH - Accept ratio is: 0.332 
    #> iter 100000
    #> iter 200000
    #> iter 300000
    #> RWMH - Accept ratio is: 0.339 
    #> RWMH Time:  96.99639  seconds
    #>        A1     A2     B1    B2     C1     C2     w1     w2     w3     w4
    #> ESSs 6155 2964.9 6097.4 210.4 5074.8 6098.9 4963.3 5657.5 1238.2 4735.3
    #> 
    #> Running RWMH for model:  DPMP2-Inv 
    #> Warm up 1 
    #> RWMH - Accept ratio is: 0.444 
    #> Warm up 2 
    #> RWMH - Accept ratio is: 0.429 
    #> Warm up 3 
    #> RWMH - Accept ratio is: 0.387 
    #> Warm up 4 
    #> RWMH - Accept ratio is: 0.352 
    #> Warm up 5 
    #> RWMH - Accept ratio is: 0.355 
    #> iter 100000
    #> iter 200000
    #> iter 300000
    #> RWMH - Accept ratio is: 0.320 
    #> RWMH Time:  77.70865  seconds
    #>          A1     A2     B1   B2     C1     C2     w1   w2     w3     w4
    #> ESSs 5256.7 5807.4 4824.5 2565 5218.7 4942.6 4990.6 4252 1197.3 5096.5
    #> 
    #> Running RWMH for model:  DPMP3-I 
    #> Warm up 1 
    #> RWMH - Accept ratio is: 0.247 
    #> Warm up 2 
    #> RWMH - Accept ratio is: 0.279 
    #> Warm up 3 
    #> RWMH - Accept ratio is: 0.298 
    #> Warm up 4 
    #> RWMH - Accept ratio is: 0.301 
    #> Warm up 5 
    #> RWMH - Accept ratio is: 0.304 
    #> Warm up 6 
    #> RWMH - Accept ratio is: 0.245 
    #> iter 100000
    #> iter 200000
    #> iter 300000
    #> RWMH - Accept ratio is: 0.271 
    #> RWMH Time:  46.36389  seconds
    #>        A1     A2   A3     B1    B2     B3     C1     C2     w1     w2
    #> ESSs 3911 2571.3 5959 3210.5 504.3 2833.4 3829.4 4008.4 2856.1 4044.2
    #>          w3     w4
    #> ESSs 3424.5 3899.8
    #> 
    #> Running RWMH for model:  DPMP3-H 
    #> Warm up 1 
    #> RWMH - Accept ratio is: 0.199 
    #> Warm up 2 
    #> RWMH - Accept ratio is: 0.406 
    #> Warm up 3 
    #> RWMH - Accept ratio is: 0.336 
    #> Warm up 4 
    #> RWMH - Accept ratio is: 0.305 
    #> Warm up 5 
    #> RWMH - Accept ratio is: 0.308 
    #> Warm up 6 
    #> RWMH - Accept ratio is: 0.300 
    #> iter 100000
    #> iter 200000
    #> iter 300000
    #> RWMH - Accept ratio is: 0.294 
    #> RWMH Time:  106.914  seconds
    #>          A1     A2     A3     B1  B2     B3     C1     C2     w1     w2
    #> ESSs 3239.2 3136.7 5326.4 3270.4 441 4055.9 3948.3 4381.7 3063.1 3921.6
    #>         w3     w4
    #> ESSs 930.8 4652.8
    #> 
    #> Running RWMH for model:  DPMP3-Inv 
    #> Warm up 1 
    #> RWMH - Accept ratio is: 0.205 
    #> Warm up 2 
    #> RWMH - Accept ratio is: 0.609 
    #> Warm up 3 
    #> RWMH - Accept ratio is: 0.317 
    #> Warm up 4 
    #> RWMH - Accept ratio is: 0.296 
    #> Warm up 5 
    #> RWMH - Accept ratio is: 0.321 
    #> Warm up 6 
    #> RWMH - Accept ratio is: 0.297 
    #> iter 100000
    #> iter 200000
    #> iter 300000
    #> RWMH - Accept ratio is: 0.305 
    #> RWMH Time:  86.54065  seconds
    #>          A1     A2     A3     B1    B2     B3     C1     C2     w1     w2
    #> ESSs 4452.1 3493.4 3518.1 5409.5 554.6 3928.9 4721.8 4111.7 3225.9 4730.3
    #>         w3   w4
    #> ESSs 748.9 3620

Model comparisons
-----------------

    IC <- function(logl, npar, k = log(npar)){
      IC <- -2 * logl + k * npar
      return(IC)
    }
    modelScores <- data.frame(matrix(
      0, 
      nrow = 9, 
      ncol = 3, 
      dimnames = list(names(modelsML), c("MarginalLikelihood", "LogLikelihood", "BIC"))
    ))
    burn <- 1000

### Compute marginals and Bayesian Information criteria (BICs)

    upperBoundVecs <- list(
      c(A = Inf, B = 1., C = rep(Inf, 3) , w = rep(Inf, 4)),
      c(A = rep(Inf, 2), B = rep(1., 2), C = rep(Inf, 2), w = rep(Inf, 4)),
      c(A = rep(Inf, 3), B = rep(1., 3), C = rep(Inf, 2), w = rep(Inf, 4))
    )
    marginalsLst <- list()
    logLikList <- list()
    for (i in 1:3) {
      numF <- kNumFactorSpecs[i]
      initParams <- initParamVecs[[i]]
      f1 <- f1s[[i]]
      priorStack <- priorStacks[[i]]
      numParams <- kNumParamsVec[i]
      lb <- lowerBoundVecs[[i]]
      ub <- upperBoundVecs[[i]]
      for (s in kScalings) {
        dpmp <- new(DPMP, initParams, priorStack, numF, s)
        logPosterior <- function(pars, data, printErrors = FALSE) {
          out <- tryCatch(
            {dpmp$LogPosteriorLWPar(pars, y, f1)},
            error = function(cond) {
              if (printErrors) message(cond)
              return(-Inf)
            }
          )
          return(out)
        }
        marginal <- bridgesampling::bridge_sampler(
          drawsRWMHLst[[dpmp$Name]][-(1:burn), ], 
          log_posterior = logPosterior,
          data = NULL,
          method = "warp3",
          lb = lb,
          ub = ub
        )
        modelScores[dpmp$Name, "MarginalLikelihood"] <- marginal$logml
        modelScores[dpmp$Name, "LogLikelihood"] <- modelsML[[dpmp$Name]]$LogLValML
        modelScores[dpmp$Name, "BIC"] <- IC(modelsML[[dpmp$Name]]$LogLValML, numParams)
      }
    }
    #> Warning: 19500 of the 19500 log_prob() evaluations on the warp-transformed
    #> posterior draws produced -Inf/Inf.
    #> Warning: 19500 of the 19500 log_prob() evaluations on the warp-transformed
    #> proposal draws produced -Inf/Inf.
    #> Error in out[!from@positive] <- -out[!from@positive]: NAs are not allowed in subscripted assignments

<table>
<caption>Marginals, log-likelihoods and BICS</caption>
<thead>
<tr class="header">
<th></th>
<th align="right">MarginalLikelihood</th>
<th align="right">LogLikelihood</th>
<th align="right">BIC</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>DPMP1-I</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td>DPMP1-H</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td>DPMP1-Inv</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td>DPMP2-I</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td>DPMP2-H</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td>DPMP2-Inv</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td>DPMP3-I</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td>DPMP3-H</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td>DPMP3-Inv</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
</tbody>
</table>

### Compute Bayes Factors (BFs)

    marginals <- modelScores["MarginalLikelihood"]
    bf1H1I <- bridgesampling::bayes_factor(
      marginals["DPMP1-H", ], marginals["DPMP1-I", ], TRUE)
    bf1Inv1H <- bridgesampling::bayes_factor(
      marginals["DPMP1-Inv", ], marginals["DPMP1-H", ], TRUE)
    bf2H2I <- bridgesampling::bayes_factor(
      marginals["DPMP2-H", ], marginals["DPMP2-I", ], TRUE)
    bf2Inv2H <- bridgesampling::bayes_factor(
      marginals["DPMP2-Inv", ], marginals["DPMP2-H", ], TRUE)
    bf3H3I <- bridgesampling::bayes_factor(
      marginals["DPMP3-H", ], marginals["DPMP3-I", ], TRUE)
    bf3Inv3H <- bridgesampling::bayes_factor(
      marginals["DPMP3-Inv", ], marginals["DPMP3-H", ], TRUE)
    bf2I1I <- bridgesampling::bayes_factor(
      marginals["DPMP2-I", ], marginals["DPMP1-I", ], TRUE)
    bf2H1H <- bridgesampling::bayes_factor(
      marginals["DPMP2-H", ], marginals["DPMP1-H", ], TRUE)
    bf2Inv1Inv <- bridgesampling::bayes_factor(
      marginals["DPMP2-Inv", ], marginals["DPMP1-Inv", ], TRUE)
    bf3I1I <- bridgesampling::bayes_factor(
      marginals["DPMP3-I", ], marginals["DPMP1-I", ], TRUE)
    bf3H1H <- bridgesampling::bayes_factor(
      marginals["DPMP3-H", ], marginals["DPMP1-H", ], TRUE)
    bf3Inv1Inv <- bridgesampling::bayes_factor(
      marginals["DPMP3-Inv", ], marginals["DPMP1-Inv", ], TRUE)

1-H | 1-I : 0  
1-Inv | 1-H : 0  
2-H | 2-I : 0  
2-Inv | 2-H : 0  
3-H | 3-I : 0  
3-Inv | 3-H : 0  
2-I | 1-I : 0  
2-H | 1-H : 0  
2-Inv | 1-Inv : 0  
3-I | 1-I : 0  
3-H | 1-H : 0  
3-Inv | 1-Inv : 0

Some Parameter Statistics
-------------------------

    selectedModels <- c("DPMP1-Inv", "DPMP2-Inv", "DPMP3-Inv")
    for (model in selectedModels) {
      summary_ <- summary(coda::mcmc(drawsRWMHLst[[model]]))$statistics
      print(knitr::kable(summary_, caption = model))
    }

<table>
<caption>DPMP1-Inv</caption>
<thead>
<tr class="header">
<th></th>
<th align="right">Mean</th>
<th align="right">SD</th>
<th align="right">Naive SE</th>
<th align="right">Time-series SE</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>A</td>
<td align="right">0.0482675</td>
<td align="right">0.0045976</td>
<td align="right">0.0000230</td>
<td align="right">0.0000591</td>
</tr>
<tr class="even">
<td>B</td>
<td align="right">0.9667727</td>
<td align="right">0.0147314</td>
<td align="right">0.0000737</td>
<td align="right">0.0001921</td>
</tr>
<tr class="odd">
<td>C1</td>
<td align="right">0.3352045</td>
<td align="right">0.0361157</td>
<td align="right">0.0001806</td>
<td align="right">0.0004466</td>
</tr>
<tr class="even">
<td>C2</td>
<td align="right">0.7369001</td>
<td align="right">0.3041936</td>
<td align="right">0.0015210</td>
<td align="right">0.0040755</td>
</tr>
<tr class="odd">
<td>C3</td>
<td align="right">-0.0007197</td>
<td align="right">0.0359887</td>
<td align="right">0.0001799</td>
<td align="right">0.0004251</td>
</tr>
<tr class="even">
<td>w1</td>
<td align="right">-5.5503212</td>
<td align="right">0.1156935</td>
<td align="right">0.0005785</td>
<td align="right">0.0016458</td>
</tr>
<tr class="odd">
<td>w2</td>
<td align="right">-10.1940548</td>
<td align="right">0.4553736</td>
<td align="right">0.0022769</td>
<td align="right">0.0068392</td>
</tr>
<tr class="even">
<td>w3</td>
<td align="right">-5.9510990</td>
<td align="right">0.0401341</td>
<td align="right">0.0002007</td>
<td align="right">0.0004373</td>
</tr>
<tr class="odd">
<td>w4</td>
<td align="right">-6.1885865</td>
<td align="right">0.3325450</td>
<td align="right">0.0016627</td>
<td align="right">0.0047810</td>
</tr>
</tbody>
</table>

<table>
<caption>DPMP2-Inv</caption>
<thead>
<tr class="header">
<th></th>
<th align="right">Mean</th>
<th align="right">SD</th>
<th align="right">Naive SE</th>
<th align="right">Time-series SE</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>A1</td>
<td align="right">0.0475394</td>
<td align="right">0.0046382</td>
<td align="right">0.0000232</td>
<td align="right">0.0000640</td>
</tr>
<tr class="even">
<td>A2</td>
<td align="right">0.0089956</td>
<td align="right">0.0028407</td>
<td align="right">0.0000142</td>
<td align="right">0.0000373</td>
</tr>
<tr class="odd">
<td>B1</td>
<td align="right">0.9679876</td>
<td align="right">0.0145260</td>
<td align="right">0.0000726</td>
<td align="right">0.0002091</td>
</tr>
<tr class="even">
<td>B2</td>
<td align="right">0.9643678</td>
<td align="right">0.0331507</td>
<td align="right">0.0001658</td>
<td align="right">0.0006546</td>
</tr>
<tr class="odd">
<td>C1</td>
<td align="right">0.3359623</td>
<td align="right">0.0361130</td>
<td align="right">0.0001806</td>
<td align="right">0.0004999</td>
</tr>
<tr class="even">
<td>C2</td>
<td align="right">0.7429017</td>
<td align="right">0.3037809</td>
<td align="right">0.0015189</td>
<td align="right">0.0043210</td>
</tr>
<tr class="odd">
<td>w1</td>
<td align="right">-5.5686462</td>
<td align="right">0.1172325</td>
<td align="right">0.0005862</td>
<td align="right">0.0016595</td>
</tr>
<tr class="even">
<td>w2</td>
<td align="right">-10.2254297</td>
<td align="right">0.4519668</td>
<td align="right">0.0022598</td>
<td align="right">0.0069312</td>
</tr>
<tr class="odd">
<td>w3</td>
<td align="right">-5.9537674</td>
<td align="right">0.1363085</td>
<td align="right">0.0006815</td>
<td align="right">0.0039393</td>
</tr>
<tr class="even">
<td>w4</td>
<td align="right">-6.2411430</td>
<td align="right">0.3369760</td>
<td align="right">0.0016849</td>
<td align="right">0.0047202</td>
</tr>
</tbody>
</table>

<table>
<caption>DPMP3-Inv</caption>
<thead>
<tr class="header">
<th></th>
<th align="right">Mean</th>
<th align="right">SD</th>
<th align="right">Naive SE</th>
<th align="right">Time-series SE</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>A1</td>
<td align="right">0.0191196</td>
<td align="right">0.0030494</td>
<td align="right">0.0000152</td>
<td align="right">0.0000457</td>
</tr>
<tr class="even">
<td>A2</td>
<td align="right">0.0082260</td>
<td align="right">0.0026635</td>
<td align="right">0.0000133</td>
<td align="right">0.0000451</td>
</tr>
<tr class="odd">
<td>A3</td>
<td align="right">0.0352054</td>
<td align="right">0.0029403</td>
<td align="right">0.0000147</td>
<td align="right">0.0000496</td>
</tr>
<tr class="even">
<td>B1</td>
<td align="right">0.9675531</td>
<td align="right">0.0166330</td>
<td align="right">0.0000832</td>
<td align="right">0.0002261</td>
</tr>
<tr class="odd">
<td>B2</td>
<td align="right">0.9669220</td>
<td align="right">0.0363751</td>
<td align="right">0.0001819</td>
<td align="right">0.0015446</td>
</tr>
<tr class="even">
<td>B3</td>
<td align="right">0.9701455</td>
<td align="right">0.0130379</td>
<td align="right">0.0000652</td>
<td align="right">0.0002080</td>
</tr>
<tr class="odd">
<td>C1</td>
<td align="right">0.2888558</td>
<td align="right">0.9755994</td>
<td align="right">0.0048780</td>
<td align="right">0.0141976</td>
</tr>
<tr class="even">
<td>C2</td>
<td align="right">0.6290629</td>
<td align="right">0.4331183</td>
<td align="right">0.0021656</td>
<td align="right">0.0067545</td>
</tr>
<tr class="odd">
<td>w1</td>
<td align="right">-5.6385708</td>
<td align="right">0.1697388</td>
<td align="right">0.0008487</td>
<td align="right">0.0029885</td>
</tr>
<tr class="even">
<td>w2</td>
<td align="right">-10.2349726</td>
<td align="right">0.4615622</td>
<td align="right">0.0023078</td>
<td align="right">0.0067110</td>
</tr>
<tr class="odd">
<td>w3</td>
<td align="right">-6.0010384</td>
<td align="right">0.1772896</td>
<td align="right">0.0008864</td>
<td align="right">0.0064786</td>
</tr>
<tr class="even">
<td>w4</td>
<td align="right">-6.2243190</td>
<td align="right">0.3105999</td>
<td align="right">0.0015530</td>
<td align="right">0.0051624</td>
</tr>
</tbody>
</table>

Plots
-----

### Posterior of intensity

    intensityDraws <- 
      array(0, dim = c(numDraws - burn, numObs, kNumTransitions))
    meanLogIntensitiesLst <- list()
    for (i in c(1, 3)) {
      numF <- kNumFactorSpecs[i]
      f1 <- f1s[[i]]
      for (s in kScalings[c(1,3)]) {
        dpmp <- new(DPMP, numF, s)
        for (i in 1:(numDraws - burn)) {
          dpmp$SetParams(as.vector(drawsRWMHLst[[dpmp$Name]][burn + i,]))
          intensityDraws[i, , ] <- dpmp$IntensityFilter(y, f1, TRUE)
          if ((i > 0) && ((i %% 1e4) == 0)) {
            cat(sprintf("iter %i\n", i));
          }
        }
        meanLogIntensities <- colMeans(intensityDraws)
        colnames(meanLogIntensities) <- kTransitionTypes
        meanLogIntensitiesLst <- 
          append(meanLogIntensitiesLst, list(meanLogIntensities))
        names(meanLogIntensitiesLst) <- c(
          names(meanLogIntensitiesLst)[1:length(meanLogIntensitiesLst) - 1], dpmp$Name)
      }
    }
    #> iter 10000
    #> iter 20000
    #> iter 30000
    #> iter 10000
    #> iter 20000
    #> iter 30000
    #> iter 10000
    #> iter 20000
    #> iter 30000
    #> iter 10000
    #> iter 20000
    #> iter 30000
    intensityDraws <- exp(intensityDraws)

### Mean log intensity plots

    selectedModels <- c("DPMP1-I", "DPMP3-I", "DPMP1-Inv", "DPMP3-Inv")
    labY <- c(
      expression("Log intensity" ~ IG %->% SIG), 
      expression("Log intensity" ~ IG %->% D), 
      expression("Log intensity" ~ SIG %->% IG), 
      expression("Log intensity" ~ SIG %->% D)
    )
    names(labY) <- kTransitionTypes
    legMeans <- c(
      "1 factor Identity scaling",
      "3 factor Identity scaling",
      "1 factor Inverse FI scaling",
      "3 factor Inverse FI scaling"
    )
    names(legMeans) <- selectedModels
    yLims <- list(c(-6.5, -4), c(-12.5, -8), c(-6.5, -5.4), c(-9, -2.9))
    names(yLims) <- kTransitionTypes
    par(mfcol = c(2, 2))
    for (tt in kTransitionTypes) {
      for (model in selectedModels) {
        if (grepl("3", model)) lty = 3 else lty = 1
        if (grepl("Inv", model)) col = "blue" else col = "red"
        if (model == selectedModels[1]) {
          plot(dates, meanLogIntensitiesLst[[model]][,tt], type = "l", xaxt = "n", 
               xlab = "Years",  ylab = labY[tt], lty = lty, ylim = yLims[[tt]],
               col = col)
        }else{
          lines(dates, meanLogIntensitiesLst[[model]][,tt], lty = lty, col = col)
        }
        axis.Date(1, at = seq(min(dates), max(dates), by = "1 years"), format = "%Y")
      }
      legend("topleft", bty = "n", col = c("red", "red", "blue", "blue"), 
             lty=c(1, 3, 1, 3), legend = legMeans)
    }

<img src="bayesgas-modelcomp_files/figure-markdown_strict/unnamed-chunk-8-1.png" alt="Mean Log Intensities"  />
<p class="caption">
Mean Log Intensities
</p>

### Joint distribution plots

    selectedDraws <- drawsRWMHLst[["DPMP3-Inv"]][-(1:burn), ]
    par(mfrow = c(2, 1), mar = c(4.2, 4.2, 1, 2))
    plot(
      selectedDraws[, "B2"],
      selectedDraws[, "A2"],
      cex = 0.2,
      cex.axis = 1,
      xlab = expression(b[2]),
      ylab = expression(a[2])
    )
    abline(h = coda::HPDinterval(coda::mcmc(selectedDraws[, "A2"])), col = "red", 
           lty = 2)
    leg <- expression("95% HPD bounds for" ~ a[2])
    legend("topleft", legend = leg, col = "red", lty = 2,  bty = "n")

    ix <- which((selectedDraws[, "C1"] <= 0) & (selectedDraws[, "C2"] <= 0))
    plot(
      selectedDraws[, "C1"],
      selectedDraws[, "C2"],
      cex = 0.2,
      cex.axis = 1,
      xlab = expression(Z[1]),
      ylab = expression(Z[2])
    )
    lines(c(-5, 0), c(0, 0), lty = 2, lwd = 1)
    lines(c(0, 0), c(-2, 0), lty = 2, lwd = 1)
    points(selectedDraws[, "C1"][ix], selectedDraws[, "C2"][ix], col = "red",
           cex = 0.2, pch = 19)

<img src="bayesgas-modelcomp_files/figure-markdown_strict/unnamed-chunk-9-1.png" style="display: block; margin: auto;" />

### Highest Posterior Density (HPD) intensity plots

    par(mfcol = c(2, 2))
    ttIdx <- 0
    for (tt in kTransitionTypes) {
      ttIdx <- ttIdx + 1
      PlotHPDOverTime(
        intensityDraws[, , ttIdx] * max(transitionData$nummonth),
        transitions[, tt],
        dates,
        ylab = "Intensity",
        statStr = "Intensity",
        obsStr = paste(tt, "counts"),
        modeCol = rgb(0,0,1,1), 
        fillCol = rgb(0,0,1,1/4),
        borderCol = rgb(0,0,1,1/2),
        newPlot = FALSE,
        ylim = NULL,
        dateAxisStep = "1 year"
      )
    }

<img src="bayesgas-modelcomp_files/figure-markdown_strict/unnamed-chunk-10-1.png" alt="HPD plots"  />
<p class="caption">
HPD plots
</p>

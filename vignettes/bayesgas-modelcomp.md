---
title: "Model Comparisons: Dynamic Pooled Marked Point Process Models"
author: "Robin Niesert"
date: "April 2nd, 2019"
output: 
  rmarkdown::html_vignette:
    keep_md: yes
vignette: >
  %\VignetteIndexEntry{Model Comparisons: Dynamic Pooled Marked Point Process Models}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
  %\VignetteDepends{}
---



# Model Comparisons: Dynamic Pooled Marked Point Process Models

This document contains the R code to reproduce the plots and statistical analysis presented in section 4.2 of 
> Niesert, R. "Bayesian Inference for Generalized Autoregressive Score Models." (2017).


```r
library(BayesianGAS)
```


```r
set.seed(100)
kTransitionTypes <- c("IGtoSIG", "IGtoD", "SIGtoIG", "SIGtoD")
kScalings <- c(0., -0.5, -1.0)
kNumFactorSpecs <- c(1, 2, 3)
kNumParamsVec <- c(9, 10, 12)
kNumTransitions <- 4
```

## Simulate Data 


```r
simParams <- c(
  A = c(0.05, 0.01),
  B = c(0.96, 0.72),
  C = c(0.39, 0.89),
  W = c(-5.42, -10.07, -5.45, -5.71)
)
dpmpSimModel <- new(DPMP, simParams, 2, -1.)
sims <- dpmpSimModel$Simulate(750, 1000, 4000, rep(0, 2))
tauRaw <- as.integer(cumsum(sims[,5]))
transitionsRaw <- data.frame(sims[,c(1:4)])
possibleTranstionsRaw <- data.frame(sims[,c(6:9)])
colnames(possibleTranstionsRaw) <- kTransitionTypes
colnames(transitionsRaw) <- kTransitionTypes
transitions <- as.matrix(aggregate(transitionsRaw, by = list(tauRaw), FUN = sum))
tau <- transitions[, 'Group.1']
diffTau <- diff(tau)
transitions <- transitions[-1, kTransitionTypes]
endTau <- tail(tau, n = 1)
dates <- seq(as.Date("1986-02-01"), by = "month", length.out = endTau)[tau] - 1
possibleTranstions <- as.matrix(
  aggregate(possibleTranstionsRaw, by = list(tauRaw), 
            FUN = function(x){tail(x, n = 1)}))[-1, kTransitionTypes]
y <- cbind(transitions, diffTau, possibleTranstions)
```

### Plot Transitions

```r
par(mfcol = c(2, 2))
for (tt in kTransitionTypes) {
  leg <- paste(tt, "counts")
  plot(dates, transitions[,tt], type = "l", xaxt = "n", xlab = "Years", 
       ylab = "Intensity")
  legend("topleft", legend = leg, bty = "n", lty = 1)
  axis.Date(1, at = seq(min(dates), max(dates), by = "1 years"), format = "%Y")
}
```

<div class="figure" style="text-align: center">
<img src="bayesgas-modelcomp_files/figure-html/unnamed-chunk-15-1.png" alt="Transitions"  />
<p class="caption">Transitions</p>
</div>
### Set data attributes

```r
numObs <- dim(y)[1]
```

## Maximum Likelihood (ML) estimation 


```r
initParamVecs <- list(
  c(A = 0.05, B = 0.95, C = c(0.5, 0.5, -0.5) , w = c(-5, -10, -5, -5)),
  c(A = rep(0.025, 2), B = rep(0.8, 2), C = rep(0.5, 2) , w = c(-5, -10, -5, -5)),
  c(A = rep(0.025, 3), B = rep(0.8, 3), C = rep(0.5, 2) , w = c(-5, -10, -5, -5))
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
  c(A = rep(Inf, 2), B = rep(0.999, 2), C = rep(Inf, 2), w = rep(Inf, 4)),
  c(A = rep(Inf, 3), B = rep(0.999, 3), C = rep(Inf, 2), w = rep(Inf, 4))
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
#> ML Log-Likelihood:  -25997.16 
#> ML parameter estimates:  0.0931982 0.8037777 0.2057308 0.5010474 0.1168047 -5.422278 -10.0544 -5.47315 -5.738657 
#> ML standard errors:  0.0931982 0.8037777 0.2057308 0.5010474 0.1168047 -5.422278 -10.0544 -5.47315 -5.738657 
#> 
#> Fitting model:  DPMP1-H 
#> ML Log-Likelihood:  -25996.99 
#> ML parameter estimates:  0.05080459 0.8014219 0.2092616 0.4957341 0.114561 -5.422582 -10.05288 -5.473332 -5.740393 
#> ML standard errors:  0.05080459 0.8014219 0.2092616 0.4957341 0.114561 -5.422582 -10.05288 -5.473332 -5.740393 
#> 
#> Fitting model:  DPMP1-Inv 
#> ML Log-Likelihood:  -25997.01 
#> ML parameter estimates:  0.0271067 0.8000182 0.2130024 0.4976545 0.1102174 -5.42284 -10.05446 -5.472889 -5.741045 
#> ML standard errors:  0.0271067 0.8000182 0.2130024 0.4976545 0.1102174 -5.42284 -10.05446 -5.472889 -5.741045 
#> 
#> Fitting model:  DPMP2-I 
#> ML Log-Likelihood:  -25997.38 
#> ML parameter estimates:  0.09376317 0.01158158 0.7962677 0.7968437 0.2081367 0.5046863 -5.422421 -10.05542 -5.470769 -5.739639 
#> ML standard errors:  0.09376317 0.01158158 0.7962677 0.7968437 0.2081367 0.5046863 -5.422421 -10.05542 -5.470769 -5.739639 
#> 
#> Fitting model:  DPMP2-H 
#> ML Log-Likelihood:  -25997.15 
#> ML parameter estimates:  0.05078986 0.0069546 0.7979974 0.7926063 0.208281 0.5008425 -5.422708 -10.05269 -5.471185 -5.741572 
#> ML standard errors:  0.05078986 0.0069546 0.7979974 0.7926063 0.208281 0.5008425 -5.422708 -10.05269 -5.471185 -5.741572 
#> 
#> Fitting model:  DPMP2-Inv 
#> ML Log-Likelihood:  -25997.1 
#> ML parameter estimates:  0.02698514 0.004079196 0.7950276 0.7987161 0.215911 0.5033953 -5.422882 -10.0537 -5.471811 -5.741671 
#> ML standard errors:  0.02698514 0.004079196 0.7950276 0.7987161 0.215911 0.5033953 -5.422882 -10.0537 -5.471811 -5.741671 
#> 
#> Fitting model:  DPMP3-I
#> Warning in sqrt(diag(solve(-optimModel$hessian))): NaNs produced
#> ML Log-Likelihood:  -26000.05 
#> ML parameter estimates:  0.006020387 0.01139102 0.09474762 0.7983493 0.7994143 0.8046904 0.5017987 0.4982762 -5.421453 -10.05366 -5.47083 -5.739326 
#> ML standard errors:  0.006020387 0.01139102 0.09474762 0.7983493 0.7994143 0.8046904 0.5017987 0.4982762 -5.421453 -10.05366 -5.47083 -5.739326 
#> 
#> Fitting model:  DPMP3-H
#> Warning in sqrt(diag(solve(-optimModel$hessian))): NaNs produced
#> ML Log-Likelihood:  -26000.04 
#> ML parameter estimates:  0.003761458 0.007045877 0.04944199 0.7563749 0.7833416 0.807977 0.5707379 0.1571744 -5.421588 -10.04823 -5.471177 -5.740302 
#> ML standard errors:  0.003761458 0.007045877 0.04944199 0.7563749 0.7833416 0.807977 0.5707379 0.1571744 -5.421588 -10.04823 -5.471177 -5.740302 
#> 
#> Fitting model:  DPMP3-Inv 
#> ML Log-Likelihood:  -26000.46 
#> ML parameter estimates:  0.001120895 0.004211854 0.02473002 0.7938715 0.7943959 0.8118626 0.4984743 0.3991225 -5.422084 -10.04447 -5.471421 -5.739309 
#> ML standard errors:  0.001120895 0.004211854 0.02473002 0.7938715 0.7943959 0.8118626 0.4984743 0.3991225 -5.422084 -10.04447 -5.471421 -5.739309
```

## MCMC using RWMH

I deviate slightly here from the analysis presented in the thesis, by thinning the posterior sample by a factor of 10 (i.e. I keep only 1 out of 10 draws). This is done to reduce memory usage.


```r
iter <- 1e5
thinning <- 10
numDraws <- floor(iter / thinning)
warmUpRounds <- c(3, 4, 4)
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
```

### Run RWMH

```r
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
      stepsize = 0.006,  # stepsize1,
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
        iter = 5e4,
        stepsize = .2,
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
```

```
#> Running RWMH for model:  DPMP1-I 
#> Warm up 1 
#> RWMH - Accept ratio is: 0.698 
#> Warm up 2 
#> RWMH - Accept ratio is: 0.515 
#> Warm up 3 
#> RWMH - Accept ratio is: 0.690 
#> RWMH - Accept ratio is: 0.446 
#> RWMH Time:  35.21421  seconds
#>           A      B     C1     C2     C3    w1     w2    w3    w4
#> ESSs 2042.8 1066.9 1740.9 2079.3 2049.9 133.8 1607.6 379.8 143.4
#> 
#> Running RWMH for model:  DPMP1-H 
#> Warm up 1 
#> RWMH - Accept ratio is: 0.623 
#> Warm up 2 
#> RWMH - Accept ratio is: 0.547 
#> Warm up 3 
#> RWMH - Accept ratio is: 0.676 
#> RWMH - Accept ratio is: 0.344 
#> RWMH Time:  54.76792  seconds
#>           A    B     C1     C2     C3     w1     w2     w3     w4
#> ESSs 1330.7 1326 1256.8 1424.9 1465.6 1911.8 1693.5 1793.3 3788.3
#> 
#> Running RWMH for model:  DPMP1-Inv 
#> Warm up 1 
#> RWMH - Accept ratio is: 0.491 
#> Warm up 2 
#> RWMH - Accept ratio is: 0.546 
#> Warm up 3 
#> RWMH - Accept ratio is: 0.681 
#> RWMH - Accept ratio is: 0.384 
#> RWMH Time:  54.85778  seconds
#>           A      B     C1     C2     C3     w1     w2   w3     w4
#> ESSs 1999.5 1238.5 1628.4 1567.7 1634.1 2410.9 1563.7 2655 3033.7
#> 
#> Running RWMH for model:  DPMP2-I 
#> Warm up 1 
#> RWMH - Accept ratio is: 0.659 
#> Warm up 2 
#> RWMH - Accept ratio is: 0.581 
#> Warm up 3 
#> RWMH - Accept ratio is: 0.491 
#> Warm up 4 
#> RWMH - Accept ratio is: 0.402 
#> RWMH - Accept ratio is: 0.279 
#> RWMH Time:  41.10792  seconds
#>        A1    A2   B1   B2   C1    C2  w1   w2  w3 w4
#> ESSs 64.2 379.9 18.7 54.1 37.8 296.2 2.7 30.4 2.1  2
#> 
#> Running RWMH for model:  DPMP2-H 
#> Warm up 1 
#> RWMH - Accept ratio is: 0.308 
#> Warm up 2 
#> RWMH - Accept ratio is: 0.402 
#> Warm up 3 
#> RWMH - Accept ratio is: 0.556 
#> Warm up 4 
#> RWMH - Accept ratio is: 0.376 
#> RWMH - Accept ratio is: 0.423 
#> RWMH Time:  103.7775  seconds
#>          A1   A2     B1     B2     C1     C2  w1     w2    w3    w4
#> ESSs 1729.2 1906 1052.6 1782.6 1291.7 1920.4 195 1745.6 328.7 293.7
#> 
#> Running RWMH for model:  DPMP2-Inv 
#> Warm up 1 
#> RWMH - Accept ratio is: 0.122 
#> Warm up 2 
#> RWMH - Accept ratio is: 0.560 
#> Warm up 3 
#> RWMH - Accept ratio is: 0.680 
#> Warm up 4 
#> RWMH - Accept ratio is: 0.203 
#> RWMH - Accept ratio is: 0.132 
#> RWMH Time:  74.93827  seconds
#>         A1    A2    B1    B2  C1    C2     w1    w2     w3    w4
#> ESSs 444.6 269.6 365.2 206.9 259 363.4 1462.6 387.8 1368.9 766.9
#> 
#> Running RWMH for model:  DPMP3-I 
#> Warm up 1 
#> RWMH - Accept ratio is: 0.397 
#> Warm up 2 
#> RWMH - Accept ratio is: 0.420 
#> Warm up 3 
#> RWMH - Accept ratio is: 0.002 
#> Warm up 4 
#> RWMH - Accept ratio is: 0.845 
#> RWMH - Accept ratio is: 0.372 
#> RWMH Time:  43.11761  seconds
#>         A1     A2     A3   B1   B2    B3    C1     C2    w1     w2    w3
#> ESSs 466.8 1017.9 1232.6 76.9 36.8 964.7 108.7 1528.3 750.2 1034.6 226.2
#>         w4
#> ESSs 450.2
#> 
#> Running RWMH for model:  DPMP3-H 
#> Warm up 1 
#> RWMH - Accept ratio is: 0.065 
#> Warm up 2 
#> RWMH - Accept ratio is: 0.694 
#> Warm up 3 
#> RWMH - Accept ratio is: 0.512 
#> Warm up 4 
#> RWMH - Accept ratio is: 0.486 
#> RWMH - Accept ratio is: 0.001 
#> RWMH Time:  80.32848  seconds
#>       A1   A2  A3   B1    B2   B3 C1  C2  w1  w2   w3  w4
#> ESSs 4.2 13.2 5.8 35.6 643.5 15.7 11 2.7 5.9 3.1 19.3 4.8
#> 
#> Running RWMH for model:  DPMP3-Inv 
#> Warm up 1 
#> RWMH - Accept ratio is: 0.156 
#> Warm up 2 
#> RWMH - Accept ratio is: 0.462 
#> Warm up 3 
#> RWMH - Accept ratio is: 0.393 
#> Warm up 4 
#> RWMH - Accept ratio is: 0.370 
#> RWMH - Accept ratio is: 0.001 
#> RWMH Time:  58.46403  seconds
#>      A1  A2   A3   B1   B2   B3  C1  C2  w1  w2  w3  w4
#> ESSs  4 8.1 14.3 15.4 29.9 37.2 6.8 3.2 5.9 7.7 1.3 3.5
```
## Model comparisons



```r
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
```

### Compute marginals and Bayesian Information criteria (BICs)

```r
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
#> Warning: 4500 of the 4500 log_prob() evaluations on the warp-transformed
#> posterior draws produced -Inf/Inf.
#> Warning: 4500 of the 4500 log_prob() evaluations on the warp-transformed
#> proposal draws produced -Inf/Inf.
#> Error in out[!from@positive] <- -out[!from@positive]: NAs are not allowed in subscripted assignments
```


Table: Marginals, log-likelihoods and BICS

             MarginalLikelihood   LogLikelihood   BIC
----------  -------------------  --------------  ----
DPMP1-I                       0               0     0
DPMP1-H                       0               0     0
DPMP1-Inv                     0               0     0
DPMP2-I                       0               0     0
DPMP2-H                       0               0     0
DPMP2-Inv                     0               0     0
DPMP3-I                       0               0     0
DPMP3-H                       0               0     0
DPMP3-Inv                     0               0     0

### Compute Bayes Factors (BFs)

```r
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
```

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

## Some Parameter Statistics

```r
selectedModels <- c("DPMP1-Inv", "DPMP2-Inv", "DPMP3-Inv")
for (model in selectedModels) {
  summary_ <- summary(coda::mcmc(drawsRWMHLst[[model]]))$statistics
  print(knitr::kable(summary_, caption = model))
}
```



Table: DPMP1-Inv

             Mean          SD    Naive SE   Time-series SE
---  ------------  ----------  ----------  ---------------
A       0.0273618   0.0044397   0.0000444        0.0000993
B       0.7805572   0.0724936   0.0007249        0.0020600
C1      0.2277380   0.0934401   0.0009344        0.0023156
C2      0.4370118   0.8335829   0.0083358        0.0210535
C3      0.1195851   0.1074581   0.0010746        0.0026582
w1     -5.4222511   0.0277406   0.0002774        0.0005650
w2    -10.0996400   0.2666992   0.0026670        0.0067445
w3     -5.4739584   0.0297156   0.0002972        0.0005767
w4     -5.7395821   0.0601529   0.0006015        0.0010921


Table: DPMP2-Inv

             Mean          SD    Naive SE   Time-series SE
---  ------------  ----------  ----------  ---------------
A1      0.0271803   0.0041726   0.0000417        0.0001979
A2      0.0062192   0.0054551   0.0000546        0.0003323
B1      0.7794280   0.0672685   0.0006727        0.0035198
B2      0.3116239   0.3843931   0.0038439        0.0267252
C1      0.2327466   0.0966190   0.0009662        0.0060036
C2      0.5153912   0.8145735   0.0081457        0.0427312
w1     -5.4260844   0.0281035   0.0002810        0.0007349
w2    -10.0956840   0.2648900   0.0026489        0.0134516
w3     -5.4725983   0.0330265   0.0003303        0.0008926
w4     -5.7372557   0.0590894   0.0005909        0.0021337


Table: DPMP3-Inv

            Mean          SD    Naive SE   Time-series SE
---  -----------  ----------  ----------  ---------------
A1     0.0109051   0.0048049   0.0000480        0.0024170
A2     0.0089534   0.0044570   0.0000446        0.0015648
A3     0.0200323   0.0017893   0.0000179        0.0004726
B1     0.9880700   0.0100798   0.0001008        0.0025648
B2     0.9906126   0.0111670   0.0001117        0.0020415
B3     0.9772344   0.0139665   0.0001397        0.0022900
C1    -0.4163860   1.7647617   0.0176476        0.6772413
C2     1.2543645   0.7175671   0.0071757        0.4018811
w1    -5.0940337   0.0230581   0.0002306        0.0095300
w2    -9.5191345   0.1032546   0.0010325        0.0372098
w3    -5.0783797   0.0460255   0.0004603        0.0402249
w4    -5.1766137   0.0489235   0.0004892        0.0262874

## Plots

### Posterior of intensity

```r
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
intensityDraws <- exp(intensityDraws)
```

### Mean log intensity plots

```r
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
```

<div class="figure" style="text-align: center">
<img src="bayesgas-modelcomp_files/figure-html/unnamed-chunk-22-1.png" alt="Mean Log Intensities"  />
<p class="caption">Mean Log Intensities</p>
</div>


### Joint distribution plots

```r
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
```

<img src="bayesgas-modelcomp_files/figure-html/unnamed-chunk-23-1.png" style="display: block; margin: auto;" />

### Highest Posterior Density (HPD) intensity plots

```r
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
```

<div class="figure" style="text-align: center">
<img src="bayesgas-modelcomp_files/figure-html/unnamed-chunk-24-1.png" alt="HPD plots"  />
<p class="caption">HPD plots</p>
</div>

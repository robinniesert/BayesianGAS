context("Test DPMP model functions")

# Simulate data
set.seed(100)
kTransitionTypes <- c("IGtoSIG", "IGtoD", "SIGtoIG", "SIGtoD")
simParams <- c(
  A = c(0.07, 0.04),
  B = c(0.99, 0.96),
  C = c(0.9, 1.2),
  W = c(-5.42, -10.07, -5.45, -5.71)
)
dpmpSimModel <- new(DPMP, simParams, 2, -1.)
sims <- dpmpSimModel$Simulate(750, 1000, 400, rep(0, 2))
tauRaw <- as.integer(cumsum(sims[,5]))
transitionsRaw <- data.frame(sims[,c(1:4)])
possibleTranstionsRaw <- data.frame(sims[,c(6:9)])
colnames(possibleTranstionsRaw) <- kTransitionTypes
colnames(transitionsRaw) <- kTransitionTypes
transitions <- as.matrix(
  aggregate(transitionsRaw, by = list(tauRaw), FUN = sum))
tau <- transitions[, 'Group.1']
diffTau <- diff(tau)
transitions <- transitions[-1, kTransitionTypes]
endTau <- tail(tau, n = 1)
possibleTranstions <- as.matrix(
  aggregate(possibleTranstionsRaw, by = list(tauRaw),
            FUN = function(x){tail(x, n = 1)}))[-1, kTransitionTypes]
y <- cbind(transitions, diffTau, possibleTranstions)
numObs <- dim(y)[1]

# Set up test model
testParams1 <- c(
  A = .05,
  B = .95,
  C = rep(.5, 3),
  W = rep(-2, 4)
)
testParams2 <- c(
  A = rep(.025, 2),
  B = rep(.95, 2),
  C = rep(.5, 2),
  W = rep(-5, 4)
)
testParams3 <- c(
  A = rep(.025, 3),
  B = rep(.95, 3),
  C = rep(.5, 2),
  W = rep(-5, 4)
)
dpmpOneH <- new(DPMP, testParams1, 1, -.5)
dpmpTwoI <- new(DPMP, testParams2, 2, 0.)
dpmpThreeInv <- new(DPMP, testParams3, 3, -1.)

# Set up tests
test_that("Simulate function works", {
  sims <- dpmpTwoI$Simulate(10, 10, 10, rep(0, 2))
  testSims <- c(0., 1., 0., 0., 5.964221, 8., 8., 8., 8., -4.986267, -4.986267,
                -4.998853, -4.972534)
  expect_equal(sims[10,], testSims, tolerance = 1e-6)
})

test_that("Log-Likelihood function produces correct output.", {
  expect_equal(dpmpOneH$LogL(y, 0), -3853.512374)
  expect_equal(dpmpTwoI$LogL(y, rep(0, 2)), -2541.813887)
  expect_equal(dpmpThreeInv$LogL(y, rep(0, 3)), -2537.172501)
})

test_that("Filter function produces correct outputs.", {
  predOut1 <- dpmpOneH$IntensityFilter(y, 0, FALSE)[-(1:(numObs - 10)),1]
  predOut2 <- dpmpTwoI$IntensityFilter(y, rep(0, 2), TRUE)[-(1:(numObs - 10)),2]
  predOut3 <-
    dpmpThreeInv$IntensityFilter(y, rep(0, 3), TRUE)[-(1:(numObs - 10)),4]
  testOut1 <- c(0.007905318, 0.009250078, 0.011724207, 0.008381782, 0.010253551,
                0.009144963, 0.006895741, 0.007156569, 0.007479168, 0.008228214)
  testOut2 <- c(-5.189722, -5.169974, -5.122544, -5.185187, -5.159093,
                -5.169257, -5.252315, -5.275162, -5.299946, -5.307735)
  testOut3 <- c(-5.492686, -5.432632, -5.252010, -5.436352, -5.303152,
                -5.351618, -5.693915, -5.630986, -5.782185, -5.825152)
  expect_equal(predOut1, testOut1, tolerance = 1e-6)
  expect_equal(predOut2, testOut2, tolerance = 1e-6)
  expect_equal(predOut3, testOut3, tolerance = 1e-6)
})

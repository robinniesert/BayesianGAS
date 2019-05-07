context("Test DPMP model functions")

# Import data
data(CreditRatings)
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
test_that("Log-Likelihood function produces correct output.", {
  expect_equal(dpmpOneH$LogL(y, 0), -19840.38781674)
  expect_equal(dpmpTwoI$LogL(y, rep(0, 2)), -17191.38684532)
  expect_equal(dpmpThreeInv$LogL(y, rep(0, 3)), -16679.66676794)
})

test_that("Filter function produces correct outputs.", {
  predOut1 <- dpmpOneH$IntensityFilter(y, 0, FALSE)[-(1:(numObs - 10)),1]
  predOut2 <- dpmpTwoI$IntensityFilter(y, rep(0, 2), TRUE)[-(1:(numObs - 10)),2]
  predOut3 <-
    dpmpThreeInv$IntensityFilter(y, rep(0, 3), TRUE)[-(1:(numObs - 10)),4]
  testOut1 <- c(0.012594198, 0.012725473, 0.010828525, 0.009321664, 0.009015872,
                0.009059863, 0.007856416, 0.008531783, 0.009438477, 0.008463666)
  testOut2 <- c(-5.167010, -5.131307, -5.156764, -5.184628, -5.208392,
                -5.216268, -5.254251, -5.243432, -5.221405, -5.239631)
  testOut3 <- c(-4.851033, -4.795654, -4.882212, -5.017942, -5.118858,
                -5.238489, -5.353825, -5.353975, -5.280178, -5.331489)
  expect_equal(predOut1, testOut1, tolerance = 1e-6)
  expect_equal(predOut2, testOut2, tolerance = 1e-6)
  expect_equal(predOut3, testOut3, tolerance = 1e-6)
})

test_that("Simulate function works", {
  set.seed(100)
  sims <- dpmpTwoI$Simulate(10, 10, 10, rep(0, 2))
  testSims <- c(0., 0., 1., 0., 2.195443, 6., 6., 9., 9., -4.996169, -4.996169,
                -5.000163, -4.992337)
  expect_equal(sims[10,], testSims, tolerance = 1e-6)
})

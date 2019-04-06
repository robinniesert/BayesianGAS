context("Test Beta-Gen-t-EGARCH model functions")

# Set constants
kFiveYrsIdx <- 3769
kNumParams <- 6

# Import data
data(SP500)
times <- as.Date(rev(spData$Date)[-(1:kFiveYrsIdx)])
returns <- ts(
  diff(rev(spData$Adj.Close))[-(1:(kFiveYrsIdx - 1))],
  start = c(2012, 4, 16),
  frequency = 254
)

# Set data attributes
exKurt <- moments::kurtosis(returns) - 3
ltScale <- log(sd(returns) * sqrt((6 / exKurt + 2) / (6 / exKurt + 4)))  # Based
# on student-t second and fourth order moments
numObs <- length(returns)

# Set up test model
testParams <- c(
  omega = 2.42996575,
  A = 0.09306954,
  B = 0.90342348,
  mu = 1.15188494,
  etaBar = 0.03713750,
  upsilon = 1.42561595
)
betaGen <- new(BetaGenTEGARCH, testParams)

# Set up tests
test_that("Log-Likelihood function produces correct output.", {
  predOut <- betaGen$LogL(returns, ltScale)
  testOut <- -5070.29285803
  expect_equal(predOut, testOut)
})

test_that("Filter function produces correct outputs.", {
  predOut <- betaGen$Filter(returns, ltScale)[-(1:(numObs - 10))]
  testOut <- c(2.117198, 2.077816, 2.020446, 2.012800, 2.098203, 2.288550,
               2.432358, 2.395860, 2.338682, 2.429323)
  expect_equal(predOut, testOut, tolerance = 1e-6)
})

test_that("Volatility filter function produces correct outputs.", {
  predOut <- betaGen$VolFilter(returns, ltScale)[-(1:(numObs - 10))]
  testOut <- c(10.037776, 9.650151, 9.112099, 9.042694, 9.848902, 11.913924,
               13.756548, 13.263517, 12.526411, 13.714863)
  expect_equal(predOut, testOut, tolerance = 1e-6)
})

test_that("Gradient of log likelihood is computed correctly.", {
  numDer <- function(x, eps=1e-5){
    grad <- rep(0, kNumParams)
    for (i in 1:kNumParams) {
      tmpX <- x
      tmpX[i] <- x[i] - eps
      fDown <- betaGen$LogLWPar(tmpX, returns, ltScale)
      tmpX[i] <- x[i] + eps
      fUp <- betaGen$LogLWPar(tmpX, returns, ltScale)
      grad[i] <- (fUp - fDown) / (2 * eps)
    }
    grad
  }
  predOut <- betaGen$GradLogLWPar(testParams, returns, ltScale)
  testOut <- numDer(testParams)
  expect_equal(predOut, testOut, tolerance = 1e-5)
})

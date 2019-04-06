context("Test Beta-t-EGARCH model functions")

# Set constants
kFiveYrsIdx <- 3769
kNumParams <- 5

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
  omega = ltScale,
  A = 0.05,
  B = 0.95,
  mu = 0.7746872,
  NuBar = 0.2
)
betaT <- new(BetaTEGARCH, testParams)

# Set up tests
test_that("Log-Likelihood function produces correct output.", {
  predOut <- betaT$LogL(returns, ltScale)
  testOut <- -5084.65012852
  expect_equal(predOut, testOut)
})

test_that("Filter function produces correct outputs.", {
  predOut <- betaT$Filter(returns, ltScale)[-(1:(numObs - 10))]
  testOut <- c(2.182745, 2.154593, 2.122868, 2.106291, 2.140993, 2.240420,
               2.340992, 2.328328, 2.299895, 2.369732)
  expect_equal(predOut, testOut, tolerance = 1e-6)
})

test_that("Volatility filter function produces correct outputs.", {
  predOut <- betaT$VolFilter(returns, ltScale)[-(1:(numObs - 10))]
  testOut <- c(11.45192, 11.13403, 10.78635, 10.60902, 10.98363, 12.13183,
               13.41542, 13.24660, 12.87526, 13.80657)
  expect_equal(predOut, testOut, tolerance = 1e-6)
})

test_that("Gradient of log likelihood is computed correctly.", {
  numDer <- function(x, eps=1e-5){
    grad <- rep(0, kNumParams)
    for (i in 1:kNumParams) {
      tmpX <- x
      tmpX[i] <- x[i] - eps
      fDown <- betaT$LogLWPar(tmpX, returns, ltScale)
      tmpX[i] <- x[i] + eps
      fUp <- betaT$LogLWPar(tmpX, returns, ltScale)
      grad[i] <- (fUp - fDown) / (2 * eps)
    }
    grad
  }
  predOut <- betaT$GradLogLWPar(testParams, returns, ltScale)
  testOut <- numDer(testParams)
  expect_equal(predOut, testOut, tolerance = 1e-5)
})

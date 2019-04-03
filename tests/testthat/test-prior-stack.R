context("Test prior stack Rcpp class")

priorStrs <- c("ImproperUniform", "ImproperUniform")
priorParams <- list(c(-10, Inf), c(-Inf, -100))

priorStack <- new(
  PriorStack,
  priorStrs,
  priorParams
)

priorStackWIndex <- new(
  PriorStack,
  priorStrs,
  priorParams,
  matrix(c(0, 2, 3, 4), nrow = 2, ncol = 2, byrow = TRUE)
)

test_that("Simple PriorStack methods and fields work properly.", {
  expect_equal(priorStack$PriorParams, priorParams)
  expect_equal(priorStack$LogPriors(c(0, -101)), 0.)
  expect_equal(priorStack$LogPriors(c(-10, -101)), -Inf)
  expect_equal(priorStack$LogPriorsWPar(c(-10, -101), priorParams), -Inf)
  expect_equal(priorStack$GradLogPriors(c(0, -101)), c(0., 0.))
  expect_error(
    priorStack$LogPriors(c(-10)),
    "Params for evaluation is not congruent with the PriorStack spec."
  )
})

test_that("PriorStack methods and fields work properly w priorToParamIndex.", {
  expect_equal(priorStackWIndex$PriorParams, priorParams)
  expect_equal(priorStackWIndex$LogPriors(c(0, 10, -9.99,  -101, -1e6)), 0.)
  expect_equal(priorStackWIndex$LogPriors(c(0, 10, -9.99,  -100, -1e6)), -Inf)
  expect_equal(
    priorStackWIndex$GradLogPriors(c(0, 10, -9.99,  -101, -1e6)),
    rep(0, 5)
  )
  expect_error(
    new(
      PriorStack,
      priorStrs,
      priorParams,
      matrix(c(0, 2, 3, 4, 5, 6), nrow = 3, ncol = 2, byrow = TRUE)
    ),
    "Specify proper format for param index, see doc."
  )
  expect_error(
    priorStack$LogPriors(c(0, 10, -9.99,  -100, -1e6, -10)),
    "Params for evaluation is not congruent with the PriorStack spec."
  )
})

context("Test Create Model functions from utils.R")

test_that("CreateModelFromStr function works",{
  expect_error(CreateModelFromStr(new(BetaGenTEGARCH)),
               "Pass string for model argument.")
  expect_error(CreateModelFromStr("not a model"),
               "Specify an implemented model, see doc for available models.")
})

test_that("CreateModel function works",{
  expect_is(CreateModel("BetaGenTEGARCH"), "Rcpp_GASModel")
  expect_equal(CreateModel("BetaGenTEGARCH")$Name, "BetaGenTEGARCH")
  expect_is(CreateModel(new(BetaGenTEGARCH)), "Rcpp_BetaGenTEGARCH")
  expect_is(CreateModel("BetaTEGARCH"), "Rcpp_GASModel")
  expect_equal(CreateModel("BetaTEGARCH")$Name, "BetaTEGARCH")
  expect_is(CreateModel(new(BetaTEGARCH)), "Rcpp_BetaTEGARCH")
  expect_error(CreateModel("not a model"), "Specified model is invalid.")
  expect_error(CreateModel(NULL), "Specified model is invalid.")

  testParams <- c(1, 2, 3, 4, 5, 6)
  modelFromStr <- CreateModel("BetaGenTEGARCH", params = testParams)
  model <- CreateModel(new(BetaGenTEGARCH), params = testParams)
  expect_equal(modelFromStr$Params, testParams)
  expect_equal(model$Params, testParams)
})

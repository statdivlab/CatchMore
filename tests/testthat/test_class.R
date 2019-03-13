library(CatchAll)
library(breakaway)
library(testthat)
context("All models have the correct format")

data(apples)

test_that("The model estimators have correct classes", {

  expect_is(Poisson_model(apples, 50), c("alpha_estimate", "list"))
  expect_is(geometric_model(apples, 50), c("alpha_estimate", "list"))
  expect_is(two_geometric_model(apples, 50), c("alpha_estimate", "list"))
  expect_is(three_geometric_model(apples, 70), c("alpha_estimate", "list"))


})





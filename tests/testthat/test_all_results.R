library(CatchAll)
library(breakaway)
library(testthat)
context("All models have the correct format")

data(apples)

test_that("Analysis of the apples dataset gives similar results", {

  apple_results_new <- all_parametric_model(apples)
  apple_results_old <- read.csv("./tests/apple_Analysis.csv")

  apple_results_new <- apple_results_new[, c("Model", "tau", "Est", "SE", "AICc", "GOF5")]
  apple_results_old <- apple_results_old[, c("Model", "Cutoff", "Estimated.Total.Sp", "SE",
                                             "AICc", "GOF5")]

  names(apple_results_new) <- c("Model", "tau", "New.Est", "New.SE", "New.AICc", "New.GOF5")
  apple_results_new$Model <- as.character(apple_results_new$Model)

  names(apple_results_old) <- c("Model", "tau", "Old.Est", "Old.SE", "Old.AICc", "Old.GOF5")
  apple_results_old$Model <- as.character(apple_results_old$Model)

  apples_merged <- merge(apple_results_new, apple_results_old,
                         by = c("Model", "tau"))


  ## test whether the estimates, standard errors, AICc and GOF5 are equal within a precision limit.
  for (ii in 1:nrow(apples_merged)) {
    expect_lte(abs(apples_merged$New.Est[ii] - apples_merged$Old.Est[ii]),
               0.05 * apples_merged$Old.Est[ii])
    expect_lte(abs(apples_merged$New.SE[ii] - apples_merged$Old.SE[ii]),
               0.05 * apples_merged$Old.SE[ii])
    expect_lte(abs(apples_merged$New.AICc[ii] - apples_merged$Old.AICc[ii]),
               0.05 * apples_merged$Old.AICc[ii])
    expect_lte(abs(apples_merged$New.GOF5[ii] - apples_merged$Old.GOF5[ii]),
               max(0.05 * apples_merged$Old.GOF5[ii], 0.01))
  }
})

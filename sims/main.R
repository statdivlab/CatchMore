# This is the main simulator file

library(simulator) # this file was created under simulator version 0.2.0
library(CatchAll)
source("model_functions.R")
source("method_functions.R")
source("eval_functions.R")

## @knitr init

name_of_simulation <- "catchall-with-poisson-or-geometric-model"

## @knitr main

poisson_sim <- new_simulation(name = "catchall-with-poisson-model",
                      label = "Richness estimation with CatchAll") %>%
  generate_model(poisson_counts, seed = 123,
                 n = as.list(1:10 * 100),
                 lambda = 0.7,
                 vary_along = "n") %>%
  simulate_from_model(nsim = 50) %>%
  run_method(list(catchall_best)) %>%
  evaluate(list(bias_richness, se_richness))

## @knitr plots

plot_eval_by(poisson_sim, "bias_richness", varying = "n")

## @knitr tables

tabulate_eval(poisson_sim, "bias_richness", output_type = "markdown",
              format_args = list(digits = 1))

## geometric model
geom_sim <- new_simulation(name = "catchall-with-geometric-model",
                              label = "Richness estimation with CatchAll") %>%
  generate_model(geom_counts, seed = 123,
                 n = as.list(1:10 * 100),
                 lambda = 0.7,
                 vary_along = "n") %>%
  simulate_from_model(nsim = 50) %>%
  run_method(list(catchall_best)) %>%
  evaluate(list(bias_richness, se_richness))

## @knitr plots

plot_eval_by(geom_sim, "bias_richness", varying = "n")

## @knitr tables

tabulate_eval(geom_sim, "bias_richness", output_type = "markdown",
              format_args = list(digits = 1))
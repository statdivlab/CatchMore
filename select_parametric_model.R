source("empty_functions.R")
source("two_mixed_exp.R")

library(breakaway)
library(dplyr)
library(magrittr)
library(parallel)
input_data <- apples
select_parametric_model <- function(input_data, 
                                    parallel = T,
                                    control = list(ncores = ceiling(detectCores()/2))) {
  input_data <- breakaway::convert(input_data)
  ii <- input_data$index
  
  
 
  all_results <- lapply(ii[3:max(length(ii) - 3, 10)], function (tau) {
    poisson_tau <- Poisson_model(input_data, cutoff = tau)
    geometric_tau <- geometric_model(input_data, cutoff = tau)
    two_mixed_geom_tau <- two_geometric_model(input_data, cutoff = tau)
    
    tau_tab <- tibble(tau = rep(tau, 3), 
                      Model = c("Poisson", "Geometric", "Two-component Mixture Exponential"),
                      Est = c(poisson_tau$estimate, geometric_tau$estimate, two_mixed_geom_tau$estimate),
                      SE = c(poisson_tau$error, geometric_tau$error, two_mixed_geom_tau$error),
                      AICc = c(poisson_tau$AICc, geometric_tau$AICc, two_mixed_geom_tau$AICc),
                      GOF0 = c(poisson_tau$GOF0, geometric_tau$GOF0, two_mixed_geom_tau$GOF0),
                      GOF5 = c(poisson_tau$GOF5, geometric_tau$GOF5, two_mixed_geom_tau$GOF5)) %>%
      filter(!is.na(Est))
    return(tau_tab)
  }) 
  
  all_results_tib <- all_results %>% 
    data.table::rbindlist() %>%
    arrange(Model, tau)
  
  ## Model selection
  all_results_tib %>% 
   # filter(GOF5 >= 0.01) %>%
    group_by(tau) %>%
    filter(AICc == min(AICc))
}
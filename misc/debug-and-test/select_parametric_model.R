source("empty_functions.R")
source("two_mixed_exp.R")
source("three_mixed_exp.R")
library(breakaway)


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
    three_mixed_geom_tau <- three_geometric_model(input_data, cutoff = tau)
    tau_tab <- data.frame(tau = tau, 
                      Model = c("Poisson", "SingleExp", "TwoMixedExp", "ThreeMixedExp"),
                      Est = c(poisson_tau$estimate, geometric_tau$estimate, two_mixed_geom_tau$estimate, three_mixed_geom_tau$estimate),
                      SE = c(poisson_tau$error, geometric_tau$error, two_mixed_geom_tau$error, three_mixed_geom_tau$error),
                      AICc = c(poisson_tau$AICc, geometric_tau$AICc, two_mixed_geom_tau$AICc, three_mixed_geom_tau$AICc),
                      GOF0 = c(poisson_tau$GOF0, geometric_tau$GOF0, two_mixed_geom_tau$GOF0, three_mixed_geom_tau$GOF0),
                      GOF5 = c(poisson_tau$GOF5, geometric_tau$GOF5, two_mixed_geom_tau$GOF5, three_mixed_geom_tau$GOF5)) 

    tau_tab <- tau_tab[!is.na(tau_tab$Est),]
    return(tau_tab)
  }) 
  
  all_results_tib <- all_results[[1]]
  
  for(i in 2:length(all_results)) {
    all_results_tib <- rbind(all_results_tib, all_results[[i]])
  }
  
  all_results_tib <- all_results_tib[order(all_results_tib$Model, all_results_tib$tau),]
  
  ## Model selection

}

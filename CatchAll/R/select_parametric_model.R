

#' Multiple parametric models for estimating species richness with various cutoff values
#'
#' This function generates the species richness estimates with Poisson model, geometric model,
#' two-component geometric mixture model and three-component geometric model.
#'
#'
#' @param input_data An input type that can be processed by convert()
#' @param parallel A logic scalar that decides whether the richness estimates are computed with parallelization.
#' Currently an error is returned for Windows.
#' @param control A list containing an integer \code{ncores} that specifies the number of cores to use
#' in parallelization when \code{parallel == T}
#'
#' @return A data frame displaying the point estimates, standard errors, AICc, GOF0 and GOF5 for different
#' parametric models and cutoffs.
#'
#' @export
#'
#' @examples
#' library(breakaway)
#' data(apples)
#' all_parametric_model(apples)

#' @export
all_parametric_model <- function(input_data,
                                    parallel = F,
                                    control = list(ncores = ceiling(detectCores()/2))) {
  ii <- input_data$index
  input_data <- convert(input_data)
  if(parallel == T) {
    cl <- makeCluster(control$ncores)
    clusterExport(cl, c("input_data", "GOF", "geometric_model", "Poisson_model",
                        "two_geometric_model", "two_mixed_exp_lld", "two_mixed_exp_init",
                        "two_mixed_exp_se", "two_mixed_exp_EM", "three_geometric_model",
                        "three_mixed_exp_lld", "three_mixed_exp_init", "three_mixed_exp_se",
                        "three_mixed_exp_EM", "ii"))

    all_results <- mclapply(ii[3:max(length(ii) - 3, 10)], function (tau) {
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
    }, mc.cores = control$ncores)

    stopCluster(cl)
  } else{
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
  }


  all_results_tib <- all_results[[1]]

  for(i in 2:length(all_results)) {
    all_results_tib <- rbind(all_results_tib, all_results[[i]])
  }

  all_results_tib <- all_results_tib[order(all_results_tib$Model, all_results_tib$tau),]

  ## Model selection

}



#' Multiple parametric models for estimating species richness with various cutoff values
#'
#' This function generates the species richness estimates with Poisson model, geometric model,
#' two-component geometric mixture model and three-component geometric model.
#'
#'
#' @param input_data An input type that can be processed by convert()
#' @param parallel A logic scalar that decides whether the richness estimates are computed with parallelization.
#' Currently an error is returned for Windows.
#' @param tau_range A two-dimensional vector of positive intergers specifing the range of the cutoffs you'd like to vary within
#' @param control A list containing an integer \code{ncores} that specifies the number of cores to use
#' in parallelization when \code{parallel == T}
#'
#' @return A data frame displaying the point estimates, standard errors, AICc, GOF0 and GOF5 for different
#' parametric models and cutoffs.
#'
#' @import dplyr
#'
#'
#' @examples
#' library(breakaway)
#' data(apples)
#' all_parametric_model(apples)

#' @export
all_parametric_model <- function(input_data,
                                 parallel = F,
                                 tau_range = NULL,
                                 control = list(ncores = ceiling(detectCores()/2))) {
  options(warn = -1)

  ii <- input_data$index
  input_data <- convert(input_data)

  if (is.null(tau_range)) {
    tau_range <- c(3, max(length(ii) - 3, 10))
  }
  if(parallel == T) {
    cl <- makeCluster(control$ncores)
    clusterExport(cl, c("input_data", "GOF", "geometric_model", "Poisson_model",
                        "two_geometric_model", "two_mixed_exp_lld", "two_mixed_exp_init",
                        "two_mixed_exp_se", "two_mixed_exp_EM", "three_geometric_model",
                        "three_mixed_exp_lld", "three_mixed_exp_init", "three_mixed_exp_se",
                        "three_mixed_exp_EM", "ii"))

    all_results <- mclapply(tau_range[1]:tau_range[2], function (tau) {
      cat("...")
      cat(tau)
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
                            GOF5 = c(poisson_tau$GOF5, geometric_tau$GOF5, two_mixed_geom_tau$GOF5, three_mixed_geom_tau$GOF5),
                            LwrCB = c(poisson_tau$ci[1], geometric_tau$ci[1], two_mixed_geom_tau$ci[1], three_mixed_geom_tau$ci[1]),
                            UprCB = c(poisson_tau$ci[2], geometric_tau$ci[2], two_mixed_geom_tau$ci[2], three_mixed_geom_tau$ci[2]))

      tau_tab <- tau_tab[!is.na(tau_tab$Est),]
      return(tau_tab)
    }, mc.cores = control$ncores)

    stopCluster(cl)
  } else{
    all_results <- lapply(tau_range[1]:tau_range[2], function (tau) {
      cat(".")
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
                            GOF5 = c(poisson_tau$GOF5, geometric_tau$GOF5, two_mixed_geom_tau$GOF5, three_mixed_geom_tau$GOF5),
                            LwrCB = c(poisson_tau$ci[1], geometric_tau$ci[1], two_mixed_geom_tau$ci[1], three_mixed_geom_tau$ci[1]),
                            UprCB = c(poisson_tau$ci[2], geometric_tau$ci[2], two_mixed_geom_tau$ci[2], three_mixed_geom_tau$ci[2]))

      tau_tab <- tau_tab[!is.na(tau_tab$Est),]
      return(tau_tab)
    })
  }


  all_results_tib <- all_results[[1]]

  for(i in 2:length(all_results)) {
    all_results_tib <- rbind(all_results_tib, all_results[[i]])
  }

  all_results_tib <- all_results_tib[order(all_results_tib$Model, all_results_tib$tau),]
  options(warn = 0)
  return(all_results_tib)

}

#' "Best" parametric models for estimating species richness with various cutoff values
#'
#' This function gives the "best" models among Poisson model, geometric model,
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
#' @import dplyr
#' @examples
#' library(breakaway)
#' data(apples)
#' select_best_models(apples)
#'
#' data(hawaii)
#' select_best_models(hawaii)
#' @export

select_best_models <- function(input_data,
                               parallel = F,
                               tau_range = NULL,
                               control = list(ncores = ceiling(detectCores()/2))) {
  options(warn = -1)
  ii <- input_data$index
  input_data <- convert(input_data)

  all_results_tib <- all_parametric_model(input_data, parallel = parallel, tau_range = tau_range, control = control) %>%
    .[complete.cases(.),]
  ## For each cutoff, evaluate all the models
  # all_results <- lapply(ii[3:max(length(ii) - 3, 10)], function (tau) {
  #   poisson_tau <- Poisson_model(input_data, cutoff = tau)
  #   geometric_tau <- geometric_model(input_data, cutoff = tau)
  #   two_mixed_geom_tau <- two_geometric_model(input_data, cutoff = tau)
  #   three_mixed_geom_tau <- three_geometric_model(input_data, cutoff = tau)
  #   tau_tab <- tibble(tau = tau,
  #                     Model = c("Poisson", "SingleExp", "TwoMixedExp", "ThreeMixedExp"),
  #                     Est = c(poisson_tau$estimate, geometric_tau$estimate, two_mixed_geom_tau$estimate, three_mixed_geom_tau$estimate),
  #                     SE = c(poisson_tau$error, geometric_tau$error, two_mixed_geom_tau$error, three_mixed_geom_tau$error),
  #                     AICc = c(poisson_tau$AICc, geometric_tau$AICc, two_mixed_geom_tau$AICc, three_mixed_geom_tau$AICc),
  #                     GOF0 = c(poisson_tau$GOF0, geometric_tau$GOF0, two_mixed_geom_tau$GOF0, three_mixed_geom_tau$GOF0),
  #                     GOF5 = c(poisson_tau$GOF5, geometric_tau$GOF5, two_mixed_geom_tau$GOF5, three_mixed_geom_tau$GOF5))
  #
  #   tau_tab <- tau_tab[!is.na(tau_tab$Est),]
  #   return(tau_tab)
  # })
  #
  # ## combine all the list
  # all_results_tib <- all_results[[1]]
  #
  # for(i in 2:length(all_results)) {
  #   all_results_tib <- bind_rows(all_results_tib, all_results[[i]])
  # }
  #
  # all_results_tib <- arrange(all_results_tib, Model, tau)

  ## Model selection, flag indicates how stringent the selection criteria are
  flag <- 0

  if (sum(all_results_tib$GOF5 > 0.01) > 0) {
    if (sum(all_results_tib$GOF0 > 0.01) > 0) {
      flag <- 2
    } else {
      flag <- 1
    }
  }


  if (flag == 0) {
    ## all models are filtered out.
    ## In this scenario, we return 4 models with best AICc and GOF5
    bestModels <- all_results_tib %>%
      group_by(tau) %>%
      dplyr::filter(AICc == min(AICc)) %>%
      arrange(desc(GOF5)) %>%
      .[1:4,]

    output <- tibble(Description = c("Best Model 1", "Best Model 2", "Best Model 3", "Best Model 4")) %>%
      bind_cols(bestModels)

  } else if (flag == 1) {
    ## relaxed criteria
    ## Model 2A is the model with the greatest GOF0
    bestModels <- all_results_tib %>%
      dplyr::filter(GOF5 > 0.01) %>%
      group_by(tau) %>%
      dplyr::filter(AICc == min(AICc)) %>%
      group_by()

    bestModel1 <- tibble(Description = "Best Model 1", tau = NA, Model = NA, Est = NA, SE = NA, AICc = NA)

    bestModel2A <- bestModels %>%
      dplyr::filter(GOF0 == max(GOF0)) %>%
      dplyr::filter(tau == max(tau)) %>%
      bind_cols(tibble(Description = c("Best Model 2A")),
                .)

    bestModel2B <- bestModels %>%
      dplyr::filter(tau == max(tau))%>%
      bind_cols(tibble(Description = c("Best Model 2B")),
                .)

    if (any(bestModels$tau <= 10)) {
      bestModel2C <- bestModels %>%
        dplyr::filter(tau <= 10) %>%
        dplyr::filter(tau == max(tau)) %>%
        bind_cols(tibble(Description = c("Best Model 2C")),
                  .)
    } else {
      bestModel2C <- tibble(Description = "Best Model 2C", tau = NA, Model = NA, Est = NA, SE = NA, AICc = NA)
    }

    output <- bind_rows(bestModel1, bestModel2A, bestModel2B, bestModel2C)





  } else if (flag == 2) {
    ## We adopt the most stringent criteria

    bestModels <- all_results_tib %>%
      dplyr::filter(GOF5 > 0.01) %>%
      group_by(tau) %>%
      dplyr::filter(AICc == min(AICc)) %>%
      group_by()

    bestModel1 <- bestModels %>%
      dplyr::filter(GOF0 > 0.01) %>%
      dplyr::filter(tau == max(tau)) %>%
      bind_cols(tibble(Description = c("Best Model 1")),
                .)


    bestModel2A <- bestModels %>%
      dplyr::filter(GOF0 == max(GOF0)) %>%
      dplyr::filter(tau == max(tau)) %>%
      bind_cols(tibble(Description = c("Best Model 2A")),
                .)

    bestModel2B <- bestModels %>%
      dplyr::filter(tau == max(tau))%>%
      bind_cols(tibble(Description = c("Best Model 2B")),
                .)

    if (min(bestModels$tau) <= 10) {
    bestModel2C <- bestModels %>%
      dplyr::filter(tau <= 10) %>%
      dplyr::filter(tau == max(tau)) %>%
      bind_cols(tibble(Description = c("Best Model 2C")),
                .)
    } else {
    bestModel2C <- tibble(Description = "Best Model 2C", tau = NA, Model = NA, Est = NA, SE = NA, AICc = NA)
    }

    output <- bind_rows(bestModel1, bestModel2A, bestModel2B, bestModel2C)


  }

  results <- list(BestModels = output,
                  flag = flag,
                  message = switch(as.character(flag),
                                   "0" = "All GOF5 < 0.01: Return Models with smallest AICc and GOF5",
                                   "1" = "All GOF0 < 0.01: Best Model 1 not available",
                                   "2" = "All best parametric models available"))
  options(warn = 0)
  return(results)
}

#' "Best" parametric model for estimating species richness
#'
#' This function gives the "best" model among Poisson model, geometric model,
#' two-component geometric mixture model and three-component geometric model.
#'
#'
#' @param input_data An input type that can be processed by convert()
#' @param alpha_estimate If \code{alpha_estimate == T}, the result will be returned as an \code{alpha_estimate} object.
#' @param ... Additional arguments for the function \code{select_best_models}
#' @return A data frame displaying the point estimates, standard errors, AICc, GOF0 and GOF5 for different
#' parametric models and cutoffs.
#'
#' @import dplyr
#' @examples
#' library(breakaway)
#' data(apples)
#' best_model(apples)
#'
#' @export
catch_best <- function(input_data, alpha_estimate = T, ...) {
  best_models <- select_best_models(input_data, ...)

  bob <- best_models$BestModels %>%
    dplyr::filter(., complete.cases(.)) %>%
    .[1, ]

  if(alpha_estimate == T) {
    bob_model <- switch(as.character(bob$Model[1]),
                        "Poisson" = "Poisson_model",
                        "SingleExp" = "geometric_model",
                        "TwoMixedExp" = "two_geometric_model",
                        "ThreeMixedExp" = "three_geometric_model")

    bob_alpha <- do.call(bob_model, args = list(input_data = input_data,
                                                cutoff = bob$tau))
    return(bob_alpha)
  } else {
    return(bob)
  }
}
## system.time(BBM <- select_best_models(hawaii))
## user  system elapsed
## 1455.55    0.23 1456.75


#' CatchAll for estimating species richness and model selection
#'
#' This function implements a series of estimator and selects several of them for estimating species richness
#'
#' @param input_data An input type that can be processed by \code{convert()}
#' @param ... Other arguments for \code{select_best_models}
#'
#' @import dplyr
#' @examples
#' library(breakaway)
#' data(apples)
#' catch_all(apples)
#'
#' @export
catch_all <- function(input_data, ...) {

  param_tab <- select_best_models(input_data, ...)
  ## chao1
  chao1_result <- try(chao1(input_data), silent = T)
  if (class(chao1_result)[1] !="try-error") {
    chao1_tab <- tibble(Description = "NonP1", tau = 1, Model = "Chao1", Est = chao1_result$estimate,
                        SE = chao1_result$error, AICc = NA, GOF0 = NA, GOF5 = NA,
                        LwrCB = chao1_result$ci[1], UprCB = chao1_result$ci[2])
  } else {
    chao1_tab <- tibble(Description = "NonP1", tau = 1, Model = "Chao1", Est = NA,
                        SE = NA, AICc = NA, GOF0 = NA, GOF5 = NA,
                        LwrCB = NA, UprCB = NA)
  }


  ## WLRM
  wlrm_result <- try(wlrm_untransformed(input_data), silent = T)
  if (class(wlrm_result)[1] !="try-error") {
    wlrm_tab <- tibble(Description = "WLRM", tau = wlrm_result$other$cutoff, Model = "UnTrans", Est = wlrm_result$estimate,
                        SE = wlrm_result$error, AICc = NA, GOF0 = NA, GOF5 = NA,
                        LwrCB = wlrm_result$ci[1], UprCB = wlrm_result$ci[2])
  } else {
    wlrm_tab <- tibble(Description = "WLRM", tau = wlrm_result$other$cutoff, Model = "UnTrans", Est = NA,
                        SE = NA, AICc = NA, GOF0 = NA, GOF5 = NA,
                        LwrCB = NA, UprCB = NA)
  }


  logwlrm_result <- try(wlrm_transformed(input_data), silent = T)
  if (class(logwlrm_result)[1] !="try-error") {
    logwlrm_tab <- tibble(Description = "WLRM", tau = logwlrm_result$other$cutoff, Model = "LogTrans", Est = logwlrm_result$estimate,
                       SE = logwlrm_result$error, AICc = NA, GOF0 = NA, GOF5 = NA,
                       LwrCB = logwlrm_result$ci[1], UprCB = logwlrm_result$ci[2])
  } else {
    logwlrm_tab <- tibble(Description = "WLRM", tau = logwlrm_result$other$cutoff, Model = "LogTrans", Est = NA,
                       SE = NA, AICc = NA, GOF0 = NA, GOF5 = NA,
                       LwrCB = NA, UprCB = NA)
  }

  kemp_result <- try(kemp(input_data), silent = T)
  if (class(kemp_result)[1] !="try-error") {
    kemp_tab <- tibble(Description = "KEMP", tau = NA, Model = "Kemp", Est = kemp_result$estimate,
                       SE = kemp_result$error, AICc = NA, GOF0 = NA, GOF5 = NA,
                       LwrCB = kemp_result$ci[1], UprCB = kemp_result$ci[2])
  } else {
    kemp_tab <- tibble(Description = "KEMP", tau = NA, Model = "Kemp", Est = NA,
                       SE = NA, AICc = NA, GOF0 = NA, GOF5 = NA,
                       LwrCB = NA, UprCB = NA)
  }


  suppressWarnings( res_tab <- bind_rows(param_tab$BestModels %>% dplyr::filter(!is.na(Est)),
                                         chao1_tab,
                                         wlrm_tab,
                                         logwlrm_tab,
                                         kemp_tab))

  return(res_tab)
}


## convert: transform phyloseq object into frequency count data

## alpha_estimate: return an "alpha_estimate" object

## each function first convert a phyloseq object or otu table into a frequency 
## count table, then give an answer.

## GOF


## O_c: observed counts;
## E_c: theoretic counts;
## param: number of parameters
## bin.tol: bin counts
GOF <- function(O_c, E_c, param, bin.tol = 0) {
  ## bin the counts that are less than five
  ## start with the first cell. bin the cells one by one
  if (bin.tol > 0) {
    for (i in 1:(length(E_c) - 1)) {
      while(length(E_c) > i & E_c[i] < bin.tol) {
        E_c[i] <- E_c[i] + E_c[i+1]
        E_c <- E_c[-(i+1)]
        O_c[i] <- O_c[i] + O_c[i+1]
        O_c <- O_c[-(i+1)]
      }
    }
    
    ## if the last cell contains less than 5, bin the last two together
    if (last(E_c) < bin.tol) {
      E_c[length(E_c) - 1] <- E_c[length(E_c) - 1] + E_c[length(E_c)]
      E_c <- E_c[-length(E_c)]
      O_c[length(O_c) - 1] <- O_c[length(O_c) - 1] + O_c[length(O_c)]
      O_c <- O_c[-length(O_c)]
    }
  }
  ## calculate the goodness-of-fit statistic
  WW <- sum((O_c - E_c)^2/E_c)
  
  ## degree of freedom of the GOF of the poisson model
  vv <- length(E_c) - param - 1
  
  return(list(GOF = WW, df = vv))
}


geometric_model <- function(input_data, cutoff = 10) {
  input_data <- breakaway::convert(input_data)
  included <- input_data %>% filter(index <= cutoff)
  excluded <- input_data %>% filter(index > cutoff)
  if (nrow(included) == 0) {
    included <- input_data
    excluded <- list(index = Inf, frequency = 0)
  }
  
  ii <- included$index
  ff <- included$frequency
  
  c_tau <- ff %>% sum
  c_excluded <- excluded$frequency %>% sum
  
  n_tau <- crossprod(ii, ff) %>% as.numeric
  
  ## parameter in the geometric model 
  theta_hat <- n_tau/c_tau - 1
  
  ## estimated richness
  C_hat_truncated <- (n_tau * c_tau) / (n_tau - c_tau)
  
  C_hat <- C_hat_truncated + c_excluded
  
  ## estimated standard error
  se_hat <- C_hat_truncated / sqrt(n_tau - c_tau)
  
  ## TODO: asymetric confidence interval
  f0 <- C_hat_truncated - c_tau
  d <- ifelse(f0 == 0, 1, exp(1.96 * sqrt(log(1 + se_hat^2/f0^2))))
  C_confint <- c(c_tau + c_excluded + f0/d, c_tau + c_excluded + f0*d)
  
  ## TODO: GOODNESS-OF-FIT
  O_c <- c(included$frequency, 0)
  E_c <- c(C_hat_truncated * dgeom(included$index, prob = 1/(1 + theta_hat)),
           C_hat_truncated * pgeom(max(included$index), 
                                   prob = 1/(1 + theta_hat), lower.tail = F))
  
  GOF0.geom <- GOF(O_c = O_c, E_c = E_c, param = 1, bin.tol = 0)
  GOF5.geom <- GOF(O_c = O_c, E_c = E_c, param = 1, bin.tol = 5)
  
  
  ## AICc
  AICc <- 2 * c_tau / (c_tau - 2) - 2 * sum(log(1:c_tau)) + 
    2 * sapply(1:length(ff), function(i) {
      sum(log(1:ff[i]))
    }) %>% sum -
    2 * sapply(1:length(ff), function(i) {
      ff[i] * (-log(theta_hat) + ii[i] * log(theta_hat/(1+theta_hat)))
    }) %>% sum
  
  alpha_estimate(estimate = C_hat,
                 error = se_hat,
                 estimand = "richness",
                 name = "Geometric Model",
                 type = "parametric",
                 model = "Geometric",
                 frequentist = TRUE,
                 parametric = TRUE,
                 reasonable = FALSE,
                 interval = C_confint,
                 interval_type = "Approximate: log-normal",
                 GOF0 = pchisq(GOF0.geom$GOF, GOF0.geom$df, lower.tail = F) %>% signif(4),
                 GOF5 = pchisq(GOF5.geom$GOF, GOF5.geom$df, lower.tail = F) %>% signif(4),
                 AICc = AICc,
                 other = list(theta_hat = theta_hat,
                              cutoff = cutoff))
  

}



three_geometric_model <- function() {
  
}


##--------------------------------------------------------------
## Weighted least square type estimators


wlrm_transformed <- function() {
  
}

wlrm_untransformed <- function() {
  
}


select_wlrm <- function(input_data) {
  
}


## -----------------------
## kemp estimator


kemp <- function() {
  
}

## ----------------------
## good_turing_richness
## as chao1, it's a lower bound

good_turing_richness <- function() {
  
}


## --------------------------
## chao1

chao1 <- function() {
  
}

## --------------------------
## ace

Poisson_model <- function(input_data, cutoff = 10) {
  input_data <- breakaway::convert(input_data)
  
  ## truncate the data by the cutoff
  included <- input_data %>% filter(index <= cutoff)
  excluded <- input_data %>% filter(index > cutoff)
  
  ## estimate the parameters
  if (nrow(included) == 0) {
    included <- input_data
    excluded <- list(index = Inf, frequency = 0)
  }
  
  ii <- included$index
  ff <- included$frequency
  
  c_tau <- sum(ff)
  c_excluded <- excluded$frequency %>% sum
  
  ## total number of individuals in the truncated data
  nn <- included %$% crossprod(index, frequency) %>% 
    as.numeric
  
  poisson_fn <- function(lambda) {
    (1 - exp(-lambda))/lambda - c_tau/nn
  }
  
  lambda_hat <- uniroot(poisson_fn, c(1e-04, 1e+06))$root
  
  ## estimated number of classes after truncation
  c_hat_truncated <- c_tau/(1 - exp(-lambda_hat))
  C_hat <- c_hat_truncated + c_excluded
  
  ## estimated standard error
  cc_se <- sqrt(c_hat_truncated / (exp(lambda_hat) - 1 - lambda_hat))
  
  ## confidence interval
  f0 <- c_hat_truncated - c_tau
  dd <- ifelse(f0 == 0, 1, exp(1.96 * sqrt(log(1 + cc_se^2/f0))))
  C_confint <- c(c_tau + c_excluded + f0/dd, c_tau + c_excluded + f0 * dd)
  
  ## get expected counts
  E_c <- c(c_hat_truncated * dpois(included$index, lambda = lambda_hat),
           c_hat_truncated * ppois(max(included$index), 
                                   lambda = lambda_hat, lower.tail = F))
  O_c <- c(included$frequency, 0)
  
  ## bin the counts that are less than five
  ## start with the first cell. bin the cells one by one
  GOF0.poisson <- GOF(O_c = O_c, E_c = E_c, param = 1, bin.tol = 0)
  GOF5.poisson <- GOF(O_c = O_c, E_c = E_c, param = 1, bin.tol = 5)
  
  ## AICc
  
  AICc <- 2 * c_tau / (c_tau - 2) - 2 * sum(log(1:c_tau)) + 
    2 * sapply(1:length(ff), function(i) {
      sum(log(1:ff[i]))
    }) %>% sum -
    2 * sapply(1:length(ff), function(i) {
      ff[i] * (ii[i] * log(lambda_hat) - lambda_hat - sum(log(1:ii[i])) -
                 log(1 - exp(-lambda_hat)))
    }) %>% sum
  
  
  ## Return the alpha estimate
  alpha_estimate(estimate = C_hat,
                 error = cc_se,
                 estimand = "richness",
                 name = "Poisson Model",
                 type = "parametric",
                 model = "Poisson",
                 frequentist = TRUE,
                 parametric = TRUE,
                 reasonable = FALSE,
                 GOF0 = pchisq(GOF0.poisson$GOF, GOF0.poisson$df, lower.tail = F) %>% signif(4),
                 GOF5 = pchisq(GOF5.poisson$GOF, GOF5.poisson$df, lower.tail = F) %>% signif(4),
                 AICc = AICc,
                 interval = C_confint,
                 interval_type = "Approximate: log-normal",
                 other = list(lambda_hat = lambda_hat,
                              cutoff = cutoff))
}

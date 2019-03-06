## Two-component mixture of exponential mixed Poisson model
library(breakaway)
library(dplyr)
source("empty_functions.R")
data(apples)
input_data <- apples

three_geometric_model <- function(input_data, cutoff = 10, ...) {
  input_data <- breakaway::convert(input_data)
  
  ## convert the data to frequency table
  included <- input_data[input_data$index <= cutoff,]
  excluded <- input_data[input_data$index > cutoff,]
  
  ii <- included$index
  ff <- included$frequency
  
  
  
  if (nrow(included) == 0) {
    included <- input_data
    excluded <- list(index = Inf, frequency = 0)
  }
  
  c_tau <- sum(included$frequency)
  c_excluded <- sum(excluded$frequency)
  
  
  
  if (length(ii) <= 4) {
    warning("Not enough different frequency counts. We recommend increase the cutoff.")
    return(alpha_estimate(estimate = NA,
                          error = NA,
                          estimand = "richness",
                          name = "Two-component Geometric mixture Model",
                          type = "parametric",
                          model = "Two-component Geometric mixture Model",
                          frequentist = TRUE,
                          parametric = TRUE,
                          reasonable = TRUE,
                          interval = NA,
                          interval_type = "Approximate: log-normal",
                          GOF0 = NA,
                          GOF5 = NA,
                          AICc = NA,
                          other = list(t1 = NA,
                                       t2 = NA,
                                       t3 = NA,
                                       u = NA,
                                       cutoff = cutoff)))
  }
  
  ## EM algorithm for obtaining the parameter estimates
  em_params <- three_mixed_exp_EM(ii, ff)
  u1 <- em_params$u1
  u2 <- em_params$u2
  t1 <- em_params$t1
  t2 <- em_params$t2
  t3 <- em_params$t3
  
  t4 <- (1 + t1) * u1 * t2 * t3 / (t1 * t2 * t3 + u2 * t1 * t3 + u1 * t2 * t3 -
                                     u2 * t1 * t2 + t1 * t2 - u1 * t1 * t2)
  t5 <- (1 + t2) * u2 * t1 * t3 / (t1 * t2 * t3 + u2 * t1 * t3 + u1 * t2 * t3 -
                                     u2 * t1 * t2 + t1 * t2 - u1 * t1 * t2)
  
  
  C_tau <- c_tau / (1 - t4 / (1 + t1) - t5 / (1 + t2) - (1 - t4 - t5) / (1 + t3))
  f0 <- C_tau - c_tau
  C_hat <- C_tau + c_excluded
  ## se
  se_hat <- three_mixed_exp_se(C_tau, t1, t2, t3, t4, t5)
  
  ## confidence interval
  dd <- ifelse(f0 == 0, 1, exp(1.96 * sqrt(log(1 + se_hat^2/f0))))
  
  C_confint <- c(c_tau + c_excluded + f0 / dd, c_tau + c_excluded + f0 * dd)
  ## AICc
  AICc <- 10 * c_tau / (c_tau - 6) -2 *three_mixed_exp_lld(ii, ff, t1, t2, t3, u1, u2, full = T)
  
  ## GOF0 and GOF5
  O_c <- c(ff, 0)
  pp_E <- t4 / t1 * (t1 / (1 + t1))^ii + t5 / t2 * (t2 / (1 + t2))^ii +
    (1 - t4 - t5) / t3 * (t3 / (1 + t3))^ii
  E_c <- c_tau * c(pp_E, 1-sum(pp_E))
  
  GOF0_threemixedexp <- GOF(O_c = O_c, E_c = E_c, param = 3, bin.tol = 0)
  GOF5_threemixedexp <- GOF(O_c = O_c, E_c = E_c, param = 3, bin.tol = 5)
  
  
  alpha_estimate(estimate = C_hat,
                 error = se_hat,
                 estimand = "richness",
                 name = "Two-component Geometric mixture Model",
                 type = "parametric",
                 model = "Two-component Geometric mixture Model",
                 frequentist = TRUE,
                 parametric = TRUE,
                 reasonable = TRUE,
                 interval = C_confint,
                 interval_type = "Approximate: log-normal",
                 GOF0 = pchisq(GOF0_threemixedexp$GOF, GOF0_threemixedexp$df, lower.tail = F) %>% signif(4),
                 GOF5 = pchisq(GOF5_threemixedexp$GOF, GOF5_threemixedexp$df, lower.tail = F) %>% signif(4),
                 AICc = AICc,
                 other = list(t1 = t1,
                              t2 = t2,
                              t3 = t3,
                              t4 = t4,
                              t5 = t5,
                              cutoff = cutoff))
}

# u <- 0.5
# t1 <- 3
# t2 <- 6

## compute the likelihood of the parameters
## if (full == T) return the complete likelihood, with an additional constant term 
three_mixed_exp_lld <- function(ii, ff, t1, t2, t3, u1, u2, full = F) {
  lld <- sum(ff * log(u1 / t1 * (t1 / (1 + t1))^ii + u2 / t2 * (t2 / (1 + t2))^ii +
                 (1 - u1 - u2) / t3 * (t3 / (1 + t3))^ii))
  if (full == T) {
    c_tau <- sum(ff)
    lld <- lld + sum(log(1:c_tau)) - 
      sum(sapply(ff, function(fi) {
        sum(log(1:fi))
      }))
  }
  return(lld)           
}

## Initial values of the EM algorithm
three_mixed_exp_init <- function(ii, ff) {
  init_values <- data.frame(t1 = 0, t2 = 0, t3 = 0, u1 = 0, u2 = 0)
  ss <- length(ii)
  c_tau <- sum(ff)
  for (ss1 in 2:(ss-2)) {
    for (ss2 in (ss1+1):(ss-1)) {
      init_values <- 
        rbind(init_values, 
              c(t1 = sum(ii[1:ss1] * ff[1:ss1]) / sum(ff[1:ss1]) - 1,
                t2 = sum(ii[(ss1 + 1):ss2] * ff[(ss1 + 1):ss2]) / sum(ff[(ss1 + 1):ss2]) - 1,
                t3 = sum(ii[(ss2 + 1):ss] * ff[(ss2 + 1):ss]) / sum(ff[(ss2 + 1):ss]) - 1,
                u1 = sum(ff[1:ss1]) / c_tau,
                u2 = sum(ff[(ss1 + 1):ss2]) / c_tau))
    }
  }
  init_values <- init_values[-1, ]
  init_values$lld <- apply(init_values, 1, function(rr) {
    three_mixed_exp_lld(ii, ff, rr[1], rr[2], rr[3], rr[4], rr[5])
  })

  
  init_values_mle <- init_values[init_values$lld == max(init_values$lld),]
  return(list(t1 = init_values_mle$t1,
              t2 = init_values_mle$t2,
              t3 = init_values_mle$t3,
              u1 = init_values_mle$u1,
              u2 = init_values_mle$u2,
              lld = init_values_mle$lld))
}

three_mixed_exp_EM <- function(ii, ff, MaxIter = 1000, tol = 1e-7) {
  
  c_tau <- sum(ff)
  
  ## initial values for EM algorithm
  inits <- three_mixed_exp_init(ii, ff)
  t1 <- inits$t1
  t2 <- inits$t2
  t3 <- inits$t3
  u1 <- inits$u1
  u2 <- inits$u2

  
  ## begin updating
  k <- 1
  flag <- 0
  
  while (k <= MaxIter) {
    #print(paste("k = ", k, ", t1 = ", t1, ", t2 = ", t2,
    #            ", u1 = ", u1, ", lld = ", lld, sep = ""))
    #print(lld)
    ## update parameters
    z1 <- (u1 / t1 * (t1 / (1 + t1))^ii) / 
      (u1 / t1 * (t1 / (1 + t1))^ii + 
         u2 / t2 * (t2 / (1 + t2))^ii +
         (1- u1 - u2) / t3 * (t3 / (1 + t3))^ii)
    z2 <- (u2 / t2 * (t2 / (1 + t2))^ii) / 
      (u1 / t1 * (t1 / (1 + t1))^ii + 
         u2 / t2 * (t2 / (1 + t2))^ii +
         (1- u1 - u2) / t3 * (t3 / (1 + t3))^ii)
    
    u1_new <- sum(ff * z1) / c_tau 
    u2_new <- sum(ff * z2) / c_tau
    
    t1_new <- sum(ff * ii * z1) / sum(ff * z1) - 1
    t2_new <- sum(ff * ii * z2) / sum(ff * z2) - 1
    t3_new <- sum(ff * ii * (1 - z1 - z2)) / sum(ff * (1 - z1 - z2)) - 1
    
    lld_new <- three_mixed_exp_lld(ii, ff, t1_new, t2_new, t3_new, u1_new, u2_new)
    ## exit iteration if convergent
    if(u1_new / u1 > 1 - tol & u1_new / u1 < 1 + tol &
       u2_new / u2 > 1 - tol & u2_new / u2 < 1 + tol &
       t1_new / t1 > 1 - tol & t1_new / t1 < 1 + tol &
       t2_new / t2 > 1 - tol & t2_new / t2 < 1 + tol &
       t3_new / t3 > 1 - tol & t3_new / t3 < 1 + tol) {
      u1 <- u1_new
      u2 <- u2_new
      t1 <- t1_new
      t2 <- t2_new
      t3 <- t3_new
      lld <- lld_new
      flag <- 1
      break()
    }
    u1 <- u1_new
    u2 <- u2_new
    t1 <- t1_new
    t2 <- t2_new
    t3 <- t3_new
    lld <- lld_new
    k <- k + 1
  }
  return(list(u1 = u1, u2 = u2, t1 = t1, t2 = t2, t3 = t3, iteration = min(k, MaxIter),
              flag = flag))
}

three_mixed_exp_se <- function(cc, t1, t2, t3, t4, t5, MaxIter = 500, tol = 1e-7) {
  a00 <- (1 - t4 / (1 + t1) - t5 / (1 + t2) - (1 - t4 - t5) / (1 + t3)) / 
    (t4 / (1 + t1) + t5 / (1 + t2) + (1 - t4 - t5) / (1 + t3))
  a0 <- - matrix(c(-t4 / (1 + t1) ^ 2, -t5 / (1 + t2) ^ 2, -(1 - t4 - t5) / (1 + t3) ^ 2,
                   1 / (1 + t1) - 1 / (1 + t3), 1 / (1 + t2) - 1 / (1 + t3)), nrow = 5) /
    (t4 / (1 + t1) + t5 / (1 + t2) + (1 - t4 - t5) / (1 + t3))
  A <- 0
  k <- 0
  while(k < MaxIter) {
    p_k <- t4 / (1 + t1) * (t1 / (1 + t1)) ^ k + t5 / (1 + t2) * (t2 / (1 + t2)) ^ k +
      (1 - t4 - t5) / (1 + t3) * (t3 / (1 + t3)) ^ k
    p_k1 <- (t1 / (1 + t1)) ^ k * (k - t1) * t4 / (t1 * (1 + t1) ^ 2)
    p_k2 <- (t2 / (1 + t2)) ^ k * (k - t2) * t5 / (t2 * (1 + t2) ^ 2)
    p_k3 <- (t3 / (1 + t3)) ^ k * (k - t3) * (1 - t4 - t5) / (t3 * (1 + t3) ^ 2)
    
    p_k4 <- 1 / (1 + t1) * (t1 / (1 + t1)) ^ k - 1 / (1 + t3) * (t3 / (1 + t3)) ^ k
    p_k5 <- 1 / (1 + t2) * (t2 / (1 + t2)) ^ k - 1 / (1 + t3) * (t3 / (1 + t3)) ^ k
    p_kj <- matrix(c(p_k1, p_k2, p_k3, p_k4, p_k5), ncol = 1)
    A_k <- p_kj %*% t(p_kj) / p_k
    A <- A + A_k
    if (max(abs(A_k)) < tol) break()
    k <- k + 1
  }
  if (any(svd(A)$d < 1e-15)) {
    warning("stardard error unattainable due to matrix singularity")
    return(NA)
  } else {
    return(as.numeric(sqrt(cc) / sqrt(a00 - t(a0) %*% solve(A) %*% a0)))
  }
}

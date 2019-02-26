#' Two-component geometric mixture model for estimating species richness
#'
#' This function implements the species richness estimation with a two-component mixture
#' exponential-mixed Poisson model (or three-component mixture exponential-mixed geometric model)
#' for estimating species richness
#'
#' @param input_data An input type that can be processed by convert()
#' @param cutoff Maximal frequency count of the data used to estimating the species richness. Default 10.
#'
#' @return An object of class \code{alpha_estimate}
#'
#' @importFrom breakaway convert alpha_estimate
#' @export
#'
#' @examples
#' library(breakaway)
#' data(apples)
#' two_geometric_model(apples, cutoff = 40)

#' @export
two_geometric_model <- function(input_data, cutoff = 10, ...) {
  input_data <- convert(input_data)

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
  em_params <- two_mixed_exp_EM(input_data = included)
  u <- em_params$u
  t1 <- em_params$t1
  t2 <- em_params$t2
  t3 <- u * t2 * (1 + t1) / (t1 + t1 * t2 + t2 * u - t1 * u)

  C_tau <- c_tau * (1 + t1) * (1 + t2) / (t1 * t2 + t2 - t2 * t3 + t1 * t3)
  f0 <- C_tau - c_tau
  C_hat <- C_tau + c_excluded
  ## se
  se_hat <- two_mixed_exp_se(C_tau, t1, t2, t3)

  ## confidence interval
  dd <- ifelse(f0 == 0, 1, exp(1.96 * sqrt(log(1 + se_hat^2/f0))))

  C_confint <- c(c_tau + c_excluded + f0/dd, c_tau + c_excluded + f0 * dd)
  ## AICc
  AICc <- 6 * c_tau / (c_tau - 4) - 2 * sum(log(1:c_tau)) +
    2 * sum(sapply(1:length(ff), function(i) {
      sum(log(1:ff[i]))
    })) -
    2 * sum(ff * log(u / t1 * (t1 / (1 + t1))^ii + (1 - u) / t2 * (t2 / (1 + t2))^ii))


  ## GOF0 and GOF5
  O_c <- c(ff, 0)
  pp_E <- u / t1 * (t1 / (1 + t1))^ii + (1-u) / t2 * (t2 / (1 + t2))^ii
  E_c <- c_tau * c(pp_E, 1 - sum(pp_E))

  GOF0_twomixedexp <- GOF(O_c = O_c, E_c = E_c, param = 3, bin.tol = 0)
  GOF5_twomixedexp <- GOF(O_c = O_c, E_c = E_c, param = 3, bin.tol = 5)


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
                 GOF0 = signif(pchisq(GOF0_twomixedexp$GOF, GOF0_twomixedexp$df, lower.tail = F), 4),
                 GOF5 = signif(pchisq(GOF5_twomixedexp$GOF, GOF5_twomixedexp$df, lower.tail = F), 4),
                 AICc = AICc,
                 other = list(t1 = t1,
                              t2 = t2,
                              t3 = t3,
                              u = u,
                              cutoff = cutoff))
}

# u <- 0.5
# t1 <- 3
# t2 <- 6

## compute the likelihood of the parameters

#' Internal functions for estimating species richness with a two-component geomic mixture model
#'
#' @export
two_mixed_exp_lld <- function(input_data, t1, t2, u) {
  ii <- input_data$index
  ff <- input_data$frequency
  sum(ff * log(u / t1 * (t1 / (1 + t1)) ^ ii +
                 (1 - u) / t2 * (t2 / (1 + t2)) ^ ii))
}

## Initial values of the EM algorithm
#' @export
two_mixed_exp_init <- function(input_data) {
  ii <- input_data$index
  ff <- input_data$frequency
  ss <- 2:(length(ii) - 2)
  c_tau <- sum(ff)
  init_values <- data.frame(u = sapply(ss, function(s) {
                              sum(ff[1:s]) / c_tau
                            }),
                            t1 = sapply(ss, function(s) {
                              sum(ff[1:s] * ii[1:s]) / sum(ff[1:s]) - 1
                            }),
                            t2 = sapply(ss, function(s) {
                              sum(ff[(s + 1):length(ii)] * ii[(s + 1):length(ii)]) /
                                sum(ff[(s + 1):length(ii)]) - 1
                            }))
  init_values$lld <- apply(init_values, 1, function(rr) {
    two_mixed_exp_lld(input_data, u = rr[1], t1 = rr[2], t2 = rr[3])
  })

  init_values_mle <- init_values[init_values$lld == max(init_values$lld),]
  return(list(u = init_values_mle$u,
              t1 = init_values_mle$t1,
              t2 = init_values_mle$t2,
              lld = init_values_mle$lld))
}

#' @export
two_mixed_exp_EM <- function(input_data, MaxIter = 1000, tol = 1e-7) {
  ii <- input_data$index
  ff <- input_data$frequency
  c_tau <- sum(ff)

  ## initial values for EM algorithm
  inits <- two_mixed_exp_init(input_data)
  u <- inits$u
  t1 <- inits$t1
  t2 <- inits$t2
  lld <- inits$lld

  ## begin updating
  k <- 1
  flag <- 0

  while (k <= MaxIter) {
    # print(paste("k = ", k, ", t1 = ", t1, ", t2 = ", t2,
    #             ", u = ", u, ", lld = ", lld, sep = ""))
    #print(lld)
    ## update parameters
    z <- (u / t1 * (t1 / (1 + t1))^ii) /
      (u / t1 * (t1 / (1 + t1))^ii + (1 - u) / t2 * (t2 / (1 + t2))^ii)
    u_new <- sum(ff * z) / c_tau
    t1_new <- sum(ii * ff * z) / sum(ff * z) - 1
    t2_new <- sum(ii * ff * (1 - z)) / sum(ff * (1 - z)) - 1
    lld_new <- two_mixed_exp_lld(input_data, t1_new, t2_new, u_new)

    ## exit iteration if convergent
    if(u_new / u > 1 - tol & u_new / u < 1 + tol &
       t1_new / t1 > 1 - tol & t1_new / t1 < 1 + tol &
       t2_new / t2 > 1 - tol & t2_new / t2 < 1 + tol) {
      u <- u_new
      t1 <- t1_new
      t2 <- t2_new
      lld <- lld_new
      flag <- 1
      break()
    }
    u <- u_new
    t1 <- t1_new
    t2 <- t2_new
    lld <- lld_new
    k <- k + 1
  }
  return(list(u = u, t1 = t1, t2 = t2, iteration = min(k, MaxIter),
              flag = flag))
}

#' @export
two_mixed_exp_se <- function(cc, t1, t2, t3, MaxIter = 500, tol = 1e-7) {
  a00 <- (t2 + t1 * t2 + t1 * t3 - t2 * t3) /
    (1 + t1 - t1 * t3 + t2 * t3)
  a0 <- matrix(c(t3 * (1 + t2) / ((1 + t1) * (1 + t1 + t2 * t3 - t1 * t3)),
                 (1 + t1) * (1 - t3) / ((1 + t2) * (1 + t1 - t1 * t3 + t2 * t3)),
                 (t1 - t2) / (1 + t1 - t1 * t3 + t2 * t3)), nrow = 3)
  A <- 0
  k <- 0
  while(k < MaxIter) {
    p_k <- t3 / (1 + t1) * (t1 / (1 + t1)) ^ k +
      (1 - t3) / (1 + t2) * (t2 / (1 + t2)) ^ k
    p_k1 <- (t1 / (1 + t1)) ^ (k - 1) * (k - t1) * t3 / (1 + t1) ^ 3
    p_k2 <- (t2 / (1 + t2)) ^ (k - 1) * (k - t2) * (1 - t3) / (1 + t2) ^ 3
    p_k3 <- (t1 / (1 + t1)) ^ k / (1 + t1) - (t2 / (1 + t2)) ^ k / (1 + t2)

    p_kj <- matrix(c(p_k1, p_k2, p_k3), ncol = 1)
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

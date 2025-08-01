############################################################################################################
#
# Numerical Optimization Using Probit Link (Adjusted for measurement Error)
# Updated: 19/05/2025
# Author: Vindyani Herath
# Adopted from Ariel Chernofsky (https://github.com/a-chernofsky/mediation_assay_lower_limit/blob/main/noptim.R)
#
############################################################################################################

# Using one assay limit across all ACTG studies
# Y: numeric binary outcome variable (viral rebound suppressed (Y = 1)/ not suppressed (Y = 0))
# M: continuous mediator variable (log10 scale in application): includes measurement error
# C: pre-treatment mediator-outcome common causes
# delta: Indicator of whether M is "above" limit of detection (1 if M > AL, 0 M <= AL)
# lod: limit of detection (assay limit on the log10 scale)
# shift: treatment mediator distribution shift (0.5, 1, 1.5 or 2)
# initial: (assay limit/2) on the log10 scale
# Assay limit is 92 copies per million CD4+ T cells for caRNA and 1 copy per ml for SCA
# SD of measurement error estimate for log10 SCA: 0.47
# SD of measurement error estimate for log10 CA-RNA: 0.29
library(tidyverse)
library(truncnorm)
library(numDeriv)


#######################################################################################################################

noptim_probit_delta <- function(Y, M, C, delta, lod, shift, initial, sigma.u){
  
  # Function to integrate for Y = 1
  p1 <- function(mm, C){
    m.alpha <- alpha0 + alpha1 * C
    y.beta  <- beta0_star + beta1_star * mm + beta2_star * C
    pnorm(y.beta)  * dnorm(mm, m.alpha, sqrt(s_star_sq), log = F) # m.alpha - mean, s_star - standard deviation of M_star,
    # log = F gives the density function of the normal distribution
  }
  
  # Function to integrate for Y = 0
  p0 <- function(mm, C){
    m.alpha <- alpha0 + alpha1 * C
    y.beta <- beta0_star + beta1_star * mm + beta2_star * C
    (1 - pnorm(y.beta)) * dnorm(mm, m.alpha, sqrt(s_star_sq), log = F) # s_star is the standard deviation of M_star
  }
  
  
  # Vectorized integrate functions for Y = 1 and Y = 0 (Vectorize allows R to handle vector inputs.
  # ie, compute the result for multiple values of C)
  vint1 <- Vectorize(function(C) integrate(p1, lower = -Inf, upper = lod, C = C)$value)
  vint0 <- Vectorize(function(C) integrate(p0, lower = -Inf, upper = lod, C = C)$value)
  
  
  # Log-likelihood to be used in the optim function
  # Reference: General-purpose Optimization. R: General-purpose optimization. https://stat.ethz.ch/R-manual/R-patched/RHOME/library/stats/html/optim.html
  log.likelihood <- function(par){
    alpha0 <<- par[1]
    alpha1 <<- par[2]
    s_star_sq <<- par[3]
    beta0_star <<- par[4]
    beta1_star <<- par[5]
    beta2_star <<- par[6]
    
    m.alpha_obs <- alpha0 + alpha1 * C
    y.beta_obs <- beta0_star + beta1_star * M + beta2_star * C
    p <- pnorm(y.beta_obs)
    
    l.obs <- sum(delta * (Y * log(p) + (1 - Y) * log(1 - p) + dnorm(M, m.alpha_obs, sqrt(s_star_sq), log = T)))
    l.cen <- sum( (1 - delta) * (Y * log( vint1(C) ) + (1 - Y) * log( vint0(C) ) ))
    
    return((l.obs + l.cen))
  }
  
  
  cen <- which(delta == 0) #finding obs. where delta = 0
  obs <- which(delta == 1) #finding obs. where delta = 1
  
  #impute missing values with (assay limit/2) on the log10 scale for initial
  # initial was defined as a parameter of the 'noptim_probit_delta' function
  M.init <- M
  M.init[cen] <- initial # initial is log10(assay limit/2)
  
  #initial estimates for parameters
  mlm <- lm( M.init ~ C )
  sigmam_star_sq <- (sigma( mlm ))^2 # extract residual variance
  yglm <- glm(Y ~ M.init + C, family = binomial(link = "probit"))
  
  #initial values for parameters
  pars.init <- c(mlm$coefficients, sigmam_star_sq, yglm$coefficients)
  print("Initial values for parameters:")
  print(pars.init)
  
  # numerical optimization - constrain sigmam_star^2 > 0 (par - initial values for parameters,
  # fn - function to be minimized)
  # By default, 'optim' performs minimization. To maximize a function, use 'control=list("fnscale"=-1)'
  # fnscale - An overall scaling to be applied to the value of fn. If negative, turns the problem into a maximization problem.
  
  fit <- optim(par = pars.init, fn = log.likelihood, control=list("fnscale"=-1), method = "L-BFGS-B",
               lower = c(rep(-Inf, 2), 0.0001, rep(-Inf, 3)), hessian = TRUE)
  
  
  # parameters from optimization
  theta_hat <- fit$par
  
  # Extract estimated parameters
  print("Estimated alpha, sigma.star^2, and beta.star:")
  print( alpha.0_hat <- fit$par[1] )
  print( alpha.1_hat <- fit$par[2] )
  print( sigma_M.star_hat_sq <- fit$par[3] )
  print( beta0.star_hat <- fit$par[4] )
  print( beta1.star_hat <- fit$par[5] )
  print( beta2.star_hat <- fit$par[6] )
  
  # Check if the Hessian matrix is valid
  print("Hessian matrix from optim:")
  print(fit$hessian)
  
  # Information matrix
  I <- fit$hessian
  var_theta <- solve(-I)
  print("var_theta:")
  print(var_theta)
  
  
  # Hessian matrix from NumDeriv
  hes.mat <- hessian(log.likelihood, x = theta_hat, method = "Richardson")
  print("Hessian matrix from NumDeriv:")
  print(hes.mat)
  
  
  print("Variance-covariance matrix (vcov_matrix):")
  print(var_theta)
  se <- sqrt(diag(var_theta))
  print("SE for each parameter:")
  print(se)
  
  # Adjust coefficients for measurement error
  sigma_M_sq <- sigma_M.star_hat_sq - sigma.u^2
  lambda <- sigma_M_sq / sigma_M.star_hat_sq
  beta.1 <- beta1.star_hat / sqrt(abs(lambda^2 - (beta1.star_hat^2 * lambda * sigma.u^2)))
  beta.0 <- (beta0.star_hat * sqrt(1 + (beta.1^2 * lambda * sigma.u^2))) - (beta.1 * (1 - lambda) * alpha.0_hat)
  beta.2 <- (beta2.star_hat * sqrt(1 + (beta.1^2 * lambda * sigma.u^2))) - (beta.1 * (1 - lambda) * alpha.1_hat)
  
  # Calculate shifted and unshifted values
  n <- length(Y)
  cat("sigma_M squared:", sigma_M_sq, "\n")
  cat("beta0:", beta.0, "\n")
  cat("beta1:", beta.1, "\n")
  cat("beta2:", beta.2, "\n")
  
  z <- qnorm(1-0.05/2)
  ci_alpha_0 <- c(alpha.0_hat - z * se[1], alpha.0_hat + z * se[1] )
  ci_alpha_1 <- c(alpha.1_hat - z * se[2], alpha.1_hat + z * se[2])
  ci_sigma_M_sq <- c(sigma_M_sq - z * se[3], sigma_M_sq + z * se[3] )
  ci_beta_0 <- c(beta.0 - z * se[4], beta.0 + z * se[4] )
  ci_beta_1 <- c(beta.1 - z * se[5], beta.1 + z * se[5])
  ci_beta_2 <- c(beta.2 - z * se[6], beta.2 + z * se[6])
  
  cat("estimate & 95% CI for alpha_0:", alpha.0_hat, ci_alpha_0, "\n")
  cat("estimate & 95% CI for alpha_1:", alpha.1_hat, ci_alpha_1, "\n")
  cat("estimate & 95% CI for sigma_M squared:", sigma_M_sq, ci_sigma_M_sq, "\n")
  cat("estimate & 95% CI for beta_0:", beta.0, ci_beta_0, "\n")
  cat("estimate & 95% CI for beta_1:", beta.1, ci_beta_1, "\n")
  cat("estimate & 95% CI for beta_2:", beta.2, ci_beta_2, "\n")
  
  # Step 3 of the estimation procedure
  # numi.0 and numi.1 are the numerator under shift for c = 0 and c = 1 in Theorem 3.2
  # The numerator represtents mu from (A.3)
  numi.0 <- beta.0 + beta.1 * (alpha.0_hat - shift)
  numi.1 <- beta.0 + beta.1 * (alpha.0_hat + alpha.1_hat - shift) + beta.2
  
  # denominator in theorem 3.2
  denom <-  sqrt(1 + (beta.1^2 * sigma_M_sq))
  
  #phati.0 and phati.1 are probit expressions under shift for c = 0 and c = 1
  phati.0 <- pnorm(numi.0 / denom)
  phati.1 <- pnorm(numi.1 / denom)
  
  #num0.0 and num0.1 are the numerators under no shift for c = 0 and c = 1
  num0.0 <- beta.0 + (beta.1 * alpha.0_hat)
  num0.1 <- beta.0 + (beta.1 * (alpha.0_hat + alpha.1_hat)) + beta.2
  
  #phat0.0 and phat0.1 are probit expressions under no shift for c = 0 and c = 1
  phat0.0 <- pnorm( num0.0 / denom )
  phat0.1 <- pnorm( num0.1 / denom )
  
  # ie.1 and ie.0 are estimated indirect effects for c = 0 and c = 1 respectively.
  ie.0 <- phati.0 - phat0.0
  ie.1 <- phati.1 - phat0.1
  
  # indirect effect
  indirect_effect <- (ie.0 * (sum(C == 0) / n)) + (ie.1 * (sum(C == 1) / n))
  
  
  # Check parameter values before passing to the gradient function
  params_g <- c(alpha.0_hat, alpha.1_hat, sigma_M.star_hat_sq, beta0.star_hat, beta1.star_hat, beta2.star_hat)
  print("Parameters passed to gradient_function:")
  print(params_g)
  
  # Calculate asymptotic variance using the delta method
  gradient_function <- function(params_g){
    alpha.0_hat <- params_g[1]
    alpha.1_hat <- params_g[2]
    sigma_M.star_hat_sq <- params_g[3]
    beta0.star_hat <- params_g[4]
    beta1.star_hat <- params_g[5]
    beta2.star_hat <- params_g[6]
    lambda <- 1 - (sigma.u^2 / sigma_M.star_hat_sq)
    
    print("lambda:")
    print(lambda)
    ## derivatives of equation (A.12)
    
    # Denominator of eq. A.12
    denom_A12 <- lambda * sqrt( 1 + (beta1.star_hat^2 * sigma_M.star_hat_sq))
    
    # Numerator of A.12 'under shift' for c = 0 and c = 1 respectively
    num_C0.shift <- (beta0.star_hat * lambda) + beta1.star_hat * ((lambda * alpha.0_hat) - shift)
    num_C1.shift <- (beta0.star_hat * lambda) + beta1.star_hat * ((lambda * alpha.0_hat) + (lambda * alpha.1_hat) - shift) +  (beta2.star_hat * lambda)
    
    # derivative of probit function 'under shift' for c = 0 and c = 1 respectively (eq. A.12).
    phi_C0.shift <- dnorm( num_C0.shift / denom_A12 )
    phi_C1.shift <- dnorm( num_C1.shift / denom_A12 )
    
    # Numerator of A.12 'under no shift' for c = 0 and c = 1 respectively
    num_C0 <- (beta0.star_hat * lambda) + beta1.star_hat * lambda * alpha.0_hat
    num_C1 <- (beta0.star_hat * lambda) + beta1.star_hat * (lambda * alpha.0_hat + lambda * alpha.1_hat) + (beta2.star_hat * lambda)
    
    # derivative of probit function 'under no shift' for c = 0 and c = 1 respectively (eq. A.12).
    phi_C0 <- dnorm( num_C0 / denom_A12 )
    phi_C1 <- dnorm( num_C1 / denom_A12 )
    
    # Probability of C = 0 and C = 1 respectively.
    P_C_0 <- sum(C == 0) / n
    P_C_1 <- sum(C == 1) / n
    
    #derivative w.r.t. alpha.0_hat (eq. A.13)
    denom_A13 <- sqrt(1 + (beta1.star_hat^2 * sigma_M.star_hat_sq))
    grad_alpha0 <- ((beta1.star_hat/denom_A13) * (phi_C0.shift - phi_C0) * P_C_0) +
      ((beta1.star_hat/denom_A13) * (phi_C1.shift - phi_C1) * P_C_1)
    
    #derivative w.r.t. alpha.1_hat (eq. A.14)
    grad_alpha1 <- (beta1.star_hat/denom_A13) * (phi_C1.shift - phi_C1) * P_C_1
    
    #derivative w.r.t. beta0.star_hat (eq. A.15)
    grad_beta0.star_hat <- ((1/denom_A13) * (phi_C0.shift - phi_C0) * P_C_0) +
      ((1/denom_A13) * (phi_C1.shift - phi_C1) * P_C_1)
    
    #derivative w.r.t. beta1.star_hat (eq. A.16)
    denom_A16 <- lambda * (1 + beta1.star_hat^2 * sigma_M.star_hat_sq)^(3/2)
    
    grad_beta1.star_hat <- (phi_C0.shift * (((lambda * alpha.0_hat - shift) - (beta1.star_hat * sigma_M.star_hat_sq * lambda * beta0.star_hat)) / denom_A16) -
                              phi_C0 * (((lambda * alpha.0_hat) - (beta1.star_hat * sigma_M.star_hat_sq * lambda * beta0.star_hat)) / denom_A16)) * P_C_0 +
      (phi_C1.shift *(((lambda * alpha.0_hat + lambda * alpha.1_hat - shift) - (beta1.star_hat * sigma_M.star_hat_sq * lambda * (beta0.star_hat + beta2.star_hat))) / denom_A16) -
         phi_C1 * (((lambda * alpha.0_hat + lambda * alpha.1_hat) - (beta1.star_hat * sigma_M.star_hat_sq * lambda * (beta0.star_hat + beta2.star_hat))) / denom_A16) ) * P_C_1
    
    
    
    #derivative w.r.t. beta2.star_hat (eq. A.17)
    grad_beta2.star_hat <- (1/denom_A13) * (phi_C1.shift - phi_C1) * P_C_1
    
    #derivative w.r.t. sigma_M2.star (eq. A.19)
    lambda_prime <- sigma.u^2 / (sigma_M.star_hat_sq^2)
    denom_A19 <- lambda^2 * (1 + beta1.star_hat^2 * sigma_M.star_hat_sq)
    
    grad_sigma_M2.star_sq <- (phi_C0.shift * ((lambda * denom_A13 * (beta0.star_hat + beta1.star_hat * alpha.0_hat) * lambda_prime) -
                                                ((lambda * (beta0.star_hat + beta1.star_hat * alpha.0_hat) - (beta1.star_hat * shift)) *
                                                   ((denom_A13 * lambda_prime) +
                                                      ((lambda * beta1.star_hat^2)/ (2 * denom_A13) )))) * P_C_0/denom_A19) -
      (phi_C0 * ((lambda * denom_A13 * (beta0.star_hat + beta1.star_hat * alpha.0_hat) * lambda_prime) -
                   ((lambda * (beta0.star_hat + beta1.star_hat * alpha.0_hat)) *
                      ((denom_A13 * lambda_prime) +
                         ((lambda * beta1.star_hat^2)/ (2 * denom_A13) )))) * P_C_0/denom_A19) +
      (phi_C1.shift * ((lambda * denom_A13 * (beta0.star_hat + (beta1.star_hat * (alpha.0_hat + alpha.1_hat)) + beta2.star_hat) * lambda_prime) -
                         ((lambda * (beta0.star_hat + (beta1.star_hat * (alpha.0_hat + alpha.1_hat)) + beta2.star_hat) - (beta1.star_hat * shift)) *
                            ((denom_A13 * lambda_prime) +
                               ((lambda * beta1.star_hat^2)/ (2 * denom_A13) )))) * P_C_1/denom_A19) -
      (phi_C1 * ((lambda * denom_A13 * (beta0.star_hat + (beta1.star_hat * (alpha.0_hat + alpha.1_hat)) + beta2.star_hat) * lambda_prime)-
                   ((lambda * (beta0.star_hat + (beta1.star_hat * (alpha.0_hat + alpha.1_hat)) + beta2.star_hat)) *
                      ((denom_A13 * lambda_prime) +
                         ((lambda * beta1.star_hat^2)/ (2 * denom_A13) )))) * P_C_1/denom_A19)
    
    
    
    return(c(grad_alpha0, grad_alpha1, grad_sigma_M2.star_sq, grad_beta0.star_hat, grad_beta1.star_hat, grad_beta2.star_hat))
  }
  
  
  
  g_gradient <- t(matrix(gradient_function(params_g)))
  # Print the Jacobian result
  print("Jacobian (g_gradient):")
  print(g_gradient)
  
  # Compute the variance of (g(theta)) using the delta method
  variance_g <- ((g_gradient) %*% (var_theta) %*% t(g_gradient))*(88/82)
  
  # Check if variance is computed properly
  print("Variance of g(theta):")
  print(variance_g)
  
  # Standard error
  se_g <- sqrt(variance_g)
  
  # 95% confidence interval
  lower_bound <- indirect_effect - qt(0.975, 82) * se_g
  upper_bound <- indirect_effect + qt(0.975, 82) * se_g
  
  # Return results and confidence interval
  results <- list(indirect_effect = indirect_effect, confidence_interval = c(lower_bound, upper_bound))
  
  cat("Indirect Effect:", results$indirect_effect, "\n")
  cat("95% Confidence Interval:" , results$confidence_interval, "\n")
  return(results$indirect_effect)
}




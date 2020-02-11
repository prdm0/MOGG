# Title: Marshall Olkin Gamma-G (MOGG)
# Author: Pedro Rafael D. Marinho

# Loading libraries. ------------------------------------------------------
library(purrr)

# Baseline functions. -----------------------------------------------------
pdf_w <- function(x, alpha, beta) dweibull(x = x, shape = alpha, scale = beta)
cdf_w <- function(x, alpha, beta) pweibull(q = x, shape = alpha, scale = beta)

# This function creates MOGG functions. ----------------------------------
pdf_mogg <- function(g, G){
  # Using Closures.
  function(x, theta, a, ...){
    if (theta <= 0 || a <= 0) 
      warning("The \"a\" and \"theta\" parameters must be greater than zero.")
    num <- theta * (-log(1 - G(x = x, ...)))^(a - 1) * g(x = x, ...)
    den <- gamma(a) * (theta + (1 - theta) * pgamma(-log(1 - G(x = x, ...)), a, 1L))^2
    num/den
  }
}

# MOGW --------------------------------------------------------------------
rmogw <- function(n = 1L, theta, a, alpha, beta){
  
  cond_c <- function(x, theta, a, alpha, beta){
    num <- pdf_mogw(x, theta, a, alpha, beta)
    den <- dweibull(x, shape = 1, scale = 1)
    -num/den
  }
  
  x_max <- optim(fn = cond_c, method = "BFGS", par = 1, theta = theta, a = a,
                 alpha = alpha, beta = beta)$par
  
  c <- pdf_mogw(x_max, theta, a, alpha, beta)/dweibull(x_max, shape = alpha, scale = beta)

  criterion <- function(y, u){
    num <- pdf_mogw(y, theta, a, alpha, beta)
    den <- dweibull(y, shape = alpha, scale = beta)
    u < num / (c * den)
  }
  
  values <- double(n)
  i <- 1L
  repeat{
    y <- rweibull(n = 1L, shape = alpha, scale = beta)
    u <- runif(n = 1L, min = 0, max = 1)
    
    if (criterion(y, u)) {
      values[i] <- y
      i <- i + 1L
    }
    if(i > n) break
  }
  values
}

# Testing the rmogw Function ----------------------------------------------
# theta = 1.2
# a = 1.5
# alpha = 1.7
# beta = 1.4
# pdf_mogw <- pdf_mogg(g = pdf_w, G = cdf_w)
# sample_data <- rmogw(n = 250L, theta, a, alpha, beta)
# x <- seq(0, 6, length.out = 500L)
# hist(sample_data, probability = TRUE, xlab = "", main = "")
# lines(x, pdf_mogw(x, theta, a, alpha, beta))

# Monte Carlo simulations. ------------------------------------------------
mc <- function(M = 1e3L, n = 100L, method = "BFGS", theta, a, alpha, beta) {
  
  # Log-likelihood function. ------------------------------------------------
  pdf_mogw <- pdf_mogg(g = pdf_w, G = cdf_w)
  log_likelihood <- function(x, par) {
    theta <- par[1L]
    a <- par[2L]
    alpha <- par[3L]
    beta <- par[4L]
    -sum(log(pdf_mogw(x, theta = theta, a = a, alpha = alpha, beta = beta)))
  }
  
  myoptim <-
    function(...)
      tryCatch(
        expr = optim(...),
        error = function(e)
          NA
      )
  
  one_step_mc <- function(i){
    sample_data <- rmogw(n, theta, a, alpha, beta)
    
    result <- myoptim(
      fn = log_likelihood, 
      par = c(1, 1, 1, 1),
      x = sample_data,
      method = method
    )
    
    while(is.na(result) || result$convergence != 0) {
      sample_data <- rmogw(n, theta, a, alpha, beta)
      result <- myoptim(
        fn = log_likelihood, 
        par = c(1, 1, 1, 1),
        method = method,
        x = sample_data
      )
    }
    
    result$par
  }
  
  sapply(X = 1L:M, FUN = one_step_mc) 
}

result_mc <- mc(
  M = 1e3L,
  n = 250L,
  method = "BFGS",
  theta = 1.5,
  a = 1.5,
  alpha = 1.5,
  beta = 1.5
)

mean_est_par <- apply(X = result_mc, MARGIN = 1L, FUN = mean)
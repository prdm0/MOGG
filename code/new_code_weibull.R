# Title: Marshall Olkin Gamma-G (MOGG)
# Author: Pedro Rafael D. Marinho

# Loading libraries. ------------------------------------------------------
library(parallel)
library(tibble)
library(pbmcapply)
library(magrittr)
library(purrr)
library(xtable)
library(glue)
library(fs)

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

pdf_mogw <- pdf_mogg(pdf_w, cdf_w)

# MOGW --------------------------------------------------------------------
rmogw <- function(n = 1L, theta, a, alpha, beta){
  
  cond_c <- function(x, theta, a, alpha, beta){
    num <- pdf_mogw(x, theta, a, alpha, beta)
    den <- pdf_w(x, alpha = 1, beta = 1)
    -num/den
  }
  
  x_max <- optim(fn = cond_c, method = "BFGS", par = 1, theta = theta, a = a,
                 alpha = alpha, beta = beta)$par
  
  c <- pdf_mogw(x_max, theta, a, alpha, beta)/dweibull(x_max, shape = alpha, scale = beta)
  
  criterion <- function(y, u){
    num <- pdf_mogw(y, theta, a, alpha, beta)
    den <- pdf_w(y, alpha, beta)
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
theta = 1.2
a = 1.5
alpha = 1.7
beta = 1.4
pdf_mogw <- pdf_mogg(g = pdf_w, G = cdf_w)
sample_data <- rmogw(n = 250L, theta, a, alpha, beta)
x <- seq(0, 6, length.out = 500L)
hist(sample_data, probability = TRUE, xlab = "", main = "")
lines(x, pdf_mogw(x, theta, a, alpha, beta))

# Monte Carlo simulations. ------------------------------------------------
mc <- function(n = 250L, M = 1e3L, par_true, method = "BFGS") {
  
  theta <- par_true[1L]
  a <- par_true[2L]
  alpha <- par_true[3L]
  beta <- par_true[4L]
  
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
  
  result_vector <-
    unlist(pbmcapply::pbmclapply(
      X = 1L:M,
      FUN = one_step_mc,
      mc.cores = parallel::detectCores()
    ))
  
  
  result <-
    tibble::as_tibble(matrix(result_vector, byrow = TRUE, ncol = 4L))
  
  names(result) <- c("theta", "a", "alpha", "beta")
  
  result
}

bias_function <- function(x, par_true){
  x - par_true
}

mse_function <- function(x, par_true) {
  (x - par_true) ^ 2
}

simulate <- function(n) {
  # True parameters (theta, a, alpha and beta) -------------------------------------
  true_parameters <- c(1, 1, 1, 1)
  
  set.seed(1L, kind = "L'Ecuyer-CMRG")
  t0 <- Sys.time()
  result_mc <- mc(n = n, M = 1e4L, par_true = true_parameters, method = "BFGS")
  total_time <- difftime(Sys.time(), t0, units = 'mins')[[1L]] 
  
  mc.reset.stream()
  
  # Average Bias of Estimators ----------------------------------------------
  eval(parse(text = glue("bias_{n} <- apply(X = result_mc, MARGIN = 1L, FUN = bias_function, par_true = true_parameters) %>% 
        apply(MARGIN = 1L, FUN = mean)"))) 
  eval(parse(text = glue("save(file = \"rdata/weibull/bias_{n}.RData\", bias_{n})")))
  
  # Mean Square Error -------------------------------------------------------
  eval(parse(text = glue("mse_{n} <- apply(X = result_mc, MARGIN = 1L, FUN = mse_function, par_true = true_parameters) %>% 
        apply(MARGIN = 1L, FUN = mean)"))) 
  eval(parse(text = glue("save(file = \"rdata/weibull/mse_{n}.RData\", mse_{n})")))
  
  # Total Time --------------------------------------------------------------
  eval(parse(text = glue("time_{n} <- total_time")))
  eval(parse(text = glue("save(file = \"rdata/weibull/time_{n}.RData\", time_{n})")))
  
  # Result MC
  eval(parse(text = glue("result_{n} <- result_mc")))
  eval(parse(text = glue("save(file = \"rdata/weibull/result_{n}.RData\", result_{n})")))
  
}

n <- c(10, 20, 60, 100, 200, 400, 600, 1000, 2000, 5000, 10000, 20000, 30000, 50000)

walk(.x = n, .f = simulate)

# Reload saved datasets
# path_data <- glue("rdata/weibull/{dir('rdata/weibull')}")
# results <- purrr::map(.x = path_data, .f = \(x) mget(load(x)))

# bias <- rbind(bias_10, bias_20, bias_60, bias_100, bias_200, bias_400, bias_600, bias_1000, bias_5000, bias_10000,
#                 bias_20000, bias_30000, bias_50000)

# times <- c(time_10, time_20, time_60, time_100, time_200, time_400, time_600, time_1000, time_5000, time_10000,
#   time_20000, time_30000, time_50000)

# Load bias
load(file = "rdata/weibull/bias_10.RData")
load(file = "rdata/weibull/bias_20.RData")
load(file = "rdata/weibull/bias_60.RData")
load(file = "rdata/weibull/bias_100.RData")
load(file = "rdata/weibull/bias_200.RData")
load(file = "rdata/weibull/bias_400.RData")
load(file = "rdata/weibull/bias_600.RData")
load(file = "rdata/weibull/bias_1000.RData")
load(file = "rdata/weibull/bias_2000.RData")
load(file = "rdata/weibull/bias_5000.RData")
load(file = "rdata/weibull/bias_10000.RData")
load(file = "rdata/weibull/bias_20000.RData")
load(file = "rdata/weibull/bias_30000.RData")
load(file = "rdata/weibull/bias_50000.RData")

# Load times
load(file = "rdata/weibull/time_10.RData")
load(file = "rdata/weibull/time_20.RData")
load(file = "rdata/weibull/time_60.RData")
load(file = "rdata/weibull/time_100.RData")
load(file = "rdata/weibull/time_200.RData")
load(file = "rdata/weibull/time_400.RData")
load(file = "rdata/weibull/time_600.RData")
load(file = "rdata/weibull/time_1000.RData")
load(file = "rdata/weibull/time_2000.RData")
load(file = "rdata/weibull/time_5000.RData")
load(file = "rdata/weibull/time_10000.RData")
load(file = "rdata/weibull/time_20000.RData")
load(file = "rdata/weibull/time_30000.RData")
load(file = "rdata/weibull/time_50000.RData")

tabela <-
  rbind(
    bias_10,
    bias_20,
    bias_60,
    bias_100,
    bias_200,
    bias_400,
    bias_600,
    bias_1000,
    bias_2000,
    bias_5000,
    bias_10000,
    bias_20000,
    bias_30000,
    bias_50000
  )
tabela <- cbind(
  tabela,
  c(
    time_10,
    time_20,
    time_60,
    time_100,
    time_200,
    time_400,
    time_600,
    time_1000,
    time_2000,
    time_5000,
    time_10000,
    time_20000,
    time_30000,
    time_50000
  )
)
colnames(tabela) <- c("theta", "a", "alpha", "beta", "Times (mins)") 
rownames(tabela) <- NULL
tabela <- tabela |> round(digits = 4L) |> as_tibble()

tabela <- as.vector(unlist(tabela)) |> matrix(byrow = F, nrow = length(n))
colnames(tabela) <- c("theta", "a", "alpha", "beta", "Times (mins)") 
tabela <- cbind(n = n, tabela) |> round(digits = 4L) |> as_tibble()
latex <- print.xtable(xtable(tabela,  caption = "Mean bias of EMV obtained by the BFGS method in 10,000 Monte Carlo repetitions.",
                             digits = 4L), print.results = FALSE)

writeLines(
  c(
    "\\documentclass[12pt]{article}",
    "\\begin{document}",
    "\\thispagestyle{empty}",
    latex,
    "\\end{document}"
  ),
  "simulations/mc_simulation_weibull.tex"
)

tools::texi2pdf("simulations/mc_simulation_weibull.tex", clean = TRUE)

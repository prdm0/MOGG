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
library(kableExtra)

# Baseline functions. -----------------------------------------------------
pdf_dagum <- function(x, alpha, beta, p) 
  alpha * p / x * (x / beta) ^ (alpha * p) / ((x / beta) ^ alpha + 1) ^ (p + 1)
# integrate(f = pdf_dagum, lower = 0, upper = Inf, alpha = 1.2, beta = 1.6, p = 2.2)

cdf_dagum <- function(x, alpha, beta, p) 
  (1 + (x / beta) ^ (-alpha)) ^ (-p)
# cdf_dagum(x = Inf, alpha = 1, beta = 4, p = 1)

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

rdagum <- function(n = 1L, alpha, beta, p) {
  beta * (runif(n = n, min = 0, max = 1) ^ (-1 / p) - 1) ^ (-1 / alpha)
}

pdf_mogdagum <- pdf_mogg(g = pdf_dagum, G = cdf_dagum)

# MOG-Dagum --------------------------------------------------------------------
rmogdagum <- function(n = 1L, theta, a, alpha, beta, p){
  
  cond_c <- function(x, theta, a, alpha, beta, p){
    num <- pdf_mogdagum(x, theta, a, alpha, beta, p)
    den <- pdf_dagum(x, alpha = alpha, beta = beta , p = p)
    -num / den
  }
  
  x_max <- optim(fn = cond_c, method = "BFGS", par = 1, theta = theta, a = a,
                 alpha = alpha, beta = beta, p = p)$par
  
  c <- pdf_mogdagum(x_max, theta, a, alpha, beta, p)/pdf_dagum(x_max, alpha = alpha, beta = beta, p = p)

  criterion <- function(y, u){
    num <- pdf_mogdagum(y, theta, a, alpha, beta, p)
    den <- pdf_dagum(y, alpha = alpha, beta = beta, p = p)
    u < num / (c * den)
  }
  
  values <- double(n)
  i <- 1L
  repeat{
    y <- rdagum(n = 1L, alpha = alpha, beta = beta, p = p)
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
theta = 5
a = 1.4
alpha = 5
beta = 1
p = 1
pdf_mogdagum <- pdf_mogg(g = pdf_dagum, G = cdf_dagum)
sample_data <- rmogdagum(n = 250L, theta, a, alpha, beta, p)
x <- seq(0, max(sample_data), length.out = 500L)
hist(sample_data, probability = TRUE, xlab = "", main = "")
lines(x, pdf_mogdagum(x, theta, a, alpha, beta, p))

# Monte Carlo simulations. ------------------------------------------------
mc <- function(n = 250L, M = 1e3L, par_true, method = "BFGS") {
  
  theta <- par_true[1L]
  a <- par_true[2L]
  alpha <- par_true[3L]
  beta <- par_true[4L]
  p <- par_true[5L]
  
  # Log-likelihood function. ------------------------------------------------
  pdf_mogw <- pdf_mogg(g = pdf_dagum, G = cdf_dagum)
  log_likelihood <- function(x, par) {
    theta <- par[1L]
    a <- par[2L]
    alpha <- par[3L]
    beta <- par[4L]
    p <- par[5L]
    
    -sum(log(pdf_mogdagum(x, theta = theta, a = a, alpha = alpha, beta = beta, p = p)))
  }
  
  myoptim <-
    function(...)
      tryCatch(
        expr = optim(...),
        error = function(e)
          NA
      )
  
  one_step_mc <- function(i){
    sample_data <- rmogdagum(n, theta, a, alpha, beta, p)
    
    result <- myoptim(
      fn = log_likelihood, 
      par = c(1, 1, 1, 1, 1),
      x = sample_data,
      method = method
    )
    
    while(is.na(result) || result$convergence != 0) {
      sample_data <- rmogdagum(n, theta, a, alpha, beta, p)
      result <- myoptim(
        fn = log_likelihood, 
        par = c(1, 1, 1, 1, 1),
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
    tibble::as_tibble(matrix(result_vector, byrow = TRUE, ncol = 5L))
  
  names(result) <- c("theta", "a", "alpha", "beta", "p")
  
  result
}

bias_function <- function(x, par_true){
  round(x - par_true, digits = 3L)
}

mse_function <- function(x, par_true) {
  round((x - par_true) ^ 2, digits = 3L)
}

simulate <- function(n) {
  # True parameters (theta, a, alpha, beta and p) -------------------------------------
  true_parameters <- c(1, 1, 1, 1, 1)
  
  set.seed(1L, kind = "L'Ecuyer-CMRG")
  t0 <- Sys.time()
  result_mc <- mc(n = n, M = 1e4L, par_true = true_parameters, method = "BFGS")
  total_time <- difftime(Sys.time(), t0, units = 'mins')[[1L]] 
    
  mc.reset.stream()
  
  # Average Bias of Estimators ----------------------------------------------
  eval(parse(text = glue("bias_{n} <- apply(X = result_mc, MARGIN = 1L, FUN = bias_function, par_true = true_parameters) %>% 
        apply(MARGIN = 1L, FUN = mean)"))) 
  eval(parse(text = glue("save(file = \"rdata/dagum/bias_{n}.RData\", bias_{n})")))
  
  # Mean Square Error -------------------------------------------------------
  eval(parse(text = glue("mse_{n} <- apply(X = result_mc, MARGIN = 1L, FUN = mse_function, par_true = true_parameters) %>% 
        apply(MARGIN = 1L, FUN = mean)"))) 
  eval(parse(text = glue("save(file = \"rdata/dagum/mse_{n}.RData\", mse_{n})")))
  
  # Total Time --------------------------------------------------------------
  eval(parse(text = glue("time_{n} <- total_time")))
  eval(parse(text = glue("save(file = \"rdata/dagum/time_{n}.RData\", time_{n})")))
  
  # Result MC
  eval(parse(text = glue("result_{n} <- result_mc")))
  eval(parse(text = glue("save(file = \"rdata/dagum/result_{n}.RData\", result_{n})")))
  
}
n <- c(10, 20, 60, 100, 200, 400, 600, 1000, 2000, 5000, 10000, 20000, 30000, 50000)
walk(.x = n, .f = simulate)

# Load bias
load(file = "rdata/dagum/bias_10.RData")
load(file = "rdata/dagum/bias_20.RData")
load(file = "rdata/dagum/bias_60.RData")
load(file = "rdata/dagum/bias_100.RData")
load(file = "rdata/dagum/bias_200.RData")
load(file = "rdata/dagum/bias_400.RData")
load(file = "rdata/dagum/bias_600.RData")
load(file = "rdata/dagum/bias_1000.RData")
load(file = "rdata/dagum/bias_2000.RData")
load(file = "rdata/dagum/bias_5000.RData")
load(file = "rdata/dagum/bias_10000.RData")
load(file = "rdata/dagum/bias_20000.RData")
load(file = "rdata/dagum/bias_30000.RData")
load(file = "rdata/dagum/bias_50000.RData")

# Load times
load(file = "rdata/dagum/time_10.RData")
load(file = "rdata/dagum/time_20.RData")
load(file = "rdata/dagum/time_60.RData")
load(file = "rdata/dagum/time_100.RData")
load(file = "rdata/dagum/time_200.RData")
load(file = "rdata/dagum/time_400.RData")
load(file = "rdata/dagum/time_600.RData")
load(file = "rdata/dagum/time_1000.RData")
load(file = "rdata/dagum/time_2000.RData")
load(file = "rdata/dagum/time_5000.RData")
load(file = "rdata/dagum/time_10000.RData")
load(file = "rdata/dagum/time_20000.RData")
load(file = "rdata/dagum/time_30000.RData")
load(file = "rdata/dagum/time_50000.RData")

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
colnames(tabela) <- c("theta", "a", "alpha", "beta", "p", "Times (mins)") 
rownames(tabela) <- NULL
tabela <- tabela |> round(digits = 4L) |> as_tibble()

tabela <- as.vector(unlist(tabela)) |> matrix(byrow = F, nrow = length(n))
colnames(tabela) <- c("theta", "a", "alpha", "beta", "p", "Times (mins)") 
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
  "simulations/mc_simulation_dagum.tex"
)

tools::texi2pdf("simulations/mc_simulation_dagum.tex", clean = TRUE)

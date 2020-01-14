############################################################
###############   MO-Ga-Weibull DISTRIBUTION     ###########
############################################################
#PURPOSE: Marshall Olkin Gamma-Weibull (MO-Ga-Weibull) model 
#AUTHOR: Tatiane F. Ribeiro
#DATE: December 26th, 2019.
############################################################
#Clear the memory
rm(list = ls())
#Seed
set.seed(19)

#Inicializations
R = 10 #number of replications
vn<- c(100,150, 200,250) #sample sizes

#Parameters 
a = 1
theta = 4
lambda = 2
gama = 2

#Weibull
dwei <- function(x){
  lambda*gama*(lambda*x)^(gama-1)*exp(-(lambda*x)^gama)
}
pwei <- function(x){
  1-exp(-((lambda*x)^gama))
}
#Test density
integrate(dwei,0,Inf)

qwei <- function(u,par1,par2){
  1/lambda*(-log(1-u))^(1/gama)
}
#Tests quantile function
aux=pwei(c(3,.4,5,.7,6))
qwei(aux)

#cdf MO-Ga-Weibull
cdf_MOGa_Wei <- function(x){
  pgamma(-log(1-pwei(x)),a)/((1-theta)*
                               pgamma(-log(1-pwei(x)),a)+theta)
}

#pdf MO-Ga-Weibull
pdf_MOGa_Wei <- function(x){
  g = lambda*gama*(lambda*x)^(gama-1)*exp(-(lambda*x)^gama)
  G =   1-exp(-((lambda*x)^gama))
  
  num = theta*g*(-log(1-G))^(a-1)
  den = ((1-theta)*pgamma(-log(1-G),a,1)+theta)^2
  
  num/(gamma(a)*den)
}
#integrate(pdf_MOGa_Wei,0,Inf)

#Quantile function MO-Ga-Weibull
qf_MOGa_Wei <- function(u){
  quant = 1-exp(-qgamma(u*theta/(-u*(1-theta)+1),a))
  res = qwei(quant,lambda,gama)    
  return(res)
}

#Inversion method
r_MOGa_Wei <- function(size,par1,par2,par3,par4){
  u = runif(n, min = 0, max = 1)
  quant = 1-exp(-qgamma(u*par2/(-u*(1-par2)+1),par1))
  z = qwei(quant,par3,par4)  
  return(z)
}
#Log-likelihood function
l_MOGa_Wei <- function(par){
  a <- par[1]
  theta <- par[2]
  lambda <- par[3]
  gama <- par[4]
  
  n*log(theta)+n*log(gama)+n*a*gama*log(lambda)+(a*gama-1)*
    sum(log(x))-lambda^gama*sum(x^gama)-
    n*log(gamma(a))-2*sum(log(theta+(1-theta)*pgamma((lambda*x)^gama,a)))
}

#Monte Carlo Simulation
Mat_mle_a = Mat_mle_theta = matrix(NA, length(vn), 4)
Mat_mle_lambda = Mat_mle_gama = matrix(NA, length(vn), 4)

chute <- c(a,theta,lambda,gama) #true value of the parameters (just test)
#Loop n
i = 1 #index results matrix 
for(n in vn){
  mle_a = mle_theta = mle_lambda = mle_gama = rep(NA,R)
  falhas = 0
  
  #Loop MC
  k = 0
  while(k < R){
    x = r_MOGa_Wei(n,a,theta,lambda, gama)
    max_log_lik <- optim(chute,l_MOGa_Wei,
                         method = "S", 
                         control = list(fnscale=-1))
    
    if(max_log_lik$convergence == 0){
      k = k+1
      mle_a[k] = max_log_lik$par[1]
      mle_theta[k] = max_log_lik$par[2]
      mle_lambda[k] = max_log_lik$par[3]
      mle_gama[k] = max_log_lik$par[4]
    }
    
    else{
      falhas = falhas + 1
    }
  }
  
  #Mean LME's
  m_mle_a = mean(mle_a)
  m_mle_theta = mean(mle_theta)
  m_mle_lambda = mean(mle_lambda)
  m_mle_gama = mean(mle_gama)
  
  #Bias LME's
  b_mle_a <- m_mle_a - a
  b_mle_theta <- m_mle_theta - theta
  b_mle_lambda <- m_mle_lambda - lambda
  b_mle_gama <- m_mle_gama - gama
  
  #Root mean squared error (RMSE)
  rmse_mle_a = sqrt(sd(mle_a)^2+b_mle_a^2)
  rmse_mle_theta = sqrt(sd(mle_theta)^2+b_mle_theta^2)
  rmse_mle_lambda = sqrt(sd(mle_lambda)^2+b_mle_lambda^2)
  rmse_mle_gama = sqrt(sd(mle_gama)^2+b_mle_gama^2)
  
  #Results
  colnames(Mat_mle_a) <- c("n","m_mle_a","b_mle_a","rmse_mle_a")
  Mat_mle_a[i,] <- c(n,m_mle_a,b_mle_a,rmse_mle_a)
  
  colnames(Mat_mle_theta) <- c("n","m_mle_theta","b_mle_theta","rmse_mle_theta")
  Mat_mle_theta[i,] <- c(n,m_mle_theta,b_mle_theta,rmse_mle_theta)
  
  colnames(Mat_mle_lambda) <- c("n","m_mle_lambda","b_mle_lambda","rmse_mle_lambda")
  Mat_mle_lambda[i,] <- c(n,m_mle_lambda,b_mle_lambda,rmse_mle_lambda)
  
  colnames(Mat_mle_gama) <- c("n","m_mle_gama","b_mle_gama","rmse_mle_gama")
  Mat_mle_gama[i,] <- c(n,m_mle_gama,b_mle_gama,rmse_mle_gama)
  
  i = i+1
  print(i)
}

Mat_mle_a
Mat_mle_theta
Mat_mle_lambda
Mat_mle_gama


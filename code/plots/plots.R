  library(ggplot2)
  library(cowplot)
  
  # Baseline functions. -----------------------------------------------------
  pdf_weibull <- function(x, alpha, beta)
   dweibull(x, shape = alpha, scale = beta)
  
  cdf_weibull <- function(x, alpha, beta)
    pweibull(x, shape = alpha, scale = beta)
  
  # This function creates MOGW functions ----------------------------------
  pdf_mogw <- function(g, G) {
    # Using Closures.
    function(x, theta, a, ...) {
      if (theta <= 0 || a <= 0)
        warning("The \"a\" and \"theta\" parameters must be greater than zero.")
      num <-
        theta * (-log(1 - G(x = x, ...))) ^ (a - 1) * g(x = x, ...)
      den <-
        gamma(a) * (theta + (1 - theta) * pgamma(-log(1 - G(x = x, ...)), a, 1L)) ^
        2
      num / den
    }
  }
  
  pdf_mogweibull <- pdf_mogw(g = pdf_weibull, G = cdf_weibull)
  
  hazard <- function(x, theta, a, alpha, beta) {
    survival_f <- function(x, theta, a, alpha, beta)
      1 - integrate(
        f = pdf_mogweibull,
        lower = .Machine$double.eps,
        upper = x,
        theta,
        a,
        alpha,
        beta
      )$value
    
    vectorize_survival_f <-
      Vectorize(FUN = survival_f, vectorize.args = "x")
    
    pdf_mogweibull(x, theta, a, alpha, beta) / vectorize_survival_f(x, theta, a, alpha, beta)
  }
  
  #  Plots Probability densiy function  -------------------------------------
  cores <- c("Legenda 1" = "#FF6B58",
             "Legenda 2" = "#FF4373",
             "Legenda 3" = "#9239F6",
             "Legenda 4" = "#E3D26D")
  
  base <- ggplot(data.frame(x = c(0, 10)), aes(x)) +
    ggtitle(
      "Probability density function",
      subtitle = "With G following the Weibull distribution "
    ) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(face = "bold")
    )
  
  p1 <- 
    base +
    stat_function(
      fun = pdf_mogweibull,
      args = list(
        theta = 2.2,
        a = 3.2,
        alpha = 2,
        beta = 3
      ),
      aes(color = "Legenda 1"),
      size = 1.5
    ) +
    
    stat_function(
      fun = pdf_mogweibull,
      args = list(
        theta = 3,
        a = 1.5,
        alpha = 0.7,
        beta = 0.3
      ),
      aes(color = "Legenda 2"),
      size = 1.5
    ) +
    
    stat_function(
      fun = pdf_mogweibull,
      args = list(
        theta = 0.2,
        a = 4,
        alpha = 1.4,
        beta = 2.1
      ),
      aes(color = "Legenda 3"),
      size = 1.5
    ) +
    
    stat_function(
      fun = pdf_mogweibull,
      args = list(
        theta = 5.1,
        a = 4.1,
        alpha = 2.5,
        beta = 2.7
      ),
      aes(color = "Legenda 4"),
      size = 1.5
    ) +
    
    labs(x = "x",
         y = expression("f(x)"),
         color = expression(paste("MO-", Gamma, "-Weibull")),
         element_text(face = "bold"))  +
    
    scale_color_manual(values = cores, labels =
                         expression(
                           paste(theta, " = 2.2, ", "a = 3.2, ", alpha, " = 2.0, ", beta, " = 3.0"),
                           paste(theta, " = 3.0, ", "a = 1.5, ", alpha, " = 0.7, ", beta, " = 0.3"),
                           paste(theta, " = 0.2, ", "a = 4.0, ", alpha, " = 1.4, ", beta, " = 2.1"),
                           paste(theta, " = 5.1, ", "a = 4.1, ", alpha, " = 2.5, ", beta, " = 2.7")
                         )) +
    theme(
      legend.title = element_text(face = "bold"),
      legend.background = element_rect(color = "black"),
      legend.position = c(0.7, 0.84),
      text = element_text(size = 18)
    )
  
  
  # Plot hazard function ----------------------------------------------------
  
  cores <- c("Legenda 1" = "#FF6B58",
             "Legenda 2" = "#FF4373",
             "Legenda 3" = "#9239F6",
             "Legenda 4" = "#00B19D") 
  
  base <- ggplot(data.frame(x = c(0, 25)), aes(x)) +
    ggtitle(
      "Hazard Function",
      subtitle = "With G following the Weibull distribution "
    ) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(face = "bold")
    )
  
  p2 <- 
    base +
    stat_function(
      fun = hazard,
      args = list(
        theta = 1,
        a = 2.5,
        alpha = 0.4,
        beta = 0.3
      ),
      aes(color = "Legenda 1"),
      size = 1.5
    ) +
    
    stat_function(
      fun = hazard,
      args = list(
        theta = 4.5,
        a = 1.5,
        alpha = 0.65,
        beta = 0.43
      ),
      aes(color = "Legenda 2"),
      size = 1.5
    ) +
    
    stat_function(
      fun = hazard,
      args = list(
        theta = 2,
        a = 1,
        alpha = 1.5,
        beta = 5
      ),
      aes(color = "Legenda 3"),
      size = 1.5
    ) +
    
    stat_function(
      fun = hazard,
      args = list(
        theta = 1,
        a = 1,
        alpha = 1,
        beta = 7
      ),
      aes(color = "Legenda 4"),
      size = 1.5
    ) + 
    labs(x = "x",
         y = expression("h(x)"),
         color = expression(paste("MO-", Gamma, "-Weibull")),
         element_text(face = "bold"))  +
    
    scale_color_manual(values = cores, labels =
                         expression(
                           paste(theta, " = 1.0, ", "a = 2.5, ", alpha, " = 0.40, ", beta, " = 0.30"),
                           paste(theta, " = 4.5, ", "a = 1.5, ", alpha, " = 0.65, ", beta, " = 0.43"),
                           paste(theta, " = 2.0, ", "a = 1.0, ", alpha, " = 1.50, ", beta, " = 5.00"),
                           paste(theta, " = 1.0, ", "a = 1.0, ", alpha, " = 1.00, ", beta, " = 7.00")
                         )) +
    theme(
      legend.title = element_text(face = "bold"),
      legend.background = element_rect(color = "black"),
      legend.position = c(0.65, 0.39),
      text = element_text(size = 18)
    )
  
  pdf(file="pdf-hazard.pdf",width=15,height=7, paper="special",
      family="Bookman",pointsize=14)
  plot_grid(p1, p2, labels = c("A", "B"))
  dev.off()
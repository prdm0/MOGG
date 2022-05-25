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

type_lines <- c("solid",
                "dotted",
                "dotdash",
                "twodash")

base <- ggplot(data.frame(x = c(0, 10)), aes(x)) +
  # ggtitle(
  #   "Probability density function",
  #   subtitle = "With G following the Weibull distribution "
  # ) +
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
    size = 1.5,
    linetype = "solid"
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
    size = 1.5,
    linetype = "dotted"
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
    size = 1.5,
    linetype = "dotdash"
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
    size = 1.5,
    linetype = "twodash"
  )  +
  
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
                       ), 
                     guide = guide_legend(override.aes = list(
                         linetype = type_lines
                         ))) +
  theme(
    legend.title = element_text(face = "bold"),
    legend.background = element_rect(color = "black"),
    legend.position = c(0.7, 0.84),
    text = element_text(size = 15),
    legend.key.width = unit(3,"cm"),
    legend.text = element_text(size=10),
    legend.spacing.x = unit(0.7, 'cm')
  )



# Plot hazard function ----------------------------------------------------

cores <- c("Legenda 1" = "#FF6B58",
           "Legenda 2" = "#FF4373",
           "Legenda 3" = "#9239F6",
           "Legenda 4" = "#00B19D")

base <- ggplot(data.frame(x = c(0, 25)), aes(x)) +
  # ggtitle(
  #   "Hazard Function",
  #   subtitle = "With G following the Weibull distribution "
  # ) +
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
    size = 1.5,
    linetype = "solid"
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
    size = 1.5,
    linetype = "dotted"
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
    size = 1.5,
    linetype = "dotdash"
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
    size = 1.5,
    linetype = "twodash"
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
                       ), guide = guide_legend(override.aes = list(
                         linetype = type_lines
                       ))) +
  theme(
    legend.title = element_text(face = "bold"),
    legend.background = element_rect(color = "black"),
    legend.position = c(0.65, 0.39),
    text = element_text(size = 15),
    legend.key.width = unit(3,"cm"),
    legend.text = element_text(size=10),
    legend.spacing.x = unit(0.7, 'cm')
  )

pdf(file="pdf-hazard.pdf",width=15,height=7, paper="special",
    family="Bookman",pointsize=12)
  plot_grid(p1, p2, labels = c("A", "B"))
dev.off()

first_application <-
  c(
    0.4210000,
    0.4382105,
    0.5721316,
    0.1955526,
    0.2758158,
    0.3565263,
    0.6120000,
    0.4730526,
    0.3726579,
    0.2322368,
    0.3732632,
    0.5158947,
    0.3612632,
    0.5563421,
    0.4591053,
    0.2533947,
    0.3685526,
    0.2410526,
    0.4955526,
    0.2010789,
    0.3653158,
    0.2544737,
    0.5685263,
    0.2859474,
    0.7626316,
    0.2863684,
    0.4884474,
    0.3060263,
    0.4754474,
    0.3826053,
    0.5050526,
    0.6820526,
    0.7587632,
    0.4176053,
    0.3923684,
    0.2513158,
    0.6070000,
    0.3881842
  )

pdf(file="TTT_app_1.pdf",width=7,height=7, paper="special",
    family="Bookman",pointsize=16)
AdequacyModel::TTT(first_application, lwd = 4)
dev.off()

second_application <- c(
  150,
  120,
  120,
  180,
  138,
  115,
  130,
  150,
  200,
  120,
  190,
  90,
  130,
  120,
  200,
  140,
  110,
  134,
  160,
  140,
  105,
  126,
  129,
  120,
  100,
  130,
  118,
  144,
  180,
  138,
  110,
  140,
  120,
  118,
  110,
  110,
  130,
  140,
  130,
  165,
  180,
  130,
  140,
  112,
  130,
  158,
  112,
  150,
  140,
  142,
  110,
  140,
  130,
  132,
  140,
  140,
  122,
  128,
  90,
  118,
  120,
  110,
  122,
  200,
  110,
  140,
  150,
  120,
  150,
  120,
  164,
  122,
  112,
  130,
  140,
  102,
  122,
  130,
  102,
  130,
  122,
  200,
  140,
  180,
  124,
  110,
  124,
  90,
  120,
  159,
  142,
  140,
  118,
  122,
  108,
  170,
  120,
  140,
  100,
  118,
  110,
  114,
  150,
  160,
  140,
  190,
  118,
  120,
  150,
  120,
  200,
  150,
  168,
  110,
  142,
  150,
  160,
  142,
  160,
  150,
  110,
  128,
  122,
  150,
  140,
  122,
  120,
  130,
  100,
  130,
  150,
  130,
  100,
  120,
  105,
  100,
  150,
  196,
  130,
  110,
  140,
  122,
  110,
  164,
  120,
  120,
  150,
  160,
  150,
  135,
  124,
  110,
  100,
  95,
  130,
  120,
  108,
  118,
  170,
  105,
  120,
  95,
  95,
  120,
  140,
  142,
  160,
  110,
  190,
  180,
  130,
  130,
  120,
  204,
  150,
  150,
  120,
  122,
  120,
  130,
  140,
  148,
  118,
  126,
  136,
  140,
  130,
  102,
  110,
  110,
  130,
  126,
  142,
  140,
  128,
  130,
  124,
  162,
  130,
  130,
  110,
  80,
  166,
  140,
  160,
  160,
  140,
  98,
  138,
  120,
  112,
  112,
  134,
  140,
  115,
  140,
  98,
  115,
  120,
  80,
  160,
  126,
  110,
  130,
  104,
  236,
  118,
  120,
  140,
  120,
  98,
  164,
  150,
  110,
  120,
  130,
  170,
  180,
  110,
  120,
  130,
  118,
  130,
  190,
  158,
  90,
  99,
  210,
  180,
  140,
  184,
  105,
  120,
  150,
  140,
  130,
  160,
  118,
  210,
  100,
  170,
  150,
  130,
  170,
  150,
  120,
  134,
  90,
  125,
  170,
  140,
  150,
  110,
  105,
  140,
  120,
  100,
  124,
  112,
  160,
  140,
  118,
  190,
  110,
  118,
  160,
  150,
  124,
  128,
  150,
  120,
  125,
  118,
  132,
  110,
  143,
  170,
  98,
  124,
  180,
  178,
  110,
  98,
  159,
  110,
  140,
  130,
  122,
  110,
  98,
  180,
  90,
  118,
  165,
  138,
  138,
  170,
  106,
  170,
  140,
  90,
  118,
  110,
  102,
  102,
  180,
  100,
  110,
  162,
  140,
  110,
  98,
  140,
  140,
  110,
  170,
  112,
  90,
  102,
  106,
  124,
  110,
  180,
  138,
  90,
  150,
  126,
  110,
  130,
  150,
  145,
  140,
  156,
  110,
  150,
  160,
  120,
  140,
  120,
  110,
  120,
  140,
  160,
  160,
  110,
  150,
  118,
  110,
  120,
  120,
  146,
  124,
  170,
  124,
  170,
  159,
  120,
  120,
  118,
  152,
  190
)

pdf(file="TTT_app_2.pdf",width=7,height=7, paper="special",
    family="Bookman",pointsize=16)
AdequacyModel::TTT(second_application, lwd = 4)
dev.off()


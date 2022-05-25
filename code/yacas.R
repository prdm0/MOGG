library(Ryacas)

G <- "(1 + (x/beta) ^ (-alpha)) ^ (-p)"
g <- y_fn(G, "D(x)") |> yac_str()

fun <- function(x, alpha, beta, p){
  eval(yac_expr(g))
}

Fun <- function(x, alpha, beta, p){
  eval(yac_expr(G))
}

paste0("theta * (-Ln(1 - (", G, ")))^(alpha - 1) * ", g) |> 
  y_fn(fn = "Simplify") |> 
  yac_str() |> 
  y_fn(fn = "TexForm") |> 
  yac_str()

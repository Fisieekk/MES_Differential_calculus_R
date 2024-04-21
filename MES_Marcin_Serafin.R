library(Matrix)
library(ggplot2)
n<-7
h<-3/(n-2)
G<-6.674*(10^(-11))
# funkcja rho
rho <- function(x){
  if(x<=1 || x>2){return (0)}
  else{return (1)}
}
x_i<-function(i){
  return(h*i)
}
u <- function(x) {
  return(5 - x / 3)
}

# funkcje bazowe z metody Galerkina
e<-function(x, i){
  if(x<x_i(i-1) || x>x_i(i+1)){
      return (0)
    }
  else if (x<=x_i(i)){
    return ((x-x_i(i-1))/x_i(1))
  }
  else { 
    return ((x_i(i+1)-x)/x_i(1))
  }
}
# pochodne funkcji bazowych
de<-function(x, i){
  if(x<x_i(i-1) || x>x_i(i+1)){
    return (0)
  }
  else if (x<=x_i(i)){
    return ((1/(x_i(1)) ))
  }
  else {
    return (-1/(x_i(1)))
  }
}
# iloczyn pochodnych
du_dv=function(i, j){
  return (function(x){
    return(de(x, i)*de(x, j))
  }
  )
}
# Calkowanie metodą Gaussa-Legendre które jednak pomijam ze względu na wbudowaną funckcje całkującą
#integrater <- function(f, a, b) {
#  # Wagi i punkty Gaussa-Legendre dla dwóch punktów kwadratury
#  weights=1.0
#  min_node <- (0.5 * (b - a) * -sqrt(1/3) + 0.5 * (a + b))
#  plus_node<- (0.5 * (b - a) * sqrt(1/3) + 0.5 * (a + b))
#  result_nodes <- weights*(f(min_node)+f(plus_node))
#  result <- 0.5 * (b - a) * result_nodes
#  return(result)
#}
# funkcja phi
gravity=function(i){
  return (function(x){
    return(4*pi*G*rho(x)*e(x,i))
  }
  )
}


B<- function(i,j){
  from = x_i(i-1)
  to = x_i(i+1)
  if(j==-1){
    return((1/3)*(e(i,to)-e(i,from)))
  }
    return(-integrate(du_dv(i,j), from, to))
}
L<-function(i){
  from = x_i(i-1)
  to = x_i(i+1)
  return(integrate(gravity(i),from,to)-B(i,-1))
}

create_matrix <- function() {
  mat <- matrix(0, n-2, n-2)
  mat[1, 1] <- B(1, 1)
  mat[n-2, n-2] <- B(n-2, n-2)
  for (i in 1:(n-3)) {
    mat[i, i+1] <- B(i+1, i+2)
    mat[i+1, i] <- mat[i, i+1]
    if(mat[i,i]!=mat[1,1]){
      mat[i,i]=mat[1,1]
    }
  }
  return(mat)
}
vec <- function() {
  C <- vector()
  for (i in 1:(n-2)){
    C=c(C,L(i+1))
  }
  return(C)
}
solver <- function(){
  M<-create_matrix()
  C<-vec()
  vs = solve(M, C)
  return (vs)
}
plot_result <- function(){
  x=seq(h, 3-h, length.out = n-2)
  y = solver()
  u_values=mapply(u,x)
  y<-y+u_values
  vec_x=c(0,x,3)
  vec_y=c(5,y,4)
  data <- data.frame(x = vec_x, y = vec_y)
  ggplot(data, aes(vec_x, vec_y)) +
    geom_line(color = 'red') +
    labs(title = 'Gravitational Potential FEM Solution', x = 'x', y = 'φ(x)') +
    theme_minimal() +
    theme(panel.grid = element_blank())
}
plot_gradient <- function(x, y, col=NULL, ...) {
  y <- t(apply(y, 1, cumsum))
  if(is.null(col)) col <- rainbow(ncol(y))
  
  x.poly <- c(x, rev(x))
  y <- cbind(null= rep(0, length(x)), y)
  plot(y[,ncol(y)] ~ x, type="n", ...)
  
  for(i in 2:ncol(y)){
    y.poly <- c(y[,i-1], rev(y[,i]))
    polygon(x.poly, y.poly, col=col[i-1])
  }
  
}

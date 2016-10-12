plot_gradient <- function(x, y, col=NULL, ...) {
  if(is.null(dim(y))) {
    y <- cbind(null = rep(0, length(x)), y)
  } else {
    y <- t(apply(y, 1, cumsum))
    y <- cbind(null = rep(0, length(x)), y)
  }
  
  if(is.null(col)) col <- rainbow(ncol(y)-1)
  
  x.poly <- c(x, rev(x))
  
  plot(y[,ncol(y)] ~ x, type="n", ...)
  for(i in 2:ncol(y)){
    y.poly <- c(y[,i-1], rev(y[,i]))
    polygon(x.poly, y.poly, col=col[i-1])
  }
  
}


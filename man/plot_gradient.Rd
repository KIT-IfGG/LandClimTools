\name{plot_gradient}
\alias{plot_gradient}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot succession or elevation gradient
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot_gradient(x, y, col = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
x values, either elevation or decades from e.g. elevation output.
}
  \item{y}{
y values as matrix. Species columns from e.g. elevation output. 
}
  \item{col}{
Color vector of length ncol(y).
}
  \item{\dots}{
Additional arguments for function plot().
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

## The function is currently defined as
function (x, y, col = NULL, ...) 
{
    y <- t(apply(y, 1, cumsum))
    if (is.null(col)) 
        col <- rainbow(ncol(y))
    x.poly <- c(x, rev(x))
    y <- cbind(null = rep(0, length(x)), y)
    plot(y[, ncol(y)] ~ x, type = "n", ...)
    for (i in 2:ncol(y)) {
        y.poly <- c(y[, i - 1], rev(y[, i]))
        polygon(x.poly, y.poly, col = col[i - 1])
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

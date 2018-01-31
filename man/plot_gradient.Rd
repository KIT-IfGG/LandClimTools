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
Color vector of length ncol(y). Use e.g. landclim_colors().
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

dat <- read.table(system.file("elevation_biomass_out.csv", package = "LandClimTools"), sep = ",", dec = ".", header = TRUE)

species <- c("abiealba" , "piceabie", "fagusilv", "pinusilv", "querpetr")

plot_gradient(dat$decade[dat$elevation==823], dat[dat$elevation==823, c(species)], col=landclim_colors(length(species)), xlab="Decade", ylab="Biomass", ylim=c(0,400))
legend("topleft", legend = species, col = landclim_colors(length(species)), pch=15, pt.cex=2)

plot_gradient(dat$elevation[dat$decade==50], dat[dat$decade==50, c(species)], xlab="Elevation", ylab="Biomass", col=landclim_colors(length(species)), ylim=c(0,400))
legend("topleft", legend = species, col = landclim_colors(length(species)), pch=15, pt.cex=2)



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

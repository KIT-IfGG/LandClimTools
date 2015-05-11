\name{plot_elevation_gradient}
\alias{plot_elevation_gradient}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Plot elevation gradient
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Create standart figure for elevation gradient for selected decade based on elevation aggregated LandClim output file.
}
\usage{
plot_elevation_gradient(elevationBiomassOut, species, selection = 10, lty = 1, cols = rainbow(length(species)), plotlegend = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{elevationBiomassOut}{
%%     ~~Describe \code{elevationBiomassOut} here~~
}
  \item{species}{
%%     ~~Describe \code{species} here~~
}
  \item{selection}{
%%     ~~Describe \code{selection} here~~
}
  \item{lty}{
%%     ~~Describe \code{lty} here~~
}
  \item{cols}{
%%     ~~Describe \code{cols} here~~
}
  \item{plotlegend}{
%%     ~~Describe \code{plotlegend} here~~
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
\code{\link{plot_forest}}
}
\examples{
dat <- read.table(system.file("elevation_biomass_out.csv", package = "landclimtools"), sep=",", dec=".", header=T)

species <- c("abiealba" , "piceabie", "fagusylv", "pinusilv", "querpetr")
plot_elevation_gradient(elevationBiomassOut=dat, species=species, selection=30, lty=1,  cols= rainbow(length(species)))  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

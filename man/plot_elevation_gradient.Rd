\name{plot_elevation_gradient}
\alias{plot_elevation_gradient}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Plot elevation gradient
}
\description{
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
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{plot_forest}}
}
\examples{
  dat <- read.table(system.file("elevation_biomass_out.csv", package = "LandClimTools"), sep=",", dec=".", header=TRUE)

  species <- c("abiealba" , "piceabie", "fagusylv", "pinusilv", "querpetr")
  plot_elevation_gradient(elevationBiomassOut=dat, species=species, selection=30, lty=1,  cols= rainbow(length(species)))
}

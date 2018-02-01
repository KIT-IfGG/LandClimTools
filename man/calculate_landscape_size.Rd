\name{calculate_landscape_size}
\alias{calculate_landscape_size}

\title{
Calculate landscape size
}
\description{
Function to caculate the landscape size using the LandClim input digital elevation ascii-file (e.g. dem.asc).
}

\usage{
  calculate_landscape_size(dem)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dem}{
    Path to digital elevation model file (e.g. dem.asc) in the simulation input folder.
  }
}
\value{the calculated landscape size}
\examples{
 ## calculate_landscape_size(dem="simulations/input/dem.asc")
}

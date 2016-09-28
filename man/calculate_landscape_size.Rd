\name{calculate_landscape_size}
\alias{calculate_landscape_size}

\title{Calculate landscape size}
\description{Function to caculate the landscape size using the LandClim input digital elevation ascii-file (e.g. dem.asc).}

\usage{
  calculate_landscape_size(dem)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dem}{
    %%     ~~Describe \code{dem} here~~
  }
}
\value{the calculated landscape size}
\examples{
  ## The function is currently defined as
  function(dem){
    return(prod(dim(dem)) * prod(res(dem)) / (100*100))
  }
}

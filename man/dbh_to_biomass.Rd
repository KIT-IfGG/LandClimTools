\name{dbh_to_biomass}
\alias{dbh_to_biomass}

\title{LandClim allomentry for biomass - DBH conversion.}
\description{Converts DBH to biomass}

\usage{
  dbh_to_biomass(dbh, leafHabit, allometry = "SCHUMACHER")
}
\arguments{
  \item{dbh}{
    % TODO: Document this parameter
  }
  \item{leafHabit}{
    LandClim leaf habit, e.g. c("EVERGREEN", "BROADLEAFEVERGREEN", "DECIDUOUS")
  }
  \item{allometry}{
    Type of LandClim allometry to use, e.g. "SCHUMACHER", "POWER".
  }
}
\value{the value converted to biomass}

\seealso{
  \code{\link{biomass_to_dbh}}
}
\examples{
  dat <- rnorm(10, 20, 5)
  dbh_to_biomass(dbh=dat, leafHabit="DECIDUOUS", allometry = "SCHUMACHER")
}

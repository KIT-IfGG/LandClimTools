\name{biomass_to_dbh}
\alias{biomass_to_dbh}

\title{LandClim allomentry for biomass to DBH conversion.}
\description{Converts biomass to DBH}

\usage{
  biomass_to_dbh(biomass, leafHabit, allometry = "SCHUMACHER")
}
\arguments{
  \item{biomass}{
    Numeric vector containing biomass in tonns.
  }
  \item{leafHabit}{
    LandClim leaf habit, e.g. "EVERGREEN", "BROADLEAFEVERGREEN", "DECIDUOUS".
  }
  \item{allometry}{
    Type of LandClim allometry to use, e.g. "SCHUMACHER", "POWER".
  }
}
\value{the value converted to DBH}
\seealso{
  \code{\link{dbh_to_biomass}}
}
\examples{
  biomass_to_dbh(biomass=3, leafHabit="DECIDUOUS", allometry = "SCHUMACHER")
}

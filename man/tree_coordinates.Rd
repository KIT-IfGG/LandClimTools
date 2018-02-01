\name{tree_coordinates}
\alias{tree_coordinates}

\title{
Calculate tree coordinates.
}
\description{

}
\usage{
  tree_coordinates(file, a = 25, biomasslargetrees = 10, decade = 30, oldtrees = NULL, silent = TRUE)
}

\arguments{
  \item{file}{LandClim output file, e.g. fullOut_30.csv.}
  \item{a}{Resolution of LandClim cells in meters, e.g. 25.}
  \item{biomasslargetrees}{Biomass in tons. Trees larger than the given value will be placed on a wider grid than the others.}
  \item{decade}{
    
  }
  \item{oldtrees}{Optional. Tree coordinates created by the same function, but "before" to preserve ids and coordinates. 
          
  }
  \item{silent}{
    TRUE/FALSE
  }
}
\value{
  Data frame with coordinates for trees, ready to use in \code{\link{plot_forest}}.
}
\references{

}
\author{
  Klara Dolos
}


\seealso{
  \code{\link{plot_forest}}
}
\examples{
trees <- tree_coordinates(file=system.file("fullOut_50.csv", package = "LandClimTools"), a=25)

### Function also creates tree ID for su
## oldtrees <- tree_coordinates(file=paste(path, "fullOut_49.csv", sep=""), a=25, decade=39, oldtrees=NULL)
## trees <- tree_coordinates(file=paste(path, "fullOut_50.csv", sep=""), a=25, decade=40, oldtrees=oldtrees)

}

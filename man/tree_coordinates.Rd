\name{tree_coordinates}
\alias{tree_coordinates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate tree coordinates.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
  tree_coordinates(file, a = 25, biomasslargetrees = 10, decade = 30, oldtrees = NULL, silent = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{LandClim output file, e.g. fullOut_30.csv.}
  \item{a}{Resolution of LandClim cells in meters, e.g. 25.}
  \item{biomasslargetrees}{Biomass in tons. Trees larger than the given value will be placed on a wider grid than the others.}
  \item{decade}{
    % TODO: Document this parameter
  }
  \item{oldtrees}{
    % TODO: Document this parameter
  }
  \item{silent}{
    % TODO: Document this parameter
  }
}
\value{
  Data frame with coordinates for trees, ready to use in \code{\link{plot_forest}}.
}
\references{
  %% ~put references to the literature/web site here ~
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

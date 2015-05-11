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
tree_coordinates(file = "Data/Output/fullOut_30.csv", a = 25, biomasslargetrees = 10, silent=T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{LandClim output file, e.g. fullOut_30.csv.
%%     ~~Describe \code{file} here~~
}
  \item{a}{Resolution of LandClim cells, e.g. 25 m.
%%     ~~Describe \code{a} here~~
}
  \item{biomasslargetrees}{Biomass (t). Trees larger than the given value will be placed on a wider grid than the others.
%%     ~~Describe \code{biomasslargetrees} here~~
}
  \item{silent}{Default is TRUE.
%%     ~~Describe \code{biomasslargetrees} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
\note{
%%  ~~further notes~~
}


\seealso{
\code{\link{plot_forest}}
}
\examples{
trees <- tree_coordinates(file=system.file("fullOut_50.csv", package = "landclimtools"), a=25)

### Function also creates tree ID for su
## oldtrees <- tree_coordinates(file=paste(path, "fullOut_49.csv", sep=""), a=25, decade=39, oldtrees=NULL)
## trees <- tree_coordinates(file=paste(path, "fullOut_50.csv", sep=""), a=25, decade=40, oldtrees=oldtrees)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

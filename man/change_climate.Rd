\name{change_climate}
\alias{change_climate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Simulate climate change.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The funcion reads a LandClim climate file and adds a constant to temperature and precipitation.
}
\usage{
change_climate(inputPath, outputPath, dt = 0, dn = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputPath}{
%%     ~~Describe \code{inputPath} here~~
Path to LandClim climate file.
}
  \item{outputPath}{
%%     ~~Describe \code{outputPath} here~~
Path to LandClim climate changed file.
}
  \item{dt}{
%%     ~~Describe \code{dt} here~~
Temperature in K to add. Also negative values allowed.
}
  \item{dn}{
%%     ~~Describe \code{dn} here~~
Precipitation to add in mm.  Also negative values allowed.
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
Klara Dolos
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (inputPath, outputPath, dt = 0, dn = 0)
{
    header <- readLines(inputPath)
    header <- header[1:14]
    clim <- read.table(inputPath, skip = 11)
    clim[, 2:13] <- clim[, 2:13] + dt
    clim[, 14:25] <- clim[, 14:25] + dn
    clim[, 14:25][clim[, 14:25] < 0] <- 0
    writeLines(header, outputPath)
    write.table(clim, outputPath, append = TRUE, row.names = FALSE,
        col.names = FALSE)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

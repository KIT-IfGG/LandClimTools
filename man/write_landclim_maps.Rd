\name{write_LandClim_maps}
\alias{write_LandClim_maps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Write landClim maps
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Takes raster stack in the appropiate resolution and writes LandClim maps in the required format.
}
\usage{
write_LandClim_maps(LandClimRasterStack, nodata_value = "-9999", lcResolution = 25, ex)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{LandClimRasterStack}{
%%     ~~Describe \code{LandClimRasterStack} here~~
}
  \item{nodata_value}{
%%     ~~Describe \code{nodata_value} here~~
}
  \item{lcResolution}{
%%     ~~Describe \code{lcResolution} here~~
}
  \item{ex}{
%%     ~~Describe \code{ex} here~~
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (LandClimRasterStack, nodata_value = "-9999", lcResolution = 25, 
    ex) 
{
    LandClimRasterStack_list <- lapply(unstack(LandClimRasterStack), 
        function(x) crop(x, ex))
    rs <- stack(LandClimRasterStack_list)
    names(rs) <- names(LandClimRasterStack)
    writeRaster(rs, "landClimMaps.tif", overwrite = T)
    rm(rs)
    foo <- function(x) {
        sink(paste(names(x), ".asc", sep = ""))
        writeLines(c(paste("ncols", ncol(x)), paste("nrows", 
            nrow(x)), paste("xllcorner", xmin(x)), paste("yllcorner", 
            ymin(x)), paste("cellsize", lcResolution), paste("nodata_value ", 
            nodata_value)))
        sink()
        write.table(matrix(round(x[]), nrow = nrow(x), ncol = ncol(x), 
            byrow = T), file = paste(names(x), ".asc", sep = ""), 
            append = T, quote = FALSE, row.names = FALSE, col.names = FALSE)
    }
    lapply(LandClimRasterStack_list, function(x) foo(x))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

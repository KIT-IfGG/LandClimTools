\name{resample_LandClim_maps}
\alias{resample_LandClim_maps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Resample LandClim maps
}
\description{
Function to resample LandClim maps of any resolution to required resolution of 25 m (or 30 m).
}
\usage{
resample_LandClim_maps(LandClimRasterStack, targetResolution = 25)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{LandClimRasterStack}{ Raster stack (with all required input maps)
%%     ~~Describe \code{LandClimRasterStack} here~~
}
  \item{targetResolution}{
%%     ~~Describe \code{targetResolution} here~~
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
Raster stack ready to be written in LandClim format.
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
\code{\link{write_LandClim_maps}}
}
\examples{
gk_projection<-CRS("+init=epsg:31467")
require(raster)
nr <- 20
nc <- 20
res <- 45
ex <- extent(0, nc*res, 0, nr*res)
dem <- raster(nrows=nr, ncols=nc, ex)
projection(dem) <- gk_projection
dem[] <- rep(seq(500, 2200,len=nr), each=nc) 

### LandClim map "slope" and "aspect".
slope <- terrain(dem, filename="slopeAspect.tif", opt='slope', unit="degrees", overwrite = T)
slope[]<- 0

###  LandClim map "soil".
soil <- dem
soil[] <- 20
soil  ### Check min, max values

###  LandClim map "landtype".
landtype <- slope
landtype[] <- 1  

### Aspect (dummy)
aspect <- slope
aspect[] <- 0

###  LandClim map "nitrogen".
nitro <- slope
nitro[] <- 1

### Create raster-stack
maps <- stack(dem, slope, aspect, soil, landtype, nitro)
names(maps) <- c("dem", "slope", "aspect", "soil", "landtype", "nitro")

maps25 <- resample_LandClim_maps(LandClimRasterStack=maps)
plot(maps25$dem)

}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

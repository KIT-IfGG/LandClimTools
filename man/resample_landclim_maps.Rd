\name{resample_landclim_maps}
\alias{resample_landclim_maps}

\title{
Resample LandClim maps
}
\description{
Function to resample LandClim maps of any resolution to required resolution of 25 m (or 30 m).
}

\usage{
  resample_landclim_maps(landClimRasterStack, targetResolution = 25)
}
\arguments{
  \item{landClimRasterStack}{Raster stack (with all required input maps)}
  \item{targetResolution}{
    
  }
}
\value{Raster stack ready to be written in LandClim format.}

\seealso{
  \code{\link{write_landclim_maps}}
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
slope <- terrain(dem, filename="slopeAspect.tif", opt='slope', unit="degrees", overwrite = TRUE)
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

maps25 <- resample_landclim_maps(landClimRasterStack=maps)
plot(maps25$dem)

}

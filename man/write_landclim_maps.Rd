\name{write_landclim_maps}
\alias{write_landclim_maps}

\title{

Write landClim maps
}
\description{
Takes raster stack in the appropiate resolution and writes LandClim maps in the required format.
}
\usage{
write_landclim_maps(landClimRasterStack, nodata_value = "-9999", lcResolution = 25, folder)
}

\arguments{
  \item{landClimRasterStack}{
    
  }
  \item{nodata_value}{
    In LandClim often 999. Or another value that should be used.
  }
  \item{lcResolution}{
    Cell edge length (e.g. 25)
  }
  \item{folder}{
    Path to folder in which files will be written.
  }
}
\details{

}
\value{

}
\references{

}
\author{
Klara Dolos
}



\seealso{

}
\examples{
### Create and write LandClim maps ####
require(raster)
gk_projection<-CRS("+init=epsg:31467")
nr <-50
nc <- 50
res <- 40
ex <- extent(0, nc*res, 0, nr*res)
dem <- raster(nrows=nr, ncols=nc, ex)
projection(dem) <- gk_projection
dem
dem[] <- rep(seq(400, 2200,len=nr), each=nc)
x11()
plot(dem)

### Create LandClim map "slope".
slope <- dem
slope[]<- 0

###  LandClim map "soil".
soil <- dem
soil[] <- 20
soil  ### Check min, max values

###  LandClim map "landtype".
landtype <- slope
landtype[] <- 1

### Aspect
aspect <- slope
aspect[] <- 0

###  LandClim map "nitrogen".
nitro <- slope
nitro[] <- 1

### Create raster-stack
maps <- stack(dem, slope, aspect, soil, landtype, nitro)
names(maps) <- c("dem", "slope", "aspect", "soil", "landtype", "nitro")
x11()
plot(maps)

maps25 <- resample_landclim_maps(landClimRasterStack=maps)
res(maps25$dem)

### Write as LandClim files.
write_landclim_maps(landClimRasterStack=maps25, nodata_value="-9999", lcResolution=25)


}


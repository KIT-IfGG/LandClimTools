\name{LandClimTools-package}
\alias{LandClimTools-package}
\alias{LandClimTools}
\docType{package}

\title{
Package documentation for LandClimTools
}
\description{
The LandClimTools package contains several useful functions for working with the LandClim software
}

\references{
  Website at ETH Zurich giving an introduction to LandClim:\cr
  \url{https://www1.ethz.ch/fe/research/disturbance/landclim}

  Website about the LandClim software:\cr
  \url{https://uwis-server102.ethz.ch/openaccess/software/view/2}

  Style guide for R packages by Hadley Wickham:\cr
  \url{http://r-pkgs.had.co.nz/style.html}

  \href{https://dx.doi.org/10.1016/j.ecolmodel.2003.12.055}{Schumacher, S., H. Bugmann, and D. J. Mladenoff. 2004. Improving the formulation of tree growth and succession in a spatially explicit landscape model. Ecological Modelling 180:175-194.}

  \href{https://dx.doi.org/10.1111/j.1365-2486.2006.01188.x}{Schumacher, S. and H. Bugmann. 2006. The relative importance of climatic effects, wildfires and management for future forest landscape dynamics in the Swiss Alps. Global Change Biology 12:1435-1450.}

  \href{https://dx.doi.org/10.1007/s10980-005-2165-7}{Schumacher, S., B. Reineking, J. Sibold, and H. Bugmann. 2006. Modeling the impact of climate and vegetation on fire regimes in mountain landscapes. Landscape Ecology 21:539-554.}
}

\keyword{ package }


\examples{
###############################################################
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

################################################################### Plot LandClim output ####
### Elevation gradient
dat <- read.table(system.file("elevation_biomass_out.csv", package = "LandClimTools"), sep=",", dec=".", header=TRUE)
species <- c("abiealba" , "piceabie", "fagusylv", "pinusilv", "querpetr")
x11()
plot_elevation_gradient(elevationBiomassOut=dat, species=species, selection=30, lty=1,  cols= rainbow(length(species)))

### LandClim forest
trees <- tree_coordinates(file=system.file("fullOut_50.csv", package = "LandClimTools"), a=25)

stand <- trees[trees$row > 20 & trees$row <=40,]
stand$row <- stand$row - min(stand$row)
stand <- trees[trees$col > 20 & trees$col <=40,]
stand$col <- stand$col - min(stand$col)

x11(width=7, height=7)
par(mar=c(2,2,1,1))
plot_forest(trees=stand, species=unique(stand$species),  scol=rainbow(length(unique(stand$species))), plotlegend=TRUE, aspect=1, cex=sqrt(stand$biomass)/2)
}

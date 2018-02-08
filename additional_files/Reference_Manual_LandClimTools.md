

<!-- toc -->

februari 08, 2018

# DESCRIPTION

```
Package: LandClimTools
Type: Package
Title: Tools for Modelling Landscapes
Description: Some useful functions to create LandClim input files, analyse and visualize output from the LandClim software by ETH Zurich for modelling landscapes.
URL: https://www1.ethz.ch/fe/research/disturbance/landclim
  https://uwis-server102.ethz.ch/openaccess/software/view/2
Version: 1.0.0
Date: 2016-08-08
Encoding: UTF-8
Author: Klara Dolos, Caspar Groß
Maintainer: Klara Dolos <dolos@ifgg.kit.edu>
License: GPL-3
Depends:
  raster,
  rgdal,
  XML```


# `LandClimTools-package`: 
 Package documentation for LandClimTools


## Description


 The LandClimTools package contains several useful functions for working with the LandClim software


## References


 Website at ETH Zurich giving an introduction to LandClim: list() 
  [https://www1.ethz.ch/fe/research/disturbance/landclim](https://www1.ethz.ch/fe/research/disturbance/landclim) 
 
 Website about the LandClim software: list() 
  [https://uwis-server102.ethz.ch/openaccess/software/view/2](https://uwis-server102.ethz.ch/openaccess/software/view/2) 
 
 Style guide for R packages by Hadley Wickham: list() 
  [http://r-pkgs.had.co.nz/style.html](http://r-pkgs.had.co.nz/style.html) 
 
  [Schumacher, S., H. Bugmann, and D. J. Mladenoff. 2004. Improving the formulation of tree growth and succession in a spatially explicit landscape model. Ecological Modelling 180:175-194.](https://dx.doi.org/10.1016/j.ecolmodel.2003.12.055) 
 
  [Schumacher, S. and H. Bugmann. 2006. The relative importance of climatic effects, wildfires and management for future forest landscape dynamics in the Swiss Alps. Global Change Biology 12:1435-1450.](https://dx.doi.org/10.1111/j.1365-2486.2006.01188.x) 
 
  [Schumacher, S., B. Reineking, J. Sibold, and H. Bugmann. 2006. Modeling the impact of climate and vegetation on fire regimes in mountain landscapes. Landscape Ecology 21:539-554.](https://dx.doi.org/10.1007/s10980-005-2165-7) 


## Examples

```r 
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
 ``` 

# `biomass_to_dbh`: 
 LandClim allomentry conversion


## Description


 Converts biomass to DBH using LandClim allomentries


## Usage

```r
  biomass_to_dbh(biomass, leafHabit, allometry = "SCHUMACHER")
```


## Arguments

Argument      |Description
------------- |----------------
```biomass```     |      Numeric vector containing biomass in tonns. 
```leafHabit```     |      LandClim leaf habit, e.g. "EVERGREEN", "BROADLEAFEVERGREEN", "DECIDUOUS". 
```allometry```     |      Type of LandClim allometry to use, e.g. "SCHUMACHER", "POWER". 

## Value

the value converted to DBH

## Seealso


  [`dbh_to_biomass`](dbh_to_biomass.html) 


## Examples

```r 
 biomass_to_dbh(biomass=3, leafHabit="DECIDUOUS", allometry = "SCHUMACHER")
 ``` 

# `calculate_landscape_size`: 
 Calculate landscape size


## Description


 Function to caculate the landscape size using the LandClim input digital elevation ascii-file (e.g. dem.asc).


## Usage

```r
  calculate_landscape_size(dem)
```


## Arguments

Argument      |Description
------------- |----------------
```dem```     |      Path to digital elevation model file (e.g. dem.asc) in the simulation input folder. 

## Value

the calculated landscape size

## Examples

```r 
 ## calculate_landscape_size(dem="simulations/input/dem.asc")
 ``` 

# `change_climate`: 
 Simple climate change


## Description


 The funcion reads a LandClim climate file and adds a constant to temperature and precipitation.


## Usage

```r
  change_climate(inputPath, outputPath, dt = 0, dn = 0)
```


## Arguments

Argument      |Description
------------- |----------------
```inputPath```     |      Path to LandClim climate file. 
```outputPath```     |      Path to LandClim climate changed file. 
```dt```     |      Temperature in K to add. Also negative values allowed. 
```dn```     |      Precipitation to add in mm.  Also negative values allowed. 

## Author

Klara Dolos

## Examples

```r 
 ## The function is currently defined as
 function (inputPath, outputPath, dt = 0, dn = 0) {
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
 ``` 

# `create_movie`: 
 Create succession movie


## Description


 Creates a sequence of figures of the forest state from above displayed as a "movie".


## Usage

```r
  create_movie(files, decades, a, species,  scol, plotlegend=TRUE, aspect=1, silent=FALSE)
```


## Arguments

Argument      |Description
------------- |----------------
```files```     |      LandClim decade output "fullOut10" etc. One file per time step for the movie. 
```decades```     |      Vector, Decade corresponding to files, used as labels. 
```a```     |      LandClim cell size. 
```species```     |      Simulated species, names used as column names in the LandClim output. 
```scol```     |      Color vector. See also ?landclim_colors. 
```plotlegend```     |      TRUE or FALSE 
```aspect```     |      Function argument asp in plot(). 
```silent```     |      TRUE/FALSE 

## Examples

```r 
 
 ## full <- read.csv(file = "fullOut_30.csv", strip.white = TRUE)
 ## species <- unique(full$species)
 
 ## fi <- data.frame(file = c("fullOut_30.csv", "fullOut_31.csv","fullOut_32.csv","fullOut_33.csv","fullOut_34.csv"), decade = c(30, 31, 32, 33,34))
 
 ## pdf(paste(path,"movie.pdf", sep=""))
 ## par(mfrow=c(2,2))
 ## create_movie(files=fi$file, decades=fi$decade, a=25, species=species,  scol=rainbow(length(species)), plotlegend=TRUE, aspect=1)
 ## dev.off()
 
 ## library(animation)
 ## saveGIF(create_movie(files = fi$file, decades = fi$decade, a = 25, species = species,  scol = rainbow(length(species)), plotlegend=TRUE, aspect=1), movie.name = "LandClimForest.gif")
 
 ``` 

# `dbh_to_biomass`: 
 LandClim allomentry for biomass - DBH conversion


## Description


 Converts DBH to biomass


## Usage

```r
  dbh_to_biomass(dbh, leafHabit, allometry = "SCHUMACHER")
```


## Arguments

Argument      |Description
------------- |----------------
```dbh```     |       % TODO: Document this parameter  
```leafHabit```     |      LandClim leaf habit, e.g. c("EVERGREEN", "BROADLEAFEVERGREEN", "DECIDUOUS") 
```allometry```     |      Type of LandClim allometry to use, e.g. "SCHUMACHER", "POWER". 

## Value

the value converted to biomass

## Seealso


  [`biomass_to_dbh`](biomass_to_dbh.html) 


## Examples

```r 
 dat <- rnorm(10, 20, 5)
 dbh_to_biomass(dbh=dat, leafHabit="DECIDUOUS", allometry = "SCHUMACHER")
 ``` 

# `global_coordinates`: 
 Local to global coordinates


## Description


 Takes local coordinates, the numbers of rows/columns and the width/height of each part of the grid. From that the global coordinates are calculated.


## Usage

```r
global_coordinates(x.local, y.local, row, col, a)
```


## Arguments

Argument      |Description
------------- |----------------
```x.local```     |     local x coordinate
```y.local```     |     local y coordinate
```row```     |     row index
```col```     |     column index
```a```     |     size (both width and height) of each local cell

## Value


 Global coordinates data frame


# `landclim_colors`: 
 LandClim color palette


## Description


 Create a color vector for plotting.


## Usage

```r
landclim_colors(n)
```


## Arguments

Argument      |Description
------------- |----------------
```n```     |      Number of different tree species = number of colors.

## Value


 Character vector with colors.


## Note


 Please feel free to enhance this palette!


## Author


 Klara Dolos


## Examples

```r 
 n <- 10
 pie(rep(1, n), col=landclim_colors(n))
 
 ### Get number of available colors.
 landclim_colors(99)
 
 ``` 

# `plot_forest`: 
 Plot LandClim forest


## Description


 Plots a forest based on tree data from LandClim.


## Usage

```r
  plot_forest(trees, species = unique(trees$species), scol = rainbow(length(species)), plotlegend = TRUE, a = 25, aspect = 1, cex = 1)
```


## Arguments

Argument      |Description
------------- |----------------
```trees```     |     Object created by function tree_coordinates.
```species```     |     Species names, character vector.
```scol```     |     Color vector.
```plotlegend```     |     TRUE or FALSE indicating whether to plot a legend.
```a```     |     Resolution of LandClim maps used in the simulation.
```aspect```     |     Argument asp in function plot.
```cex```     |     Argument cex in function plot.

## Seealso


  [`tree_coordinates`](tree_coordinates.html) 


## Author


 Klara Dolos


## Examples

```r 
 
 trees <- tree_coordinates(file=system.file("fullOut_50.csv", package = "LandClimTools"), a=25)
 
 range(trees$col)
 stand <- trees[trees$row > 20 & trees$row <=40,]
 stand$row <- stand$row - min(stand$row)
 
 stand <- trees[trees$col > 20 & trees$col <=40,]
 stand$col <- stand$col - min(stand$col)
 
 x11(width=7, height=7)
 par(mar=c(2,2,1,1))
 plot_forest(trees=stand, species=unique(stand$species),  scol=rainbow(length(unique(stand$species))), plotlegend=TRUE, aspect=1, cex=sqrt(stand$biomass)/2)
 
 ``` 

# `plot_gradient`: 
 Gradient plot


## Description


 Plot succession or elevation gradient based on Landclim decadal output.


## Usage

```r
plot_gradient(x, y, col = NULL, ...)
```


## Arguments

Argument      |Description
------------- |----------------
```x```     |      x values, either elevation or decades from e.g. elevation output.
```y```     |      y values as matrix. Species columns from e.g. elevation output.
```col```     |      Color vector of length ncol(y). Use e.g. landclim_colors().
```list()```     |      Additional arguments for function plot().

## Author


 Klara Dolos


## Examples

```r 
 
 dat <- read.table(system.file("elevation_biomass_out.csv", package = "LandClimTools"), sep = ",", dec = ".", header = TRUE)
 
 species <- c("abiealba" , "piceabie", "fagusilv", "pinusilv", "querpetr")
 
 plot_gradient(dat$decade[dat$elevation==823], dat[dat$elevation==823, c(species)], col=landclim_colors(length(species)), xlab="Decade", ylab="Biomass", ylim=c(0,400))
 legend("topleft", legend = species, col = landclim_colors(length(species)), pch=15, pt.cex=2)
 
 plot_gradient(dat$elevation[dat$decade==50], dat[dat$decade==50, c(species)], xlab="Elevation", ylab="Biomass", col=landclim_colors(length(species)), ylim=c(0,400))
 legend("topleft", legend = species, col = landclim_colors(length(species)), pch=15, pt.cex=2)
 
 
 
 ``` 

# `profound_climate_to_landclim`: 
 Profound climate to LandClim


## Description


 Complile climate data from the Profound database to the format needed as input file for LandClim


## Usage

```r
profound_climate_to_landclim(climate, header, file = "climate.txt")
```


## Arguments

Argument      |Description
------------- |----------------
```climate```     |      Dataframe as provided by the PROFOUND database, e.g. ProfoundData::getData("Soro", "CLIMATE").
```header```     |      Header information as provided by readLines("landclim_climate.txt"). I.e. read in an existing LandClim climate input file and provide it as function argument.
```file```     |      Output filename

## Value


 Written text file with landclim climate.


## Author


 Klara Dolos


## References


 The Profound Database is available here: https://github.com/COST-FP1304-PROFOUND/TG2 (NOT YET, Feb 2018, contact Christopher Reyer, PIK).
 
 http://cost-profound.eu/site
 
 Examples for simulations for cost sites (Geoecology Master course 2017/18, Karlsruhe)
 https://github.com/KIT-IfGG/climate_impact_research_2
 


## Examples

```r 
 ### See https://github.com/KIT-IfGG/climate_impact_research_2 for example code
 
 ``` 

# `read_species_xml`: 
 Read species XML


## Description


 Read LandClim species parameters from input file (species.xml).


## Usage

```r
read_species_xml(file)
```


## Arguments

Argument      |Description
------------- |----------------
```file```     |       File name

## Details


 


## Value


 
 Data frame with one row for each species.


## Seealso


 write_species_xml


## Author


 Björn Reineking


## References


 


## Examples

```r 
 ### Read example species file
 species_parameters  <- read_species_xml(system.file("species.xml", package = "LandClimTools"))
 
 ### Write example species file
 write_species_xml(x = species_parameters, file = "piceabi.xml")
 ``` 

# `resample_landclim_maps`: 
 Resample LandClim maps


## Description


 Function to resample LandClim maps of any resolution to required resolution of 25 m (or 30 m).


## Usage

```r
  resample_landclim_maps(landClimRasterStack, targetResolution = 25)
```


## Arguments

Argument      |Description
------------- |----------------
```landClimRasterStack```     |     Raster stack (with all required input maps)
```targetResolution```     |       

## Value

Raster stack ready to be written in LandClim format.

## Seealso


  [`write_landclim_maps`](write_landclim_maps.html) 


## Examples

```r 
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
 
 ``` 

# `run_landclim`: 
 Run LandClim


## Description


 


## Usage

```r
run_landclim(control_file)
```


## Arguments

Argument      |Description
------------- |----------------
```control_file```     |     Name of the control file (e.g. control.xml). 

## Details


 The working directory needs to be set to the level above "Input" and "Output" folders of LandClim, thus at the "site" level. Until now, the folder names are hard-coded, thus need to be "Input" and "Output", not e.g. input! Feel free to fix this!
 
 You need to set the path to your LandClim executable before you run LandClim using this function.
 


## Seealso


 set_landclim_path


## Author


 Klara Dolos


## References


 ### Add LandClim refrences


## Examples

```r 
 
 ``` 

# `tree_coordinates`: 
 Calculate tree coordinates.


## Description


 


## Usage

```r
  tree_coordinates(file, a = 25, biomasslargetrees = 10, decade = 30, oldtrees = NULL, silent = TRUE)
```


## Arguments

Argument      |Description
------------- |----------------
```file```     |     LandClim output file, e.g. fullOut_30.csv.
```a```     |     Resolution of LandClim cells in meters, e.g. 25.
```biomasslargetrees```     |     Biomass in tons. Trees larger than the given value will be placed on a wider grid than the others.
```decade```     |       
```oldtrees```     |     Optional. Tree coordinates created by the same function, but "before" to preserve ids and coordinates.  
```silent```     |      TRUE/FALSE 

## Value


 Data frame with coordinates for trees, ready to use in [`plot_forest`](plot_forest.html) .


## Seealso


  [`plot_forest`](plot_forest.html) 


## Author


 Klara Dolos


## References


 


## Examples

```r 
 trees <- tree_coordinates(file=system.file("fullOut_50.csv", package = "LandClimTools"), a=25)
 
 ### Function also creates tree ID for su
 ## oldtrees <- tree_coordinates(file=paste(path, "fullOut_49.csv", sep=""), a=25, decade=39, oldtrees=NULL)
 ## trees <- tree_coordinates(file=paste(path, "fullOut_50.csv", sep=""), a=25, decade=40, oldtrees=oldtrees)
 
 ``` 

# `write_landclim_maps`: 
 
 Write landClim maps


## Description


 Takes raster stack in the appropiate resolution and writes LandClim maps in the required format.


## Usage

```r
write_landclim_maps(landClimRasterStack, nodata_value = "-9999", lcResolution = 25, folder)
```


## Arguments

Argument      |Description
------------- |----------------
```landClimRasterStack```     |       
```nodata_value```     |      In LandClim often 999. Or another value that should be used. 
```lcResolution```     |      Cell edge length (e.g. 25) 
```folder```     |      Path to folder in which files will be written. 

## Details


 


## Value


 


## Seealso


 


## Author


 Klara Dolos


## References


 


## Examples

```r 
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
 
 
 ``` 

# `write_species_xml`: Write LandClim species parameter file

## Description

Writes the given data into the given file in the XML-format for species of LandClim

## Usage

```r
  write_species_xml(x, file)
```


## Arguments

Argument      |Description
------------- |----------------
```x```     |      Data frame with species parameters and tag names. 
```file```     |      Output file name including full path. 

## Seealso


 read_species_xml


## Examples

```r 
 ### Read example species file
 species_parameters  <- read_species_xml(system.file("species.xml", package = "LandClimTools"))
 
 ### Write example species file
 write_species_xml(x = species_parameters, file = "piceabie.xml")
 ``` 


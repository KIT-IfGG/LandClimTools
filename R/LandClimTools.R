resample_landclim_maps <- function(landClimRasterStack, targetResolution=25){

  r <- raster::raster(landClimRasterStack, layer=1)
  rt <- raster::raster(extent(r), crs=projection(r))
  res(rt) <- targetResolution

  foo <- function(x){
    rre <- resample(x, rt)
    rre[is.na(rre)] <- -9999
    rre
  }
  res <- lapply(unstack(landClimRasterStack), foo)
  stack(res)
}

write_landclim_maps <- function(landClimRasterStack, nodata_value="-9999", lcResolution=25, folder=getwd()) {
  ex <- (extent(landClimRasterStack))
  landClimRasterStack_list <- lapply(unstack(landClimRasterStack), function(x) crop(x, ex))
  rs <- stack(landClimRasterStack_list)
  names(rs) <- names(landClimRasterStack)
  writeRaster(rs, paste0(folder, "/landclim_maps.tif"), overwrite=TRUE)
  rm(rs)
  foo <- function(x){
    sink(paste(folder, "/", names(x), ".asc", sep=""))
    writeLines(c(paste("ncols", ncol(x)), paste("nrows", nrow(x)), paste("xllcorner", xmin(x)), paste("yllcorner", ymin(x)), paste("cellsize", lcResolution), paste("nodata_value ", nodata_value)))
    sink()
    write.table(matrix(round(x[]), nrow=nrow(x), ncol=ncol(x), byrow=TRUE), file=paste(folder, "/", names(x), ".asc", sep=""), append=TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  lapply(landClimRasterStack_list, function(x) foo(x))
}

change_climate <- function(inputPath, outputPath, dt=0, dn=0){
  header<- readLines(inputPath)
  header <- header[1:14]
  clim<- read.table(inputPath, skip=11)
  clim[,2:13] <- clim[,2:13] + dt
  clim[,14:25] <- clim[,14:25] + dn
  clim[,14:25][clim[,14:25] < 0] <- 0
  writeLines(header, outputPath)
  write.table(clim, outputPath, append=TRUE, row.names=FALSE, col.names=FALSE)
}

biomass_to_dbh <- function(biomass, leafHabit, allometry="SCHUMACHER") {
  if(allometry=="SCHUMACHER") {
    if(!all(leafHabit %in% c("EVERGREEN", "BROADLEAFEVERGREEN", "DECIDUOUS"))) stop ("unknown leaf habit")
    dbh <- ifelse(leafHabit == "EVERGREEN", exp(3.800 + 0.451 * log(biomass)), exp(3.708 + 0.475 * log(biomass)))
  }

  if(allometry=="POWER") {
    carbon_kg <- 500 * biomass # assumption: 50% of dry weight is carbon
    dbh <- exp(1.3481648 + 0.3977240 * log(carbon_kg))
  }
  return(dbh)
}

dbh_to_biomass <- function(dbh, leafHabit, allometry="SCHUMACHER") {
  if(allometry=="SCHUMACHER") {
    if(!all(leafHabit %in% c("EVERGREEN", "BROADLEAFEVERGREEN", "DECIDUOUS"))) stop ("unknown leaf habit")
    biomass <- ifelse(leafHabit == "EVERGREEN", exp((log(dbh) - 3.800)/0.451), exp((log(dbh)- 3.708)/0.475))
  }

  if(allometry=="POWER") {
    carbon_kg <- exp((log(dbh) - 1.3481648) / 0.3977240)
    biomass <- carbon_kg/500   # assumption: 50% of dry weight is carbon

  }
  return(biomass)
}

calculate_landscape_size <- function(dem){
  return(prod(dim(dem)) * prod(res(dem)) / (100*100))
}

plot_elevation_gradient <- function(elevationBiomassOut, species, selection=10,   lty=1,  cols= rainbow(length(species)), plotlegend=TRUE){
  a <-elevationBiomassOut$decade==selection
  matplot(elevationBiomassOut$elevation[a], elevationBiomassOut[a,colnames(elevationBiomassOut) %in% species], type="l", lty=lty, col=cols, xlab="Elevation (m a.s.l.)", ylab="Biomass (t/ha)")
  if(plotlegend) legend("topright", legend=species, lty=lty, col=cols, bg="white")
}

read_species_xml <- function(file) {
  data1 <- XML::xmlParse(file)
  XML::xmlToDataFrame(nodes=getNodeSet(data1, "//species//set"))
}

write_species_xml <- function(x, file) {
  names <- colnames(x)
  suppressWarnings(tr <- XML::xmlTree("species"))
  for (i in 1:NROW(x)) {
    tr$addTag("set", close=FALSE)
    for (j in names) { tr$addTag(j, as.character(x[i, j])) }
    tr$closeTag()
  }
  invisible(XML::saveXML(tr$value(), file=file))
}

global_coordinates <- function(x.local, y.local, row, col, a){
  x <- a * row + x.local
  y <- a * col + y.local
  data.frame(x,y)
}

tree_coordinates <- function(file, a=25, biomasslargetrees=10, decade=30, oldtrees=NULL, silent=TRUE){
  full <- read.csv(file, strip.white=TRUE)
  full$birth <- decade - full$age/10
  full$cohortID <- paste(full$row, full$col, full$birth, full$species, sep="_")

  ln <- 10   ### ln^2 = Number of positions within a cell.
  trees <- full[rep(1:nrow(full), full$stems),]
  trees$ID <- paste(trees$cohortID, unlist(lapply(full$stems, function(x) 1:x)), sep="_")
  if(length(oldtrees)>1){
    m <- match(trees$ID, oldtrees$ID)
    trees$x <- oldtrees$x[m]
    trees$y <- oldtrees$y[m]
  } else{
    trees$x <- NA
    trees$y <- NA
  }
  gmax <- global_coordinates(x.local=a,y.local=a,row=max(full$row), col=max(full$col), a=a)

  localgrid <- expand.grid(x.local=seq(1,a, len=ln), y.local=seq(1,a, len=ln))

  cells <- unique(trees$cell)
  for(i in 1:length(cells)){
    if (i %in% seq(100, length(cells), by=100) & !silent) cat(paste(i , " ", sep=""))
    if (sum(trees$cell == i)) {
      stand <- trees[trees$cell == i & is.na(trees$x),]
      if (nrow(stand)>0) {

        stand$x.local <- NA
        stand$y.local <- NA

        largetrees <- stand$biomass > biomasslargetrees

        pos.largetrees <- sample(seq(1,nrow(localgrid), by=4), sum(largetrees), replace=FALSE)
        stand$x.local[largetrees] <- localgrid[pos.largetrees, "x.local"]
        stand$y.local[largetrees] <- localgrid[pos.largetrees, "y.local"]

        pos.smalltrees <- sample((1:nrow(localgrid))[!1:nrow(localgrid) %in% pos.largetrees], nrow(stand)-sum(largetrees), replace=TRUE)

        stand$x.local[!largetrees] <- localgrid[pos.smalltrees, "x.local"]
        stand$y.local[!largetrees] <- localgrid[pos.smalltrees, "y.local"]
        global <- global_coordinates(stand$x.local, stand$y.local, row=stand$row, col=stand$col, a=25)
        trees[trees$cell==i & is.na(trees$x),"x"] <- global$x
        trees[trees$cell==i & is.na(trees$y),"y"] <- global$y

      }

    } else {
      ### For empty cells.
    }
  }
  trees
}


plot_forest <- function(trees, species = unique(trees$species),  scol = rainbow(length(species)), plotlegend=TRUE, a=25, aspect=1, cex=1){
  trees$colors <- scol[match(trees$species, species)]
  plot(y ~ x, trees, pch=16, cex=cex, col=trees$colors, type="p", xlab="", ylab="", asp=aspect)
  if(plotlegend) legend("topright", legend=species, pch=16, col=scol, bg="white")
}


create_movie <- function(files, decades, a, species,  scol, plotlegend=TRUE, aspect=1, silent=FALSE){
  oldtrees <- tree_coordinates(file=files[1], a=a, decade=decades[1], oldtrees=NULL)
  plot_forest(trees=oldtrees, species=species,  scol=scol, plotlegend=plotlegend, aspect=aspect, cex=sqrt(oldtrees$biomass)/2)
  mtext(paste("Decade = ", decades[1], sep=""), side=3, line=0.5, adj = 0)
  for(i in 2:length(files)){
    cat(paste(i , " ", sep=""))
    if (i %in% seq(1, length(files), by=5) & !silent) cat(Sys.time())
    trees <- tree_coordinates(file=files[i], a=a, decade=decades[i], oldtrees=oldtrees)
    plot_forest(trees=trees, species=species,  scol=scol, plotlegend=plotlegend, aspect=aspect, cex=sqrt(trees$biomass)/2)
    mtext(paste("Decade = ", decades[i], sep=""), side=3, line=0.5, adj = 0)
    oldtrees <- trees
  }
}


profound_climate_to_landclim <- function(climate, header, file="climate.txt"){
  temp <- aggregate(climate$tmean_degC, by=list(mo=climate$mo,year= climate$year), mean)
  precip <- aggregate(climate$p_mm, by=list(mo=climate$mo, year=climate$year), sum)

  temp <- round(reshape(temp[,c("year", "mo", "x")],  timevar= "mo", idvar="year", direction = "wide"), 2)
  precip <- round(reshape(precip[,c("year", "mo", "x")],  timevar= "mo", idvar="year", direction = "wide"), 0)

  climate <- cbind(temp[,grep("x", names(temp))], precip[,grep("x", names(temp))])

  writeLines(header, con=file)
  write.table(climate, file=file, append=TRUE, col.names = FALSE, quote=FALSE)
}

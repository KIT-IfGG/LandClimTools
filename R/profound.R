profound_climate_to_landclim <- function(climate, header, file="climate.txt", shuffle = FALSE, n_years=NULL){
  ### Function related to the profound database created in cost action profound 2017. Contact Chr. Reyer, F. Hartig
  ### for access.
  temp <- aggregate(climate$tmean_degC, by=list(mo=climate$mo,year= climate$year), mean)
  precip <- aggregate(climate$p_mm, by=list(mo=climate$mo, year=climate$year), sum)
  
  temp <- round(reshape(temp[,c("year", "mo", "x")],  timevar= "mo", idvar="year", direction = "wide"), 2)
  precip <- round(reshape(precip[,c("year", "mo", "x")],  timevar= "mo", idvar="year", direction = "wide"), 0)
  
  climate <- cbind(temp[,grep("x", names(temp))], precip[,grep("x", names(temp))])
  
  if(shuffle) climate <- climate[sample(1:nrow(climate), size = n_years, replace=TRUE),]
  
  writeLines(header, con=file)
  write.table(climate, file=file, append=TRUE, col.names = FALSE, quote=FALSE)
}

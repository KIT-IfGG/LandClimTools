
lc_path=NULL # assign lc_path a value. Needed to avoid that `R CMD check` complains about "no visible binding for global variable ‘lc_path’"

set_landclim_path <- function(landclim_path){
  if(file.exists(landclim_path)) {
    assign("lc_path", landclim_path, envir=globalenv())
  } else {
    print("Invalid path to LandClim executable.")
    }
}


simulate <- function(control_file){
  ### Todo: Prototype!
  oldwd <- getwd()
  setwd(paste(getwd(), "/Input/", sep=""))
  if(file.exists(control_file)) {
    system(paste(lc_path, control_file, sep=" "))
    setwd(oldwd)
  } else {
    print("Invalid path to LandClim control file.")
  }
}


clean_output_ubuntu <- function(){
  fis <- list.files()
  file.copy(fis[grep("Output", fis)], paste("Output/",fis[grep("Output", fis)], sep=""))
  file.remove(fis[grep(c("csv"), fis)])
  file.remove(fis[grep(c("txt"), fis)])
}

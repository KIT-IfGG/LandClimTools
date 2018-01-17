
lc_path=NULL # assign lc_path a value. Needed to avoid that `R CMD check` complains about "no visible binding for global variable ‘lc_path’"

set_landclim_path <- function(landclim_path){
  if(file.exists(landclim_path)) {
    assign("lc_path", landclim_path, envir=globalenv())
  } else {
    print("Invalid path to LandClim executable.")
    }
}

simulate <- function(control_file){
  print("This is an old function. Use run_landclim() instead.")
  ### Todo: Prototype!
  oldwd <- getwd()
  setwd(paste(getwd(), "/Input/", sep=""))
  if(file.exists(control_file)) {
    system(paste(lc_path, control_file, sep=" "))

  } else {
    print("Invalid path to LandClim control file.")
  }
  setwd(oldwd)
}

run_landclim <- function(control_file = "control.xml") {
  print("Works only for Ubuntu until now!")
  print("Working directory must be at the site level.")
  if (file.exists(paste0("Input/", control_file, sep = "")) & file.exists(lc_path)) {
    oldwd <- getwd()
    dir.create("Output")
    file.remove(list.files("Output", full=TRUE))
    setwd(paste(oldwd, "/Input/", sep = ""))
    system(paste0(lc_path, " ", control_file))
    setwd(oldwd)
    clean_output_ubuntu()
  }
  else {
    print("Invalid path to LandClim executable or to control file.")
    print(paste0("Landclim path: ", lc_path))
    print(paste0("Control file: ", paste(getwd(), "/Input/", control_file, sep = "")))
    print(paste0("Names of *INPUT* and *OUTPUT* folders in the control file unfortunately must be Input and Output:"))
    print(list.files())
  }
}



clean_output_ubuntu <- function(){
  fis <- list.files()
  file.copy(fis[grep("Output", fis)], paste("Output/",fis[grep("Output", fis)], sep=""))
  file.remove(fis[grep(c("csv"), fis)])
  file.remove(fis[grep(c("txt"), fis)])
}

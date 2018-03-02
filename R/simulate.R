run_landclim <- function(control_file = "control.xml", input_folder="Input", output_folder="Output") {
  cat("Works only for Ubuntu until now!\n")
  input_dir <- paste(getwd(), input_folder, sep="/")
  output_dir <- paste(getwd(), output_folder, sep="/")
  control_file <- paste(input_dir, "control.xml", sep="/")
  cat(paste0("Landclim executable: ", landclim_path, "\nInput directory: ", input_dir, "\nOutput directory: ", output_dir, "\n"))
  
  if (!file.exists(landclim_path)) {
    cat(paste0("Could not find the LandClim executable!"))
  } else if (!file.exists(control_file)) {
    cat(paste0("The control file does not exist! Please check if the current path is all right: ", control_file, "\n"))
  } else {
    oldwd <- getwd()
    dir.create(output_dir, showWarnings = FALSE)
    file.remove(list.files("Output", full=TRUE))
    setwd(input_dir)
    system2(landclim_path, control_file)
    setwd(oldwd)
  }
}

set_landclim_path <- function(landclim_path){
  if(file.exists(landclim_path)) {
    assign("landclim_path", landclim_path, envir=globalenv())
  } else {
    print("Invalid path to LandClim file.")
  }
}

clean_output_ubuntu <- function(){
  fis <- list.files()
  file.copy(fis[grep("Output", fis)], paste("Output/",fis[grep("Output", fis)], sep="")) # Why is the `Output` directory copied into itself?
  file.remove(fis[grep(c("\\.(csv|txt)$"), fis)])
}

run_landclim <- function(control_file = "control.xml") {
  cat("Works only for Ubuntu until now!\n")
  cat("Working directory must be at the site level.\n")
  lc_path <- Sys.which("landclim")["landclim"]
  input_dir <- paste(getwd(), "Input", sep="/")
  output_dir <- paste(getwd(), "Output", sep="/")
  control_file <- paste(input_dir, "control.xml", sep="/")
  cat(paste0("Landclim executable: ", lc_path, "\nInput directory: ", input_dir, "\nOutput directory: ", output_dir, "\n"))

  if (!file.exists(lc_path)) {
    cat(paste0("Could not find the LandClim executable! Make sure you can start LandClim in any working directory only by typing 'landclim' (see manual for details).\n"))
  } else if (!file.exists(control_file)) {
    cat(paste0("The control file does not exist! It must be located at the following path: ", control_file, "\n"))
  } else {
    oldwd <- getwd()
    dir.create(output_dir, showWarnings = FALSE)
    file.remove(list.files("Output", full=TRUE))
    setwd(input_dir)
    system2(lc_path, control_file)
    setwd(oldwd)
    clean_output_ubuntu()
  }
}



clean_output_ubuntu <- function(){
  fis <- list.files()
  file.copy(fis[grep("Output", fis)], paste("Output/",fis[grep("Output", fis)], sep="")) # Why is the `Output` directory copied into itself?
  file.remove(fis[grep(c("\\.(csv|txt)$"), fis)])
}

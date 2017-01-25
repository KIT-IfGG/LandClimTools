
  landclim_colors <- function(n=10) {
    landclim_colors <- c("dodgerblue2","#E31A1C", # red
           "green4",
           "#6A3D9A", "purple",# purple
           "#FF7F00", "orange3",# orange
           "black","gold1",
           "skyblue2","#FB9A99", # lt pink
           "pink",
           "palegreen2",
           "#CAB2D6", # lt purple
           "#FDBF6F", # lt orange
           "gray70", "grey40", "khaki2",
           "khaki4",
           "maroon","orchid1","deeppink1","blue1","steelblue4",
           "darkturquoise","green1","yellow4","yellow3",
           "darkorange4","brown")
   if (n >= length(landclim_colors)) {
     print(paste("n larger than number of available colors. Max n is ", length(landclim_colors), ".",sep=""))
   } else {
   landclim_colors[1:n]
   }
}

  
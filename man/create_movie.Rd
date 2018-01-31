\name{create_movie}
\alias{create_movie}

\title{Create succession movie}
\description{
Creates a sequence of figures of the forest state from above displayed as a "movie".
}
\usage{
  create_movie(files, decades, a, species,  scol, plotlegend=TRUE, aspect=1, silent=FALSE)
}
\arguments{
  \item{files}{
    LandClim decade output "fullOut10" etc. One file per time step for the movie.
  }
  \item{decades}{
    Vector, Decade corresponding to files, used as labels.
  }
  \item{a}{
    LandClim cell size.
  }
  \item{species}{
    Simulated species, names used as column names in the LandClim output.
  }
  \item{scol}{
    Color vector. See also ?landclim_colors.
  }
  \item{plotlegend}{
    TRUE or FALSE
  }
  \item{aspect}{
    Function argument asp in plot().
  }
  \item{silent}{
    TRUE/FALSE
  }
}

\examples{

## full <- read.csv(file = "fullOut_30.csv", strip.white = TRUE)
## species <- unique(full$species)

## fi <- data.frame(file = c("fullOut_30.csv", "fullOut_31.csv","fullOut_32.csv","fullOut_33.csv","fullOut_34.csv"), decade = c(30, 31, 32, 33,34))

## pdf(paste(path,"movie.pdf", sep=""))
## par(mfrow=c(2,2))
## create_movie(files=fi$file, decades=fi$decade, a=25, species=species,  scol=rainbow(length(species)), plotlegend=TRUE, aspect=1)
## dev.off()

## library(animation)
## saveGIF(create_movie(files = fi$file, decades = fi$decade, a = 25, species = species,  scol = rainbow(length(species)), plotlegend=TRUE, aspect=1), movie.name = "LandClimForest.gif")

}

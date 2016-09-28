\name{change_climate}
\alias{change_climate}

\title{Simulate climate change.}
\description{
  The funcion reads a LandClim climate file and adds a constant to temperature and precipitation.
}
\usage{
  change_climate(inputPath, outputPath, dt = 0, dn = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputPath}{
    Path to LandClim climate file.
  }
  \item{outputPath}{
    Path to LandClim climate changed file.
  }
  \item{dt}{
    Temperature in K to add. Also negative values allowed.
  }
  \item{dn}{
    Precipitation to add in mm.  Also negative values allowed.
  }
}
\author{Klara Dolos}

\examples{
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
}

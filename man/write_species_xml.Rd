\name{write_species_xml}
\alias{write_species_xml}

\title{Write LandClim species parameter file}
\description{Writes the given data into the given file in the XML-format for species of LandClim}
\usage{
  write_species_xml(x, file)
}

\arguments{
  \item{x}{
    Data frame with species parameters and tag names.
  }
  \item{file}{
    Output file name including full path.
  }
}

\seealso{
read_species_xml
}

\examples{
  ### Read example species file
  species_parameters  <- read_species_xml(system.file("species.xml", package = "LandClimTools"))
  
  ### Write example species file
  write_species_xml(x = species_parameters, file = "piceabie.xml")
}

\name{read_species_xml}
\alias{read_species_xml}

\title{
Read species XML
}
\description{
Read LandClim species parameters from input file (species.xml).
}
\usage{
read_species_xml(file)
}

\arguments{
  \item{file}{ 

File name
}
}
\details{

}
\value{

Data frame with one row for each species.
}
\references{

}
\author{
Björn Reineking
}

\seealso{
write_species_xml
}
\examples{
  ### Read example species file
  species_parameters  <- read_species_xml(system.file("species.xml", package = "LandClimTools"))
  
  ### Write example species file
  write_species_xml(x = species_parameters, file = "piceabi.xml")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{profound_climate_to_landclim}
\alias{profound_climate_to_landclim}

\title{
Profound climate to LandClim 
}
\description{
Complile climate data from the Profound database to the format needed as input file for LandClim
}
\usage{
profound_climate_to_landclim(climate, header, file = "climate.txt")
}

\arguments{
  \item{climate}{
Dataframe as provided by the PROFOUND database, e.g. ProfoundData::getData("Soro", "CLIMATE").
}
  \item{header}{
Header information as provided by readLines("landclim_climate.txt"). I.e. read in an existing LandClim climate input file and provide it as function argument.
}
  \item{file}{
Output filename
}
}

\value{
Written text file with landclim climate.
}
\references{
The Profound Database is available here: https://github.com/COST-FP1304-PROFOUND/TG2 (NOT YET, Feb 2018, contact Christopher Reyer, PIK).

http://cost-profound.eu/site

Examples for simulations for cost sites (Geoecology Master course 2017/18, Karlsruhe)
https://github.com/KIT-IfGG/climate_impact_research_2

}
\author{
Klara Dolos
}


\examples{
 ### See https://github.com/KIT-IfGG/climate_impact_research_2 for example code

}

\name{profound_climate_to_landclim}
\alias{profound_climate_to_landclim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Complile climate from Profound database 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
profound_climate_to_landclim(climate, header, file = "climate.txt")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{climate}{
Dataframe as provided by the PROFOUND database, e.g. ProfoundData::getData("Soro", "CLIMATE").
}
  \item{header}{
Header information as provided by readLines("landclim_climate.txt").
}
  \item{file}{
Output filename.
}
}
\details{

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
 ### Todo

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

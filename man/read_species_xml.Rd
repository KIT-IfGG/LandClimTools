\name{read_species_xml}
\alias{read_species_xml}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Read LandClim species parameters.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Read LandClim species parameters from input file.(species.xml).
}
\usage{
read_species_xml(file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{ 
%%     ~~Describe \code{file} here~~
File name.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Data frame with one row for each species.
}
\references{
%% ~put references to the literature/web site here ~
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

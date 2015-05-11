\name{biomass_to_dbh}
\alias{biomass_to_dbh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
LandClim allomentry for biomass to DBH conversion.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
biomass_to_dbh(biomass, leafHabit, allometry = "SCHUMACHER")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{biomass}{
Numeric vector containing biomass in tonns.
}
  \item{leafHabit}{
LandClim leaf habit, e.g. "EVERGREEN", "BROADLEAFEVERGREEN", "DECIDUOUS".
}
  \item{allometry}{
Type of LandClim allometry to use, e.g. "SCHUMACHER", "POWER".
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
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
biomass_to_dbh(biomass=3, leafHabit="DECIDUOUS", allometry = "SCHUMACHER")
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

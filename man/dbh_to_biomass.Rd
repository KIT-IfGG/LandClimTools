\name{dbh_to_biomass}
\alias{dbh_to_biomass}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
LandClim allomentry for biomass -  DBH conversion.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dbh_to_biomass(dbh, leafHabit, allometry = "SCHUMACHER")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dbh}{
%%     ~~Describe \code{dbh} here~~
}
  \item{leafHabit}{
LandClim leaf habit, e.g. c("EVERGREEN", "BROADLEAFEVERGREEN", "DECIDUOUS")
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

}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{biomass_to_dbh}}

}
\examples{

dat <- rnorm(10, 20, 5)
dbh_to_biomass(dbh=dat, leafHabit="DECIDUOUS", allometry = "SCHUMACHER")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

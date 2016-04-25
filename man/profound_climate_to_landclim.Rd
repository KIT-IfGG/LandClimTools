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
Klara Dolos
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (climate, header, file = "climate.txt") 
{
    temp <- aggregate(climate$tmean_degC, by = list(mo = climate$mo, 
        year = climate$year), mean)
    precip <- aggregate(climate$p_mm, by = list(mo = climate$mo, 
        year = climate$year), sum)
    temp <- round(reshape(temp[, c("year", "mo", "x")], timevar = "mo", 
        idvar = "year", direction = "wide"), 2)
    precip <- round(reshape(precip[, c("year", "mo", "x")], timevar = "mo", 
        idvar = "year", direction = "wide"), 0)
    climate <- cbind(temp[, grep("x", names(temp))], precip[, 
        grep("x", names(temp))])
    writeLines(header, con = file)
    write.table(climate, file = file, append = TRUE, col.names = FALSE)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

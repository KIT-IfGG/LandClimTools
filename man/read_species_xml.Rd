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
%%  ~~who you are~~
Björn Reineking
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
write_species_xml
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (file) 
{
    require(XML)
    doc <- xmlTreeParse(file)
    daten <- t(xmlSApply(xmlRoot(doc), function(x) xmlSApply(x, 
        xmlValue)))
    rownames(daten) <- NULL
    daten <- data.frame(daten)
    tmpfile <- file()
    write.table(daten, tmpfile)
    daten <- read.table(tmpfile)
    daten
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

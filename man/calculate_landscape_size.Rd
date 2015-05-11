\name{calculate_landscape_size}
\alias{calculate_landscape_size}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~ 
Calculate landscape size
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Function ot caculate the landscape size using the LandClim input digital elevation ascii-file (e.g. dem.asc).
}
\usage{
calculate_landscape_size(dem)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dem}{
%%     ~~Describe \code{dem} here~~
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (dem) 
{
    prod(dim(dem)) * prod(res(dem))/(100 * 100)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

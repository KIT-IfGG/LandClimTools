\name{run_landclim}
\alias{run_landclim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Run LandClim on Ubuntu
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
run_landclim(control_file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{control_file}{Name of the control file (e.g. control.xml).
%%     ~~Describe \code{x} here~~
}
}
\details{
The working directory needs to be set to the level above "Input" and "Output" folders of LandClim, thus at the "site" level. Until now, the folder names are hard-coded, thus need to be "Input" and "Output", not e.g. input! Feel free to fix this!
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

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

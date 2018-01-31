\name{landclim_colors}
\alias{landclim_colors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
LandClim color palette
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
landclim_colors(n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Number of different tree species = number of colors.
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
Please feel free to enhance this palette!
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
n <- 10
pie(rep(1, n), col=landclim_colors(n))

### Get number of available colors.
landclim_colors(99)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

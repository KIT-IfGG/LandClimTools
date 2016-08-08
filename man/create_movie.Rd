\name{create_movie}
\alias{create_movie}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create succession movie
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
create_movie(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
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

## full <- read.csv(file="fullOut_30.csv", strip.white=TRUE)
## species <- unique(full$species)

## fi <- data.frame(file=c("fullOut_30.csv", "fullOut_31.csv","fullOut_32.csv","fullOut_33.csv","fullOut_34.csv"), decade=c(30, 31,32,33,34))

## pdf(paste(path,"movie.pdf", sep=""))
## par(mfrow=c(2,2))
## create_movie(files=fi$file, decades=fi$decade, a=25, species=species,  scol=rainbow(length(species)), plotlegend=TRUE, aspect=1)
## dev.off()

## library(animation)
## saveGIF(create_movie(files=fi$file, decades=fi$decade, a=25, species=species,  scol=rainbow(length(species)), plotlegend=TRUE, aspect=1), movie.name="LandClimForest.gif")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

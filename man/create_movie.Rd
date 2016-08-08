\name{create_movie}
\alias{create_movie}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create succession movie}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
  create_movie(files, decades, a, species,  scol, plotlegend=TRUE, aspect=1, silent=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{files}{
    % TODO: document this parameter
  }
  \item{decades}{
    % TODO: document this parameter
  }
  \item{a}{
    % TODO: document this parameter
  }
  \item{species}{
    % TODO: document this parameter
  }
  \item{scol}{
    % TODO: document this parameter
  }
  \item{plotlegend}{
    % TODO: document this parameter
  }
  \item{aspect}{
    % TODO: document this parameter
  }
  \item{silent}{
    % TODO: document this parameter
  }
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

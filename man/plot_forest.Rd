\name{plot_forest}
\alias{plot_forest}

\title{Plot LandClim forest}
\description{Plots a forest based on tree data from LandClim}

\usage{
  plot_forest(trees, species = unique(trees$species), scol = rainbow(length(species)), plotlegend = TRUE, a = 25, aspect = 1, cex = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trees}{Object created by function tree_coordinates.}
  \item{species}{Species names, character vector.}
  \item{scol}{Color vector.}
  \item{plotlegend}{TRUE or FALSE indicating whether to plot a legend.}
  \item{a}{Resolution of LandClim maps used in the simulation.}
  \item{aspect}{Argument asp in function plot.}
  \item{cex}{Argument cex in function plot.}
}
\author{
  Klara Dolos
}

\seealso{
  \code{\link{tree_coordinates}}
}
\examples{

  trees <- tree_coordinates(file=system.file("fullOut_50.csv", package = "LandClimTools"), a=25)

  range(trees$col)
  stand <- trees[trees$row > 20 & trees$row <=40,]
  stand$row <- stand$row - min(stand$row)

  stand <- trees[trees$col > 20 & trees$col <=40,]
  stand$col <- stand$col - min(stand$col)

  x11(width=7, height=7)
  par(mar=c(2,2,1,1))
  plot_forest(trees=stand, species=unique(stand$species),  scol=rainbow(length(unique(stand$species))), plotlegend=TRUE, aspect=1, cex=sqrt(stand$biomass)/2)

}
